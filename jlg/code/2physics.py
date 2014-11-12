#!/usr/bin/env python3

import numpy as np
import scipy as sp
import scipy.sparse as spar
import scipy.optimize as opt
class boundary_condition:
		def __init__(self, BC):
			self.ltype = BC[0]
			self.rtype = BC[2]
			if callable(BC[1]):
				self.left = BC[1]
			else:
				self.left = lambda t: BC[1]
			if callable(BC[3]):
				self.right = BC[3]
			else:
				self.right = lambda t: BC[3]
class diffusion_system:
	def __init__(self, k, Dk, Q, BC, sig=0):
		if callable(k):
			self.k = k
		else:
			self.k = lambda x,t,uph,uth: k
		if callable(Dk):
			self.Dk = Dk
		else:
			self.Dk = lambda x,t,uph,uth: Dk
		if callable(Q):
			self.Q = Q
		else:
			self.Q = lambda x,t,uph,uth: Q
		if callable(sig):
			self.sig = sig
		else:
			self.sig = lambda x,t,uph,uth: sig
		self.BC = boundary_condition(BC)
class boundary:
	'''
	Boundary conditions container
	'''
	def __init__(self, l,lv,r,rv):
		self.left = l
		self.right = r
		self.lval = lv
		self.rval = rv
	def astuple(self):
		return (self.left, self.lval, self.right, self.rval)
class FEM_SYS:
	def __init__(self, N, T, photon_system, thermal_system, order=2):
		# Discretation
		self.N = N # Number of spatial steps
		self.T = T # Number of temporal steps
		self.x = np.linspace(0,1,N+1) # X values of nodes (N+1)
		
		# Diffusion parameters
		self.ph = photon_system  # Photon equation parameters
		self.th = thermal_system # Thermal equation parameters
		
		# Quadrature points used for integration
		# Weights of quadrature
		self.qorder = order
		[self.xq, self.wq] = np.polynomial.legendre.leggauss(self.qorder)
		# Basis functions at quadrature points
		self.b = np.array([(1-self.xq)/2, (1+self.xq)/2])
		# Derivative of the basis functions at quadrature points
		self.dbdx = np.array([[-.5, .5],[-.5, .5]])
	def forward_solve(self, u0):
		u_for = np.empty((self.T,(self.N+1)**2))
		u_for[0,:] = u0
		for t in range(1, self.T+1):
			u_for[t,:] = forward_time_step(u_for[t-1,:], (t-1)/T)
		return u_for
	def forward_time_step(self, u_cur, t_cur):
		"""
		Takes a single (crank nicolson) time step forward
		Returns the flux at the next time step.
		"""
	def newton_jacobian(self, u_k, u_cur, t_cur, F0):
		"""
		Calculate the jacobian for the FEM_STATE
		Takes:
		u_k, current guess
		u_cur, flux at time t_cur
		t_cur, current time
		R0, initial residual normal
		"""
		JAC = spar.csr_matrix(((self.N+1)*2, (self.N+1)*2))
	def newton_residual_half(self, u_cur, t, d_sys, part):
		offset = (self.N+1)*part
		F = np.zeros(self.N+1)
		for k in range(self.N):
			x0 = self.x[k]
			x1 = self.x[k+1]
			jacob = (x1-x0)/2
			# Quadrature points mapped to the element
			x = (x1+x0)/2 + self.xq*jacob
			# u and du/dx values at quadrature points
			ph_uq = np.dot(self.b,    u_cur[k:k+self.qorder])
			th_uq = np.dot(self.b,    u_cur[k+self.N+1:k+self.N+1+self.qorder])
			
			uq    = np.dot(self.b,    u_cur[k + offset:k + offset + self.qorder])
			dudxq = np.dot(self.dbdx, u_cur[k + offset:k + offset + self.qorder])
			# k(x,t,u) and dk(x,t,u)/du values at quadrature points
			kq  = d_sys.k (x, t, ph_uq, th_uq)/jacob
			Dkq = d_sys.Dk(x, t, ph_uq, th_uq)
			# sigma(x,t,u) at quadrature points
			sigq = d_sys.sig(x, t, ph_uq, th_uq)
			# Q(x,t,u) at the quadrature points
			Qq = d_sys.Q(x, t, ph_uq, th_uq)
			local = np.zeros(self.qorder)
			for j in range(self.qorder):
				local[j] = sum(uq*sigq*self.wq*self.b[:,j]) + np.dot(kq*self.wq*self.dbdx[:,j], dudxq) - np.dot(Qq*self.wq, self.b[:,j])
			F[k:k+self.qorder] += local
		# Apply Boundary Conditions
		if d_sys.BC.ltype == 0: # Neumann
			F[0] -= d_sys.BC.left(t)
		elif d_sys.BC.ltype  == 1: # Robin
			F[0] += 1/2 * u_cur[offset] - 2*d_sys.BC.left(t)
		elif d_sys.BC.ltype  == 2: # Dirichlet
			F[0] = u_cur[offset] - d_sys.BC.left(t)
		if d_sys.BC.rtype == 0: # Neumann
			F[-1] -= d_sys.BC.right(t)
		elif d_sys.BC.rtype == 1: # Robin
			F[-1] += 1/2 * u_cur[N+offset] - 2*d_sys.BC.right(t)
		elif d_sys.BC.rtype == 2: # Dirichlet
			F[-1] = u_cur[N+offset] - d_sys.BC.right(t)
		return F
	def newton_residual(self, u_cur, t, F0=1):
		"""
		Calculate the residual divided by the initial residual
		Inital defaults to "1" so it returns the undivided residual if not provided.
		Takes:
		u_cur, current guess
		t, current time
		F0, normal of initial residual 
		"""
		F = np.zeros((self.N+1)*2)
		F[0:self.N+1] = self.newton_residual_half(u_cur, t, self.ph, 0)
		F[self.N+1: ] = self.newton_residual_half(u_cur, t, self.th, 1)
		return sum(F**2)/F0
	def newton_solve(self, u0, t):
		"""
		Use newton's method to solve at time=1, using the initial guess u0
		"""
		if u0 == None:
			u0 = np.zeros((N+1)*2)
		F0 = self.newton_residual(u0, t)
		if F0 == 0:
			return u0
		return opt.minimize(
				fun = self.newton_residual,
				x0 = u0,
				args = (t, F0),
				#jac = self.newton_jacobian,
				tol = 1E-12)

N = 2
T = 1
ph = diffusion_system(0,0,1,(2,1,2,1),1)
th = diffusion_system(0,0,0,(0,0,0,0))
constant = FEM_SYS(N, T, ph, th)
u_0 = np.asarray([0,0,0,0,0,0])
u_sol = np.asarray([1,1,1,0,0,0])
print(constant.newton_residual(u_0,0))
print(constant.newton_residual(u_sol,0))
print(constant.newton_solve(None,0))

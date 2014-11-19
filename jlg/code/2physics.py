#!/usr/bin/env python3

import numpy as np
import scipy as sp
import scipy.sparse as spar
import scipy.optimize as opt
from numpy.linalg import norm

np.set_printoptions(linewidth=200)

class boundary_condition:
		def __init__(self, BC):
			self.ltype = BC[0]
			self.rtype = BC[2]
			self.left  = BC[1] if callable(BC[1]) else lambda t: BC[1]
			self.right = BC[3] if callable(BC[3]) else lambda t: BC[3]
class diffusion_system:
	def __init__(self, N, k, Dk, Q, BC, offset, sig=0, order=2):
		self.N = N
		self.x = np.linspace(0,1,N+1) # X values of nodes (N+1)
		self.offset = offset
		
		self.k   = k   if callable(k)   else lambda x,t,uph,uth: k
		self.Dk  = Dk  if callable(Dk)  else lambda x,t,uph,uth: Dk
		self.Q   = Q   if callable(Q)   else lambda x,t,uph,uth: Q
		self.sig = sig if callable(sig) else lambda x,t,uph,uth: sig
		
		# Quadrature points used for integration
		# Weights of quadrature
		self.qorder = order
		[self.xq, self.wq] = np.polynomial.legendre.leggauss(self.qorder)
		# Basis functions at quadrature points
		self.b = np.array([(1-self.xq)/2, (1+self.xq)/2])
		# Derivative of the basis functions at quadrature points
		self.dbdx = np.array([[-.5, .5],[-.5, .5]])
		# Boundary Conditions
		self.BC = boundary_condition(BC)
		# Generate and save the mass matrix
		self.mass_matrix = self.generate_mass_matrix(N)
	def generate_mass_matrix(self, N):
		"""
		Mass matrix with a flux = u and at time = t
		"""
		h = 1/(N)
		# Create the mass matrix
		d_main = 4/3 * np.ones(N+1) # Diagonal Entries
		d_off = 1/3 * np.ones(N+1) # Off diagonals
		M = spar.spdiags([d_off, d_main, d_off], [-1,0,1], N+1, N+1, format='csc')
		M[0,0] /= 2
		M[-1,-1] /= 2
		return M*h/2
	def stiffness_matrix(self, u0, t):
		"""
		Stiffness matrix with a flux = u and at time = t
		"""
		A_nobc = np.zeros((self.N+1,self.N+1))
		f = np.zeros(self.qorder)
		for i in range(self.qorder):
			f[i] = sum(self.wq*self.b[:,i])
		for k in range(self.N):
			# u and du/dx values at quadrature points
			ph_uq = np.dot(self.b,    u0[k:k+self.qorder])
			th_uq = np.dot(self.b,    u0[k+self.N+1:k+self.N+1+self.qorder])

			x0 = self.x[k]
			x1 = self.x[k+1]
			jacob = (x1-x0)/2
			# Quadrature points mapped to the element
			x = (x1+x0)/2 + self.xq*jacob
			# u and du/dx values at quadrature points
			uq    = np.dot(self.b,    u0[k + self.offset*(N+1):k + self.offset*(N+1) + self.qorder])
			dudxq = np.dot(self.dbdx, u0[k + self.offset*(N+1):k + self.offset*(N+1) + self.qorder])
			# D(x,u) values at quadrature points
			Dq = self.k(uq, x, ph_uq, th_uq)/jacob
			dDduq = self.Dk(uq, x, ph_uq, th_uq)
			# Sigma values at quadrature points
			Sq = self.sig(uq, x, ph_uq, th_uq)*jacob
			Alet = np.zeros((self.qorder,self.qorder))
			for i in range(self.qorder):
				for j in range(self.qorder):
					Alet[i,j] = sum(self.wq*self.dbdx[:,i]*Dq*self.dbdx[:,j]) + sum(self.wq*Sq*self.b[:,i]*self.b[:,j])
			A_nobc[k:k+self.qorder,k:k+self.qorder] += Alet
		return A_nobc
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
		self.dt = 1/T # Time step width
		self.x = np.linspace(0,1,N+1) # X values of nodes (N+1)
		
		# Diffusion systems, creates a constant zero system if it is None
		self.ph =  photon_system if  photon_system is not None else diffusion_system(N, 0, 0, 0, (2,lambda t: 0,2,lambda t: 0), 0)
		self.th = thermal_system if thermal_system is not None else diffusion_system(N, 0, 0, 0, (2,lambda t: 0,2,lambda t: 0), 0)
		
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
		raise NotImplementedError('Function not implemented')
	def newton_jacobian(self, u_k, u_cur, t_cur):
		"""
		Calculate the jacobian for the FEM_STATE
		Takes:
		u_k, current guess
		u_cur, flux at time t_cur
		t_cur, current time
		R0, initial residual normal
		"""
		JAC = spar.csr_matrix(((self.N+1)*2, (self.N+1)*2))
		# Mass + Stiffness for photon part (Top Left)
		
		# Mass + Stiffness for thermal part (Bottom Right)
		
		# Mass (Top Right)
		
		# (Bottom Left)
		
		raise NotImplementedError('Function not implemented')
		
	def newton_residual_half(self, u_cur, Auold, t, d_sys, part):
		offset = (self.N+1)*part
		F = np.zeros(self.N+1)
		
		# (M+A(k-1)*tau/2)*u(k-1)
		# Previous time step for time dependence.
		ph_A = self.ph.mass_matrix + self.dt/2 * self.ph.stiffness_matrix(u_cur,t)
		th_A = self.th.mass_matrix + self.dt/2 * self.th.stiffness_matrix(u_cur,t)
		Aucur = np.empty((N+1)*2)
		Aucur[:len(u_cur)//2] = np.dot(ph_A, u_cur[:len(u_cur)//2])
		Aucur[len(u_cur)//2:] = np.dot(th_A, u_cur[len(u_cur)//2:])
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
			
			Auoldq = np.dot(self.b, Auold[k + offset:k + offset + self.qorder])
			Aucurq = np.dot(self.b, Aucur[k + offset:k + offset + self.qorder])
			local = np.zeros(self.qorder)
			for j in range(self.qorder):
				local[j]  = sum(uq*sigq*self.wq*self.b[:,j])
				local[j] += np.dot(kq*self.wq*self.dbdx[:,j], dudxq)
				local[j] -= np.dot(Qq*self.wq, self.b[:,j])
				#local[j] += np.dot(Aucurq*self.wq, self.b[:,j]) # Time Term  ( M+A(k)  *tau/2)u(k)
				#local[j] += np.dot(Auoldq*self.wq, self.b[:,j]) # Time Term -(-M+A(k-1)*tau/2)u(k-1)
#				print('Local',k,j,local[j])
#				print(sum(uq*sigq*self.wq*self.b[:,j]))
#				print(np.dot(kq*self.wq*self.dbdx[:,j], dudxq))
#				print(-np.dot(Qq*self.wq, self.b[:,j]))
#				print(sum(self.wq*np.dot(self.b, Aucur[k + offset:k + offset + self.qorder])*self.b[:,j]))
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
	def newton_residual(self, u_cur, u_old, t):
		"""
		Calculate the residual
		Takes:
		u_cur, current guess
		u_old, solution of previous time step
		t, current time
		F0, normal of initial residual 
		"""
		F = np.zeros((self.N+1)*2)
		# (-M+A(k-1)*tau/2)*u(k-1)
		# Previous time step for time dependence.
		ph_A = -self.ph.mass_matrix + self.dt/2 * self.ph.stiffness_matrix(u_old,t - self.dt)
		th_A = -self.th.mass_matrix + self.dt/2 * self.th.stiffness_matrix(u_old,t - self.dt)
		Auold = np.empty((N+1)*2)
		Auold[:len(u_old)//2] = np.dot(ph_A, u_old[:len(u_old)//2])
		Auold[len(u_old)//2:] = np.dot(th_A, u_old[len(u_old)//2:])
		
		# Each half of the residual
		F[0:self.N+1] = self.newton_residual_half(u_cur, Auold, t, self.ph, 0)
		F[self.N+1: ] = self.newton_residual_half(u_cur, Auold, t, self.th, 1)
		print(u_old)
		print('\t',F)
		return F
	def newton_residual_norm(self, u_cur, u_old, t, normed=1):
		if u_old is None:
			u_old = u_cur
		F = self.newton_residual(u_cur, u_old, t)
		R = norm(F)/normed
		return R
	def newton_solve(self, u0, u_old, t):
		"""
		Use newton's method to solve at time=1, using the initial guess u0
		"""
		if u0 is None:
			u0 = np.zeros((N+1)*2)
		res0 = self.newton_residual_norm(u0,u_old,t)
		return opt.minimize(
				fun = self.newton_residual_norm,
				x0 = u0,
				args = (u_old,t,res0),
				tol=1E-16
				#jac = self.newton_jacobian,
				)
def test_problem(ph, T, ID):
	u = FEM_SYS(ph.N, T, ph, None)
	u_sol = np.asarray([1,1,1,0,0,0])
	print(ID, u.newton_solve(None,None,0))
	print('Sol', u.newton_residual_norm(u_sol,u_sol,0))
T = 1

N = 2
kph  = lambda x, t, ph_uq, th_uq: 1
Dkph = lambda x, t, ph_uq, th_uq: 0
Qph  = lambda x, t, ph_uq, th_uq: 0
BCph = (2,lambda t: 1,2,lambda t: 1)
ph_1 = diffusion_system(N, kph, Dkph, Qph, BCph, offset=0)
test_problem(ph_1, T, 'u=1')

'''
N = 2
kph  = lambda x, t, ph_uq, th_uq: x
Dkph = lambda x, t, ph_uq, th_uq: 0
Qph  = lambda x, t, ph_uq, th_uq: 1
BCph = (2,lambda t: 0,2,lambda t: 1)
ph_x = diffusion_system(N, kph, Dkph, Qph, BCph, offset=0)
test_problem(ph_x, T, 'u=x')

N = 2
kph  = lambda x, t, ph_uq, th_uq: 1
Dkph = lambda x, t, ph_uq, th_uq: 0
Qph  = lambda x, t, ph_uq, th_uq: 2*x
BCph = (2,lambda t: 0,2,lambda t: 1)
ph_x2 = diffusion_system(N, kph, Dkph, Qph, BCph, offset=0)
test_problem(ph_x2, T, 'u=x²')
'''

#!/usr/bin/env python3

import numpy as np
import scipy as sp
import scipy.sparse as spar
import scipy.optimize as opt
from numpy.linalg import norm
import sympy

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
		return spar.csr_matrix(M*h/2)
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
	def residual(self, u_cur, u_prev, t_cur, dt):
		#M(new-old)/dt - (Anew + Aold)/2 - (Qnew + Qold)/2
		l = len(u_cur)
		t0 = t_cur - dt
		t1 = t_cur
		if self.offset==0:
			this1 =  u_cur[0:l//2] # This Physics
			this0 = u_prev[0:l//2] # This Physics
			othr1 = u_cur[l//2:]   # Other Physics
			othr0 = u_prev[l//2:]  # Other Physics
			u0ph = this0
			u1ph = this1
			u0th = othr0
			u1th = othr1
		else:
			this1 = u_cur[l//2:]   # This Physics
			this0 = u_prev[l//2:]  # This Physics
			othr1 =  u_cur[0:l//2] # Other Physics 
			othr0 = u_prev[0:l//2] # Other Physics
			u0ph = othr0
			u1ph = othr1
			u0th = this0
			u1th = this1
		#M(Δu)/dt
		Mudt = self.mass_matrix.dot(this1-this0)/dt
		F = np.zeros(self.N+1)
		
		for i in range(self.N):
			x0 = self.x[i]
			x1 = self.x[i+1]
			jacob = (x1-x0)/2
			# Quadrature points mapped to the element
			x = (x1+x0)/2 + self.xq*jacob
			# u values at quadrature points
			u0q = np.dot(self.b, this0[i:i+self.qorder]) # This Physics
			u1q = np.dot(self.b, this1[i:i+self.qorder]) # This Physics
			ph0q = np.dot(self.b, u0ph[i:i+self.qorder])
			ph1q = np.dot(self.b, u1ph[i:i+self.qorder])
			th0q = np.dot(self.b, u0th[i:i+self.qorder])
			th1q = np.dot(self.b, u1th[i:i+self.qorder])
			# du/dx values at quadrature points
			du0dxq = np.dot(self.dbdx, this0[i:i+self.qorder])
			du1dxq = np.dot(self.dbdx, this1[i:i+self.qorder])
			dudxq = (du0dxq + du1dxq)/2
			# D(x,u) values at quadrature points
			k0q =   self.k(x, t0, ph0q, th0q)/jacob
			k1q =   self.k(x, t1, ph1q, th1q)/jacob
			# Sigma(x,u) values at quadrature points
			S0q = self.sig(x, t0, ph0q, th0q)*jacob
			S1q = self.sig(x, t1, ph1q, th1q)*jacob
			# q(x) at the quadrature points
			Q0q =   self.Q(x, t0, ph0q, th0q)*jacob
			Q1q =   self.Q(x, t1, ph1q, th1q)*jacob
			# M(Δu)/dt at the quadrature points
			Mudtq = np.dot(self.b, Mudt[i:i+self.qorder])
			uq  = (u0q+u1q)/2
			kq  = (k0q + k1q)/2
			Suq = (S0q*u0q + S1q*u1q)/2
			Qq  = (Q0q + Q1q)/2
			local = np.zeros(self.qorder)
			for j in range(self.qorder):
				local[j]	= np.dot(Mudtq*self.wq, self.b[:,j])\
							+ np.dot(  Suq*self.wq, self.b[:,j])\
							+ np.dot(   kq*self.wq*self.dbdx[:,j], dudxq)\
							- np.dot(   Qq*self.wq, self.b[:,j])
			F[i:i+2] += local
		# Apply Boundary Conditions
		if self.BC.ltype == 0: # Neumann
			F[0] -= self.BC.left(t)
		elif self.BC.ltype  == 1: # Robin
			F[0] += 1/2 * this1[0] - 2*self.BC.left(t1)
		elif self.BC.ltype  == 2: # Dirichlet
			F[0] = this1[0] - self.BC.left(t1)
		if self.BC.rtype == 0: # Neumann
			F[-1] -= self.BC.right(t)
		elif self.BC.rtype == 1: # Robin
			F[-1] += 1/2 * this1[-1] - 2*self.BC.right(t1)
		elif self.BC.rtype == 2: # Dirichlet
			F[-1] = this1[-1] - self.BC.right(t1)
		return F

class FEM_SYS:
	def __init__(self, N, T, photon_system, thermal_system, order=2):
		# Discretation
		self.N = N # Number of spatial steps
		self.T = T # Number of temporal steps
		self.dt = 1/T # Time step width
		self.x = np.linspace(0,1,N+1) # X values of nodes (N+1)
		
		# Diffusion systems, creates a constant zero system if it is None
		self.ph =  photon_system if  photon_system is not None else diffusion_system(N, 1, 0, 0, (2,lambda t: 1,2,lambda t: 1), 0)
		self.th = thermal_system if thermal_system is not None else diffusion_system(N, 1, 0, 0, (2,lambda t: 1,2,lambda t: 1), 0)
		
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
	def newton_jacobian(self, u_cur, u_old, t):
		"""
		Calculate the jacobian for the FEM_STATE
		Takes:
		u_k, current guess
		u_cur, flux at time t_cur
		t_cur, current time
		R0, initial residual normal
		"""
		if u_old is None:
			u_old = u_cur
		JAC = spar.dok_matrix(((self.N+1)*2, (self.N+1)*2))
		# Mass + Stiffness for photon part (Top Left)
		
		# Mass + Stiffness for thermal part (Bottom Right)
		
		# Mass (Top Right)
		
		# (Bottom Left)
		R0 = self.newton_residual(u_cur, u_old, t)
		ep = 1E-8
		for i in range((self.N+1)*2):
			epsilon = np.zeros((self.N+1)*2)
			epsilon[i] = ep
			R1 = self.newton_residual(u_cur+epsilon, u_old, t)
			D = (R1-R0)/ep
			for j in range((self.N+1)*2):
				if D[j] != 0:
					JAC[j,i] = D[j] 
		return JAC
	def newton_residual_half(self, u_cur, u_old, t, d_sys, part):
		offset = (self.N+1)*part
		F = np.zeros(self.N+1)
		
		# M*u_cur
		# Previous time step for time dependence.
		cur_ph_A = self.ph.mass_matrix
		cur_th_A = self.th.mass_matrix
		cur_Mu = np.zeros((N+1)*2)
		cur_Mu[:len(u_cur)//2] = cur_ph_A.dot(u_cur[:len(u_cur)//2])
		cur_Mu[len(u_cur)//2:] = cur_th_A.dot(u_cur[len(u_cur)//2:])
		# M*u_old
		# Previous time step for time dependence.
		old_ph_A = self.ph.mass_matrix# + self.dt/2 * self.ph.stiffness_matrix(u_old,t - self.dt)
		old_th_A = self.th.mass_matrix# + self.dt/2 * self.th.stiffness_matrix(u_old,t - self.dt)
		old_Mu = np.zeros((N+1)*2)
		old_Mu[:len(u_old)//2] = old_ph_A.dot(u_old[:len(u_old)//2])
		old_Mu[len(u_old)//2:] = old_th_A.dot(u_old[len(u_old)//2:])
		for k in range(self.N):
			x0 = self.x[k]
			x1 = self.x[k+1]
			jacob = (x1-x0)/2
			# Quadrature points mapped to the element
			x = (x1+x0)/2 + self.xq*jacob
			# u and du/dx values at quadrature points
			cur_ph_uq = np.dot(self.b,    u_cur[k:k+self.qorder])
			cur_th_uq = np.dot(self.b,    u_cur[k+self.N+1:k+self.N+1+self.qorder])
			old_ph_uq = np.dot(self.b,    u_old[k:k+self.qorder])
			old_th_uq = np.dot(self.b,    u_old[k+self.N+1:k+self.N+1+self.qorder])
			
			cur_uq    = np.dot(self.b,    u_cur[k + offset:k + offset + self.qorder])
			cur_dudxq = np.dot(self.dbdx, u_cur[k + offset:k + offset + self.qorder])
			old_uq    = np.dot(self.b,    u_old[k + offset:k + offset + self.qorder])
			old_dudxq = np.dot(self.dbdx, u_old[k + offset:k + offset + self.qorder])
			# k(x,t,u) and dk(x,t,u)/du values at quadrature points
			cur_kq  = d_sys.k(x, t,   cur_ph_uq, cur_th_uq)/jacob
			old_kq  = d_sys.k(x, t-1, old_ph_uq, old_th_uq)/jacob
			cur_Dkq = d_sys.Dk(x, t,   cur_ph_uq, cur_th_uq)/jacob
			old_Dkq = d_sys.Dk(x, t-1, old_ph_uq, old_th_uq)/jacob
			# sigma(x,t,u) at quadrature points
			cur_sigq = d_sys.sig(x, t,   cur_ph_uq, cur_th_uq)*jacob
			old_sigq = d_sys.sig(x, t-1, old_ph_uq, old_th_uq)*jacob
			# Q(x,t,u) at the quadrature points
			cur_Qq = d_sys.Q(x, t  , cur_ph_uq, cur_th_uq)*jacob
			old_Qq = d_sys.Q(x, t-1, old_ph_uq, old_th_uq)*jacob
			
			# Mass Matrix Terms
			cur_Muq = np.dot(self.b, cur_Mu[k + offset:k + offset + self.qorder])
			old_Muq = np.dot(self.b, old_Mu[k + offset:k + offset + self.qorder])

			local = np.zeros(self.qorder)
			for j in range(self.qorder):
				cur_A_term	 = sum(cur_uq*cur_sigq*self.wq*self.b[:,j])
				cur_A_term	-= np.dot(cur_kq*self.wq*self.dbdx[:,j], cur_dudxq)
				
				old_A_term	 = sum(old_uq*old_sigq*self.wq*self.b[:,j])
				old_A_term	-= np.dot(old_kq*self.wq*self.dbdx[:,j], old_dudxq)
				
				cur_Q_term	 = np.dot(cur_Qq*self.wq, self.b[:,j])
				old_Q_term	 = np.dot(old_Qq*self.wq, self.b[:,j])
				local[j] = -(cur_A_term+old_A_term)*self.dt/2 - (cur_Q_term+old_Q_term)*self.dt/2
				local[j] += np.dot(cur_Muq*self.wq,self.b[:,j]) # Time Term  M unew
				local[j] -= np.dot(old_Muq*self.wq,self.b[:,j]) # Time Term -M uold
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
		if u_old is None:
			u_old = np.copy(u_cur)
		F = np.zeros((self.N+1)*2)
		# Each half of the residual
		F[0:self.N+1] = self.ph.residual(u_cur, u_old, t, self.dt)#self.newton_residual_half(u_cur, u_old, t, self.ph, 0)
		F[self.N+1: ] = self.th.residual(u_cur, u_old, t, self.dt)#self.newton_residual_half(u_cur, u_old, t, self.th, 1)
		return F
	def newton_residual_norm(self, u_cur, u_old, t, normed=1):
		if u_old is None:
			u_old = np.copy(u_cur)
		F = self.newton_residual(u_cur, u_old, t)
		R = sum(F**2)
		return R
	def newton_solve(self, u0, u_old, t):
		"""
		Use newton's method to solve at time=1, using the initial guess u0
		"""
		if u0 is None:
			u0 = np.zeros((self.N+1)*2)
		res0 = self.newton_residual_norm(u0,u_old,t)
		if res0 <= 1E-16:
			return u0
		ans = opt.leastsq(
				func = self.newton_residual,
				x0 = u0,
				args = (u_old,t),
				ftol=1E-15
				#jac = self.newton_jacobian,
				)
		return ans[0]

def manufacture_solution(N, k_str, s_str, uph_str, uth_str, is_ph):
	"""
	Constructs a diffusion_system using sympy to create the source term
	Takes k, Σ, and the solution as strings
	"""
	from sympy.parsing.sympy_parser import parse_expr
	from sympy import diff
	# Declare Sympy symbols
	t, x, up, ut = sympy.symbols('t x uph uth')
	uph = parse_expr(uph_str)
	uth = parse_expr(uth_str)
	# Create Sympy expressions
	# Create  k(x,t,uph,uth) function
	k = parse_expr(k_str)
	# Create Dk(x,t,uph,uth) function
	Dk = diff(k,up) if is_ph else diff(k,ut)
	# Create  Σ(x,t,uph,uth) function
	Sig = parse_expr(s_str)
	# Create  Q(x,t,uph,uth) function
	u = uph if is_ph else uth
	Q = diff(u,t) - k*diff(diff(u,x),x)-Dk*diff(u,x)**2 + Sig*u
	# Create BC by solving u for x=0,1
	bc_left  = sympy.lambdify(t,u.subs(x,0),"numpy")
	bc_right = sympy.lambdify(t,u.subs(x,1),"numpy")
	
	# Lambdify for construction of diffusion_system
	kfun     = sympy.lambdify((x,t,up,ut),k,"numpy")
	Dkfun    = sympy.lambdify((x,t,up,ut),Dk,"numpy")
	Qfun     = sympy.lambdify((x,t,up,ut),Q,"numpy")
	sigfun   = sympy.lambdify((x,t,up,ut),Sig,"numpy")
	BC       = (2,bc_left,2,bc_right)
	
	# Create diffusion system
	offset = 0 if is_ph else 1
	d_sys = diffusion_system(N, kfun, Dkfun, Qfun, BC, offset, sigfun)
	return d_sys
def innerproduct(u, res, N, T):
	tau = 1/T
	[xq, wq] = np.polynomial.legendre.leggauss(2)
	b = np.array([(1-xq)/2, (1+xq)/2])
	qoif_t = np.zeros(T+1)
	for t in range(T+1):
		ut = u#[t,:]
		for i in range(N):
			x0 = i/N
			x1 = (i+1)/N
			jacob = (x1-x0)/2
			x = (x1+x0)/2 + xq*jacob
			qoif_t[t] += jacob * sum(wq*np.dot(b, ut[i:i+2])*res(x,t))
	return tau/2*sum(qoif_t[1::] + qoif_t[0:-1])
def test_problem(R, sol, ph, th, T, ID):
	u = FEM_SYS(ph.N, T, ph, th)
	res = u.newton_solve(None,None,1)
	qoia = innerproduct(res[0:ph.N+1], R, ph.N, T)
	print("%8s: %8f %8f"%(ID,qoia, abs(qoia-sol)),end="\t")
	print(res)

R = lambda x,t: 1
T = 1
N = 2
print("%8s  %8s %8s\t%8s"%("Problem","QOIa","Error","Vector  "))
uph = "0"
uth = "1"
ph_1 = manufacture_solution(N, "1", "0", uph, uth, True)
th_1 = manufacture_solution(N, "1", "0", uph, uth, False)
test_problem(R, 0, ph_1, th_1, T, "u="+uph)
uph = "1"
uth = "1"
ph_1 = manufacture_solution(N, "1", "0", uph, uth, True)
th_1 = manufacture_solution(N, "1", "0", uph, uth, False)
test_problem(R, 1, ph_1, th_1, T, "u="+uph)
uph = "x"
uth = "1"
ph_1 = manufacture_solution(N, "uth*2", "uph", uph, uth, True)
th_1 = manufacture_solution(N, "1", "uth", uph, uth, False)
test_problem(R, 1/2, ph_1, th_1, T, "u="+uph)
uph = "x**2"
uth = "1-x"
ph_1 = manufacture_solution(N, "1", "uph", uph, uth, True)
th_1 = manufacture_solution(N, "1", "0", uph, uth, False)
test_problem(R, 1/3, ph_1, th_1, T, "u="+uph)
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

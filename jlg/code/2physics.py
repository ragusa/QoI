#!/usr/bin/env python3
"""
2physics.py

Solves a set of coupled diffusion equations of the form
du(x,t)/dt - dk(x,t,ui)du(x,t)/dx^2 + S(x,t,ui)*u(x,t) = Q(x,t)
"""
import numpy as np
from numpy import vectorize
from numpy.linalg import norm
import scipy as sp
import scipy.sparse as spar
from scipy.sparse.linalg import spsolve
import scipy.optimize as opt
from scipy.integrate import dblquad
import sympy
import matplotlib.pyplot as plt
from copy import deepcopy

np.set_printoptions(linewidth=200)

class boundary_condition:
		def __init__(self, BC):
			self.ltype = BC[0]
			self.rtype = BC[2]
			self.left  = BC[1] if callable(BC[1]) else lambda t: BC[1]
			self.right = BC[3] if callable(BC[3]) else lambda t: BC[3]

class diffusion_system:
	def __init__(self, N, T, k, Dk, Q, BC, offset, sig=0, order=2):
		self.offset = offset
		self.N   = N
		self.dx  = 1/N
		self.T   = T
		self.dt  = 1/T
		self.x   = np.linspace(0,1,N+1) # X values of nodes (N+1)
		self.k   = vectorize(k  ) if callable(k)   else vectorize(lambda x,t,uph,uth: k  )
		self.Dk  = vectorize(Dk ) if callable(Dk)  else vectorize(lambda x,t,uph,uth: Dk )
		self.Q   = vectorize(Q  ) if callable(Q)   else vectorize(lambda x,t,uph,uth: Q  )
		self.sig = vectorize(sig) if callable(sig) else vectorize(lambda x,t,uph,uth: sig)
		
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
		self.M = self.generate_mass_matrix(N)
	def generate_mass_matrix(self, N):
		"""
		Mass matrix with a flux = u and at time = t
		"""
		h = 1/(N)
		# Create the mass matrix
		d_main = 4/3 * np.ones(N+1) # Diagonal Entries
		d_off = 1/3 * np.ones(N+1) # Off diagonals
		M = spar.spdiags([d_off, d_main, d_off], [-1,0,1], N+1, N+1, format='csr')
		M /= 2
		M[0,1] *= 2
		M[-1,-2] *= 2
		return M*h/2
	def stiffness(self, t, uph, uth):
		"""
		Generate the stiffness matrix -∇k(t,uph,uth)∇ + Σ(t,uph,uth)
		"""
		A = spar.csr_matrix((self.N+1,self.N+1))
		for i in range(self.qorder):
			f[i] = sum(self.wq*self.b[:,i])
		for k in range(self.N):
			x0 = self.x[k]
			x1 = self.x[k+1]
			jacob = (x1-x0)/2
			# Quadrature points mapped to the element
			x = (x1+x0)/2 + self.xq*jacob
			# u at quadrature points
			uphq = np.dot(self.b, uph[k:k+self.qorder])
			uthq = np.dot(self.b, uth[k:k+self.qorder])

			kq =   self.k(x,t,uphq,uthq)/jacob # k at quadrature points
			Sq = self.sig(x,t,uphq,uthq)*jacob # Σ at quadrature points
			Alet = np.empty((self.qorder,self.qorder))
			for i in range(self.qorder):
				for j in range(self.qorder):
					Alet[i,j] = sum(self.wq*self.dbdx[:,i]*kq*self.dbdx[:,j])\
					          + sum(self.wq*Sq*self.b[:,i]*self.b[:,j])
			A[k:k+self.qorder, k:k+self.qorder] += Alet
		return A
	def residual(self, u_cur, u_prev, t_cur, dt, verbose=False):
		dx = 1/self.N
		#M(new-old)/dt - (Anew + Aold)/2 - (Qnew + Qold)/2
		t0 = t_cur - dt
		t1 = t_cur
		u0ph = u_prev[:self.N+1]
		u1ph =  u_cur[:self.N+1]
		u0th = u_prev[self.N+1:]
		u1th =  u_cur[self.N+1:]
		if self.offset==0:
			this0 = u0ph
			this1 = u1ph
			othr0 = u0th
			othr1 = u1th
		else:
			othr0 = u0ph
			othr1 = u1ph
			this0 = u0th
			this1 = u1th
		#M(Δu)/dt
		Mudt = self.M.dot(this1-this0)/dt
		F = np.zeros(self.N+1)

		k0s =   self.k(self.x,t0,u0ph,u0th)
		k1s =   self.k(self.x,t1,u1ph,u1th)
		S0s = self.sig(self.x,t0,u0ph,u0th)
		S1s = self.sig(self.x,t1,u1ph,u1th)
		Q0s =   self.Q(self.x,t0,u0ph,u0th)
		Q1s =   self.Q(self.x,t1,u1ph,u1th)
		for i in range(self.N):
			x0 = self.x[i]
			x1 = self.x[i+1]
			x01 = self.x[i:i+self.qorder]
			jacob = (x1-x0)/2
			# Quadrature points mapped to the element
			x = (x1+x0)/2 + self.xq*jacob
			# u values at quadrature points
			u0q = np.dot(self.b, this0[i:i+self.qorder]) # This Physics
			u1q = np.dot(self.b, this1[i:i+self.qorder]) # This Physics
			uq  = (u0q + u1q)/2
			ph0q = np.dot(self.b, u0ph[i:i+self.qorder])
			ph1q = np.dot(self.b, u1ph[i:i+self.qorder])
			th0q = np.dot(self.b, u0th[i:i+self.qorder])
			th1q = np.dot(self.b, u1th[i:i+self.qorder])
			# du/dx values at quadrature points
			du0dxq = np.dot(self.dbdx, this0[i:i+self.qorder])
			du1dxq = np.dot(self.dbdx, this1[i:i+self.qorder])
			dudxq = (du0dxq + du1dxq)/2
			# k(x,u) values at quadrature points
			k0q  =  np.dot(self.b, k0s[i:i+self.qorder])#self.k(x, t0, ph0q, th0q)
			k1q  =  np.dot(self.b, k1s[i:i+self.qorder])#self.k(x, t1, ph1q, th1q)
			kq   = (k0q + k1q)/2
			# Sigma(x,u) values at quadrature points
			Su0q  = np.dot(self.b, S0s[i:i+self.qorder]*this0[i:i+self.qorder])#self.sig(x, t0, ph0q, th0q)*jacob
			Su1q  = np.dot(self.b, S1s[i:i+self.qorder]*this1[i:i+self.qorder])#self.sig(x, t1, ph1q, th1q)*jacob
			Suq  = (Su0q+Su1q)/2
			# q(x) at the quadrature points
			Q0q = np.dot(self.b, Q0s[i:i+self.qorder])#self.Q(x, t0, ph0q, th0q)*jacob
			Q1q = np.dot(self.b, Q1s[i:i+self.qorder])#self.Q(x, t1, ph1q, th1q)*jacob
			Qq  = (Q0q + Q1q)/2
			# M(Δu)/dt at the quadrature points
			Mudtq = np.dot(self.b, Mudt[i:i+self.qorder])
			local = np.zeros(self.qorder)
			if self.offset==1 and i==1 and verbose:
				print(Suq - Qq)
			for j in range(self.qorder):
				local[j]	= np.dot(Mudtq*self.wq, self.b[:,j])\
							+ np.dot(  Suq*self.wq, self.b[:,j])*jacob\
							- np.dot(   Qq*self.wq, self.b[:,j])*jacob\
							+ np.dot(   kq*self.wq*self.dbdx[:,j], dudxq)/jacob
				if self.offset==1 and j!=i and verbose:
					print(t_cur, i, '--------------------')
					print("t", np.dot(Mudtq*self.wq, self.b[:,j])         )
					print("S", np.dot(  Suq*self.wq, self.b[:,j])*jacob   )
					print("Q",-np.dot(   Qq*self.wq, self.b[:,j])*jacob   )
					print("k", np.dot(   kq*self.wq*self.dbdx[:,j], dudxq)/jacob)
					print(local[j])
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
	def solve_adjoint(self, uph, uth, response):
		"""
		Takes the forward solutions uph, uth and the response function R(x,t)
		"""
		R = vectorize(response) # For easy mass solving.
		us = np.empty((self.T+1,(self.N+1)*2)) # Adjoint solution [t,x]

		# Solve tcur = T-1 .. 0
		for tcur in range(T,-1,-1):
			'''
			Right hand side special cases
			'''
			if tcur == T:
				RHS = R(self.x, tcur/self.T)
			elif tcur == T-1:
				RHS = R(self.x, tcur/self.T) - 1/2 * (self.M - self.stiffness(tcur/self.T, uph[tcur,:], uth[tcur,:]).dot(u[tcur+1]))
			elif tcur == 0:
				RHS = R(self.x, tcur/self.T) -   2 * (self.M - self.stiffness(tcur/self.T, uph[tcur,:], uth[tcur,:]).dot(u[tcur+1]))
			else:
				RHS = R(self.x, tcur/self.T) -       (self.M - self.stiffness(tcur/self.T, uph[tcur,:], uth[tcur,:]).dot(u[tcur+1]))
			'''
			Left hand side special cases
			'''
			if tcur == T:
				LHS = self.M + self.stiffness(1, uph[-1,:], uth[-1,:])
			elif tcur == 0:
				LHS = spar.eye((self.N+1))
			else:
				LHS = self.M + self.stiffness(1, uph[-1,:], uth[-1,:])
			'''
			Apply boundary conditions
			'''
			
			'''
			Solve
			'''
			us[tcur,:] = spsolve(LHS,RHS)
		return us

class FEM_SYS:
	def __init__(self, N, T, photon_system, thermal_system, order=2):
		# Discretation
		self.N = N # Number of spatial steps
		self.T = T # Number of temporal steps
		self.dt = 1/T # Time step width
		self.x = np.linspace(0,1,N+1) # X values of nodes (N+1)
		
		# Diffusion systems, creates a constant zero system if it is None
		self.ph =  photon_system if  photon_system is not None \
			else diffusion_system(N, T, 1, 0, 0, (2,lambda t: 1,2,lambda t: 1), 0)
		self.th = thermal_system if thermal_system is not None \
			else diffusion_system(N, T, 1, 0, 0, (2,lambda t: 1,2,lambda t: 1), 0)
		
		# Quadrature points used for integration
		# Weights of quadrature
		self.qorder = order
		[self.xq, self.wq] = np.polynomial.legendre.leggauss(self.qorder)
		# Basis functions at quadrature points
		self.b = np.array([(1-self.xq)/2, (1+self.xq)/2])
		# Derivative of the basis functions at quadrature points
		self.dbdx = np.array([[-.5, .5],[-.5, .5]])
	def forward_solve(self, u0, lastT=False):
		if lastT:
			u = np.copy(u0)
			for t in range(1, self.T+1):
				u = self.solve(u,u,t/self.T)
			return u
		u_for = np.empty((self.T+1,(self.N+1)*2))
		u_for[0,:] = u0
		for t in range(1, self.T+1):
			u_for[t,:] = self.solve(u_for[t-1,:],u_for[t-1,:],t/self.T)
		return u_for
	def newton_residual(self, u_cur, u_old, t, verbose=False):
		"""
		Calculate the residual
		Takes:
		u_cur, current guess
		u_old, solution of previous time step
		t, current time
		Returns
		F the residual vector
		"""
		if u_old is None:
			u_old = np.copy(u_cur)
		F = np.zeros((self.N+1)*2)
		F[0:self.N+1] = self.ph.residual(u_cur, u_old, t, self.dt, verbose)
		F[self.N+1: ] = self.th.residual(u_cur, u_old, t, self.dt, verbose)
		return F
	def solve(self, u_guess, u_old, t):
		"""
		Use newton's method to solve using the initial guess u0
		"""
		if u_guess is None:
			u_guess = np.zeros((self.N+1)*2)
		ans = opt.leastsq(
				func = self.newton_residual,
				x0 = u_guess,
				args = (u_old,t),
				ftol=1E-15
				)
		return ans[0]

def manufacture_solution(N, T, uph_str, uth_str, k_str, s_str, is_ph):
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
	u = uph if is_ph else uth
	uv = up if is_ph else ut
	# Create Sympy expressions
	# Create  k(x,t,uph,uth) function
	k = parse_expr(k_str)
	# Create Dk(x,t,uph,uth) function
	Dk = diff(k, uv)
	# Create  Σ(x,t,uph,uth) function
	Sig = parse_expr(s_str)
	# Create  Q(x,t,uph,uth) function
	ksub =   k.subs(up,uph).subs(ut,uth)
	Ssub = Sig.subs(up,uph).subs(ut,uth)
	Q = (diff(u,t) - diff(ksub*diff(u,x),x) + Ssub*u)
	#print("Q",Q)
	#print("k",k,  ' : ', ksub)
	#print("S",Sig,' : ', Ssub)
	#print("u",u)
	# Create BC by solving u for x=0,1
	ufun = sympy.lambdify((x,t),u,"numpy")
	bc_left  = lambda t: ufun(0,t)
	bc_right = lambda t: ufun(1,t)
	
	# Lambdify for construction of diffusion_system
	kfun     = sympy.lambdify((x,t,up,ut),k,"numpy")
	Dkfun    = sympy.lambdify((x,t,up,ut),Dk,"numpy")
	Qfun     = sympy.lambdify((x,t,up,ut),Q,"numpy")
	sigfun   = sympy.lambdify((x,t,up,ut),Sig,"numpy")
	BC       = (2,bc_left,2,bc_right)
	
	# Create diffusion system
	offset = 0 if is_ph else 1
	d_sys = diffusion_system(N, T, kfun, Dkfun, Qfun, BC, offset, sigfun)

	u0 = np.array([ufun(x,0) for x in np.linspace(0,1,N+1)])
	return d_sys, u0, ufun
def CN_integrate(u,response,N=None,T=None,order=2):
	'''
	Integrates u[t,x] against the response R(x,t), assuming crank-nicolson in time
	'''
	if T == None:
		T = len(u[:,0])
	if N == None:
		N = len(u[0,:]) # Guesses spatial size using the first timestep
	R = vectorize(response)
	dt = 1/T
	dx = 1/N
	[xq, wq] = np.polynomial.legendre.leggauss(order)
	b = np.array([(1-xq)/2, (1+xq)/2])
	QoI = 0
	for t in range(T+1):
		# Temporal Integral
		QoI_part = 0
		for i in range(N):
			# Spatial Integral
			x0 = i/N
			x1 = x0 + dx
			jacob = dx/2
			xs = (x0+x1)/2 + xq*jacob
			QoI_part += jacob * dt * sum(R(xs,t/T)*wq*np.dot(b,u[t,i:i+order]))
		if t == 0 or t==T:
			QoI_part *= 0.5 # Adjustment for Crank Nicolson
		QoI += QoI_part
	return QoI
def L2Err(uh, ue_fun, response, N=None,T=None,order=2):
	if T == None:
		T = len(uh)-1
	if N == None:
		N = len(uh[0,:])-1 # Guesses spatial size using the first timestep
	ue = np.array([[ue_fun(x/N,t/T) for x in range(N+1)] for t in range(T+1)])
	R = vectorize(response)
	dt = 1/T
	dx = 1/N
	[xq, wq] = np.polynomial.legendre.leggauss(order)
	b = np.array([(1-xq)/2, (1+xq)/2])
	L2 = 0
	for t in range(T+1):
		# Temporal Integral
		# Adjustment for Crank Nicolson
		if t == 0 or t == T:
			mod = 0.5
		else:
			mod = 1
		for i in range(N):
			# Spatial Integral
			x0 = i/N
			x1 = x0 + dx
			jacob = dx/2
			xs = (x0+x1)/2 + xq*jacob
			I_h = jacob * dt * dx * sum(wq*R(xs,t/T)*np.dot(b,uh[t,i:i+order]))
			I_e = jacob * dt * dx * sum(wq*R(xs,t/T)*np.dot(b,ue[t,i:i+order]))
			L2 += (mod*I_h - mod*I_e)**2
	return L2

def make_sys(N, T, PH, TH):
	ph, uph0, uphfun = manufacture_solution(N, T, PH[0], TH[0], PH[1], PH[2], True)
	th, uth0, uthfun = manufacture_solution(N, T, PH[0], TH[0], TH[1], TH[2], False) 
	u = FEM_SYS(N, T, ph, th)
	u0 = np.concatenate((uph0,uth0))
	u1 = np.concatenate(([uphfun(x,1) for x in u.ph.x], [uthfun(x,1) for x in u.th.x]))
	return u,u0,u1, uphfun, uthfun
def test_problem(ID, N, T, R_ph, R_th, PH, TH):
	u,u0,u1,uphfun,uthfun = make_sys(N, T, PH, TH)
	solp = dblquad(uphfun,0,1,lambda x: 0,lambda x: 1)[0]
	solt = dblquad(uthfun,0,1,lambda x: 0,lambda x: 1)[0]
	res = u.forward_solve(u0)
	qoifph = CN_integrate(res[:,:N+1], R_ph, N, T)
	qoifth = CN_integrate(res[:,N+1:], R_th, N, T)
	errph = L2Err(res[:,:N+1], uphfun, R_ph)
	errth = L2Err(res[:,N+1:], uthfun, R_th)
	absErrph = abs(solp - qoifph)
	absErrth = abs(solt - qoifth)
	print("%16s | %8E %8E %8E %8E | %8E %8E %8E %8E"%(ID, solp, qoifph, errph, absErrph, solt, qoifth, errth, absErrth))
	#print(res)
	#print(u.newton_residual(u1,u0,1,verbose=True))
	return u, res, qoifph, errph, qoifth, errth
def Richardson(Series, Steps, p, known, disp=True):
	L = len(Series)
	R = np.zeros((L,L))
	R[:,0] = Series
	for r in range(1,L):
		for i in range(r,L):
			t = Steps[i-1]/Steps[i]
			k = p+r-1
			R[i,r] = (t**k*R[i,r-1] - R[i-1,r-1])/(t**k-1)
	if disp:
		print(R)
		print(abs(known-R))
		print(known)
	return(R[-1,-1])
print("%16s | %12s %12s %12s %12s | %12s %12s %12s %12s"%("Problem","QOI(PH)","","L2(PH)","AbsErr(PH)","QOI(TH)","","L2(TH)","AbsErr(TH)"))
print("-"*(124))
R_ph = lambda x,t: 1
R_th = lambda x,t: 1
kph = "uth"
sph = "uth"
kth = "uph"
sth = "uth"

# Constant
uph = "1"
uth = "1"
T = 2
N = 2
test_problem("Constant", N, T, R_ph, R_th, (uph, kph, sph), (uth, kth, sth))
# Linear in Space
uph = "x"
uth = "x"
T = 1
N = 2
test_problem("Linear x", N, T, R_ph, R_th, (uph, kph, sph), (uth, kth, sth))
# Linear in Time
uph = "t"
uth = "t"
T = 2
N = 2
test_problem("Linear t", N, T, R_ph, R_th, (uph, kph, sph), (uth, kth, sth))
# Linear in both
uph = "t"
uth = "x"
T = 1
N = 2
test_problem("Linear t&x", N, T, R_ph, R_th, (uph, kph, sph), (uth, kth, sth))
# Quadratic in both
uph = "t**2"
uth = "x**2"
T = 2
N = 2
test_problem("Quadradic t&x", N, T, R_ph, R_th, (uph, kph, sph), (uth, kth, sth))
# Cubic in both
uph = "t**3"
uth = "x**3"
T = 20
N = 20
test_problem("Cubic t&x", N, T, R_ph, R_th, (uph, kph, sph), (uth, kth, sth))
# Sines Cosine
uph = "x*sin(t)"
uth = "x*cos(t)"
T = 10
N = 3
test_problem("x*trig(t)", N, T, R_ph, R_th, (uph, kph, sph), (uth, kth, sth))


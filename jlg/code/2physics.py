#!/usr/bin/env python3
"""
2physics.py

Solves a set of coupled diffusion equations of the form
du(x,t)/dt - dk(x,t,ui)du(x,t)/dx^2 + S(x,t,ui)*u(x,t) = Q(x,t)
"""
import numpy as np
import scipy as sp
import scipy.sparse as spar
import scipy.optimize as opt
from numpy.linalg import norm
import sympy
import matplotlib.pyplot as plt

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
		M = spar.spdiags([d_off, d_main, d_off], [-1,0,1], N+1, N+1, format='csr')
		M /= 2
		M[0,1] *= 2
		M[-1,-2] *= 2
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
			ph_uq = np.dot(self.b, u0[k:k+self.qorder])
			th_uq = np.dot(self.b, u0[k+self.N+1:k+self.N+1+self.qorder])

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
		dx = 1/self.N
		#M(new-old)/dt - (Anew + Aold)/2 - (Qnew + Qold)/2
		t0 = t_cur - dt
		t1 = t_cur
		if self.offset==0:
			this1 =  u_cur[:self.N+1] # This Physics
			this0 = u_prev[:self.N+1] # This Physics
			othr1 =  u_cur[self.N+1:]   # Other Physics
			othr0 = u_prev[self.N+1:]  # Other Physics
			u0ph = this0
			u1ph = this1
			u0th = othr0
			u1th = othr1
		else:
			this1 =  u_cur[self.N+1:]  # This Physics
			this0 = u_prev[self.N+1:] # This Physics
			othr1 =  u_cur[:self.N+1] # Other Physics 
			othr0 = u_prev[:self.N+1] # Other Physics
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
			uq  = (u0q + u1q)/2
			ph0q = np.dot(self.b, u0ph[i:i+self.qorder])
			ph1q = np.dot(self.b, u1ph[i:i+self.qorder])
			th0q = np.dot(self.b, u0th[i:i+self.qorder])
			th1q = np.dot(self.b, u1th[i:i+self.qorder])
			# du/dx values at quadrature points
			du0dxq = np.dot(self.dbdx, this0[i:i+self.qorder])/jacob
			du1dxq = np.dot(self.dbdx, this1[i:i+self.qorder])/jacob
			dudxq = (du0dxq + du1dxq)/2
			# k(x,u) values at quadrature points
			k0q  =  self.k(x, t0, ph0q, th0q)
			k1q  =  self.k(x, t1, ph1q, th1q)
			kq  = (k0q + k1q)/2
			# Sigma(x,u) values at quadrature points
			S0q = self.sig(x, t0, ph0q, th0q)*jacob
			S1q = self.sig(x, t1, ph1q, th1q)*jacob
			Suq = (S0q*u0q + S1q*u1q)/2
			# q(x) at the quadrature points
			Q0q =   self.Q(x, t0, ph0q, th0q)*jacob
			Q1q =   self.Q(x, t1, ph1q, th1q)*jacob
			Qq  = (Q0q + Q1q)/2
			# M(Δu)/dt at the quadrature points
			Mudtq = np.dot(self.b, Mudt[i:i+self.qorder])
			local = np.zeros(self.qorder)
			for j in range(self.qorder):
				local[j]	= np.dot(Mudtq*self.wq, self.b[:,j])\
							+ np.dot(  Suq*self.wq, self.b[:,j])\
							- np.dot(   Qq*self.wq, self.b[:,j])\
							+ np.dot(   kq*self.wq*self.dbdx[:,j], dudxq)
				#print(dudxq)
				#if self.offset == 0:
					##print(j, 'T Term', np.dot(Mudtq*self.wq, self.b[:,j]))
					##print(j, 'S Term', np.dot(  Suq*self.wq, self.b[:,j]))
					#print(j, 'k Term',  np.dot(   kq*self.wq*self.dbdx[:,j], dudxq))
					#print(j, 'Q Term', -np.dot(   Qq*self.wq, self.b[:,j]))
					#print(j, local[j])
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
		self.ph =  photon_system if  photon_system is not None \
			else diffusion_system(N, 1, 0, 0, (2,lambda t: 1,2,lambda t: 1), 0)
		self.th = thermal_system if thermal_system is not None \
			else diffusion_system(N, 1, 0, 0, (2,lambda t: 1,2,lambda t: 1), 0)
		
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
	def newton_residual(self, u_cur, u_old, t):
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
		F[0:self.N+1] = self.ph.residual(u_cur, u_old, t, self.dt)
		F[self.N+1: ] = self.th.residual(u_cur, u_old, t, self.dt)
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

def manufacture_solution(N, uph_str, uth_str, k_str, s_str, is_ph):
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
	Dk = diff(k, up) if is_ph else diff(k, ut)
	# Create  Σ(x,t,uph,uth) function
	Sig = parse_expr(s_str)
	# Create  Q(x,t,uph,uth) function
	Q = diff(u,t) - k*diff(diff(u,x),x) - diff(k,uv)*diff(u,x)**2 + Sig*u
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
	d_sys = diffusion_system(N, kfun, Dkfun, Qfun, BC, offset, sigfun)

	u0 = np.array([ufun(x,0) for x in np.linspace(0,1,N+1)])
	return d_sys, u0, ufun
def innerproduct(u, res, N, T, lastT=False):
	tau = 1/T
	[xq, wq] = np.polynomial.legendre.leggauss(2)
	b = np.array([(1-xq)/2, (1+xq)/2])
	qoif_t = np.zeros(T+1)
	for t in range(T+1):
		if lastT:
			ut = u
		else:
			ut = u[t,:]
		if not lastT or t==T:
			for i in range(N):
				x0 = i/N
				x1 = (i+1)/N
				jacob = (x1-x0)/2
				x = (x1+x0)/2 + xq*jacob
				qoif_t[t] += jacob * sum(wq*np.dot(b, ut[i:i+2])*res(x,t/T))
	return sum(qoif_t[1::] + qoif_t[0:-1])
def L2Err(N, tcur, ue, u, R):
	err = 0
	dx = 1/N
	xs = np.linspace(0,1,N+1) # X values of nodes (N+1)
	# Quadrature points used for integration
	# Weights of quadrature
	qorder = 2
	[xq, wq] = np.polynomial.legendre.leggauss(qorder)
	# Basis functions at quadrature points
	b = np.array([(1-xq)/2, (1+xq)/2])
	udiff = ue - u
	for i in range(N+0):
		x0 = xs[i]
		x1 = xs[i+1]
		dx = x1 - x0
		jacob = (x1-x0)/2
		x = (x1+x0)/2 + xq*jacob
		err += dx/2 * np.dot(wq, R*(x,tcur)*np.dot(b, udiff[i:i+2]))**2
	return np.sqrt(err)
def test_problem(ID, N, T, R_ph, R_th, PH, solph, TH, solth):
	ph, uph0, uphfun = manufacture_solution(N, PH[0], TH[0], PH[1], PH[2], True)
	th, uth0, uthfun = manufacture_solution(N, PH[0], TH[0], TH[1], TH[2], False) 
	u = FEM_SYS(N, T, ph, th)
	u0 = np.concatenate((uph0,uth0))
	u1 = np.concatenate(([uphfun(x,1) for x in u.ph.x], [uthfun(x,1) for x in u.th.x]))
	res = u.forward_solve(u0,lastT=True)
	qoifph = innerproduct(res, R_ph, N, T, lastT=True)
	qoifth = innerproduct(res, R_th, N, T, lastT=True)
	errph = L2Err(N, 1, [uphfun(x,1) for x in u.ph.x], res[0:N+1], R_ph)
	errth = L2Err(N, 1, [uthfun(x,1) for x in u.th.x], res[N+1:], R_th)
	print("%12s | %8E %8E %8E %8E"%(ID, qoifph, errph, qoifth, errth))
	#print(u.newton_residual(u1,u0,1))
	return u, res, qoifph, errph, qoifth, errth


print("%12s | %12s %12s %12s %12s"%("Problem","QOI(PH)","Err(PH)","QOI(TH)","Err(TH)"))
print("-"*(81))
R_ph = lambda x,t: 0
R_th = lambda x,t: 1 if t==1 else 0
uph = "t"
kph = "1"
sph = "uth"
uth = "t"
kth = "uth"
sth = "0"
T = 1
N = 2
test_problem("Linear", N, T, R_ph, R_th, (uph, kph, sph), 0, (uph, kph, sph), 1)

print("-"*(81))
Errs = []
for T in [1,10,100,1000]:
	N = 2
	uph = "t**3"
	uth = "t**3"
	u, u1, qoifph, errph, qoifth, errth = test_problem(str(T)+" T Cubic", N, T, R_ph, R_th, (uph, kph, sph), 0, (uph, kph, sph), 1)
	Errs.append((T,errth))
conv = np.log10(Errs[-2][1]/Errs[-1][1])/np.log10(Errs[-1][0]/Errs[-2][0])
print("Rate of Convergence ", conv)

print("-"*(81))
Errs = []
for N in [1,10,100,1000]:
	T = 1
	uph = "x**3"
	uth = "x**3"
	u, u1, qoifph, errph, qoifth, errth = test_problem(str(N)+" S Cubic", N, T, R_ph, R_th, (uph, kph, sph), 0, (uph, kph, sph), 1/4)
	Errs.append((N,errth))
conv = np.log10(Errs[-2][1]/Errs[-1][1])/np.log10(Errs[-1][0]/Errs[-2][0])
print("Rate of Convergence ", conv)
print("-"*(81))

print("-"*(81))
R_ph = lambda x,t: 0
R_th = lambda x,t: 1 if t==1 else 0
kph = "1"
sph = "uth"
kth = "uth"
sth = "0"
T = 1
N = 2
u, unpert, qoifph, errph, qoifth, errth = test_problem("Linear", N, T, R_ph, R_th, (sph, kph, sph), 0, (sph, kph, sph), 1)

pert = 0.00001
pert_kph = u
pert_kph.ph.k = lambda x,t,ph,th: u.ph.k(x,t,ph,th) + pert
u_kph = pert_kph.forward_solve(u0,lastT=True)
sens_kph = innerproduct(u_kph, R_th, N, T, lastT=True)/pert

pert_sph = u
pert_sph.ph.k = lambda x,t,ph,th: u.ph.sig(x,t,ph,th) + pert
u_sph = pert_sph.forward_solve(u0,lastT=True)
sens_sph = innerproduct(u_sph, R_th, N, T, lastT=True)/pert

pert_kth = u
pert_kth.th.k = lambda x,t,ph,th: u.th.k(x,t,ph,th) + pert
u_kth = pert_kth.forward_solve(u0,lastT=True)
sens_kth = innerproduct(u_kth, R_th, N, T, lastT=True)/pert

pert_sth = u
pert_sth.th.k = lambda x,t,ph,th: u.th.sig(x,t,ph,th) + pert
u_sth = pert_sth.forward_solve(u0,lastT=True)
sens_sth = innerproduct(u_tph, R_th, N, T, lastT=True)/pert

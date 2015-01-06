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
from scipy.sparse.linalg import inv as spinv
import scipy.optimize as opt
from scipy.integrate import dblquad
import sympy
import matplotlib.pyplot as plt
from copy import deepcopy

np.set_printoptions(linewidth=200)

class boundary_condition:
		def __init__(self, BC):
			self.ltype = BC[0]
			self.left  = BC[1] if callable(BC[1]) else lambda t: BC[1]
			self.rtype = BC[2]
			self.right = BC[3] if callable(BC[3]) else lambda t: BC[3]
def quadmaker(N, p):
	'''
	Return quadrature points, weights, basis functions, and basis derivatives
	for order p
	'''
	[xq, wq] = np.polynomial.legendre.leggauss(p)
	xd = np.linspace(-1,1,p)
	b    = np.zeros((p,p))
	dbdx = np.zeros((p,p))
	gn = np.zeros((N+1,p),dtype='int')
	gn[0,:] = [i for i in range(p)]
	for i in range(1,N+1):
		gn[i,:] = gn[i-1]+1
	for i in range(p):
		NUM = 1
		DEN = 1
		for j in range(p):
			if j != i:
				NUM *= xq - xd[j]
				DEN *= xd[i] - xd[j]
		b[:,i] = NUM/DEN
	for i in range(p):
		SUM = 0
		DEN = 1
		for j in range(p):
			if j != i:
				NUM = 1
				for k in range(p):
					if k!=i and k!=j:
						NUM *= xq - xd[k]
				SUM += NUM
				DEN *= xd[i]-xd[j]
		dbdx[:,i] = SUM/DEN
	#print(xq)
	#print(wq)
	#print(b)
	#print(dbdx)
	#print(gn)
	return xq,wq,b,dbdx,gn
class diffusion_system:
	def __init__(self, N, T, k, Dk, Q, BC, offset, sig=0, order=2):
		self.offset = offset
		self.N   = N
		self.dx  = 1/N
		self.T   = T
		self.dt  = 1/T
		self.x   = np.linspace(0,1,N+1) # X values of nodes (N+1)
		self.k   = vectorize(k  ) if callable(k)   else vectorize(lambda x,t,uph,uth,duphdx: k  )
		self.Dk  = vectorize(Dk ) if callable(Dk)  else vectorize(lambda x,t,uph,uth,duphdx: Dk )
		self.Q   = vectorize(Q  ) if callable(Q)   else vectorize(lambda x,t,uph,uth,duphdx: Q  )
		self.sig = vectorize(sig) if callable(sig) else vectorize(lambda x,t,uph,uth,duphdx: sig)
		self.cond = 0 #Thermal conductivity
		# Quadrature points used for integration
		# Quadrature
		self.qorder = order
		[self.xq, self.wq, self.b, self.dbdx, self.gn] = quadmaker(self.N, self.qorder)
		# Boundary Conditions
		self.BC = boundary_condition(BC)
		# Generate and save the mass matrix
		self.M = self.generate_mass_matrix(N)
		self.iM = spinv(self.M)
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
		return M/2
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
		dup0dx = 0
		dup1dx = 0
		#k0s =   self.k(self.x,t0,u0ph,u0th,dup0dx)
		#k1s =   self.k(self.x,t1,u1ph,u1th,dup1dx)
		S0s = self.sig(self.x,t0,u0ph,u0th,dup0dx)
		S1s = self.sig(self.x,t1,u1ph,u1th,dup1dx)
		Q0s =   self.Q(self.x,t0,u0ph,u0th,dup0dx)
		Q1s =   self.Q(self.x,t1,u1ph,u1th,dup1dx)
		for i in range(self.N):
			x0 = self.x[i]
			x1 = self.x[i+1]
			x01 = self.x[i:i+self.qorder]
			jacob = (x1-x0)/2
			# Quadrature points mapped to the element
			xq = (x1+x0)/2 + self.xq*jacob
			# u values at quadrature points
			u0q = np.dot(self.b, this0[self.gn[i]])
			u1q = np.dot(self.b, this1[i:i+self.qorder])
			uq  = (u0q + u1q)/2
			ph0q = np.dot(self.b, u0ph[i:i+self.qorder])
			ph1q = np.dot(self.b, u1ph[i:i+self.qorder])
			th0q = np.dot(self.b, u0th[i:i+self.qorder])
			th1q = np.dot(self.b, u1th[i:i+self.qorder])
			# du/dx values at quadrature points
			du0dxq = np.dot(self.dbdx, this0[i:i+self.qorder])
			du1dxq = np.dot(self.dbdx, this1[i:i+self.qorder])
			dudxq = (du0dxq + du1dxq)/2
			dph0dxq = np.dot(self.dbdx, u0ph[i:i+self.qorder])
			dph1dxq = np.dot(self.dbdx, u1ph[i:i+self.qorder])
			# k(x,u) values at quadrature points
			#k0q  =  np.dot(self.b, k0s[i:i+self.qorder]*du0dxq)
			#k1q  =  np.dot(self.b, k1s[i:i+self.qorder]*du1dxq)
			#kq   = (k0q**-1 + k1q**-1)**-1/2
			# Sigma(x,u) values at quadrature points
			Su0q = np.dot(self.b, self.sig(xq,t0,ph0q,th0q,du0dxq)*(th0q**4 - ph0q))
			Su1q = np.dot(self.b, self.sig(xq,t1,ph1q,th1q,du1dxq)*(th1q**4 - ph1q))
			Suq  = (Su0q+Su1q)/2
			# q(x) at the quadrature points
			Q0q = np.dot(self.b, self.Q(xq,t0,ph0q,th0q,du0dxq))
			Q1q = np.dot(self.b, self.Q(xq,t1,ph1q,th1q,du1dxq))
			Qq  = (Q0q + Q1q)/2
			# M(Δu)/dt at the quadrature points
			Mudtq = np.dot(self.b, Mudt[i:i+self.qorder])
			local = np.zeros(self.qorder)
			# The Paper's
			if self.offset == 0:
				kq = (3*((th0q+th1q)/2)**-3 + (dph1dxq+dph0dxq)/(ph1q+ph0q))**-1#
			if self.offset == 1:
				kq = ((th0q+th1q)/2)**2.5
			for j in range(self.qorder):
				local[j]	= sum(self.wq * Mudtq * self.b[:,j])               \
							- sum(self.wq *   Suq * self.b[:,j])         *jacob\
							- sum(self.wq *    Qq * self.b[:,j])         *jacob\
							+ sum(self.wq *    kq * self.dbdx[:,j])/jacob
				if self.offset==0 and verbose:
					print(t_cur, i, '--------------------')
					print(S0s)
					print(S1s)
					print(Suq)
					print("t", sum(self.wq * Mudtq * self.b[:,j])               )
					print("S",-sum(self.wq *   Suq * self.b[:,j])         *jacob)
					print("Q",-sum(self.wq *    Qq * self.b[:,j])         *jacob)
					print("k", sum(self.wq *    kq * self.dbdx[:,j]*dudxq)/jacob)
					print(local[j])
			F[i:i+self.qorder] += local
		# Apply Boundary Conditions
		if self.BC.ltype == 0: # Neumann
			F[0] -= self.BC.left(t1)
		elif self.BC.ltype  == 1: # Robin
			F[0] += 1/4 * this1[0] - (1/6/S1s[0])*self.BC.left(t1)
		elif self.BC.ltype  == 2: # Dirichlet
			F[0] = this1[0] - self.BC.left(t1)
		if self.BC.rtype == 0: # Neumann
			F[-1] -= self.BC.right(t1)
		elif self.BC.rtype == 1: # Robin
			F[-1] += 1/4 * this1[-1] + (1/6/S1s[-1])*self.BC.right(t1)
		elif self.BC.rtype == 2: # Dirichlet
			F[-1] = this1[-1] - self.BC.right(t1)
		return F
	def res2(self, u, t, verbose=False):
		dx = 1/self.N
		limiter = 2
		m = 1
		z = 1
		F = np.zeros(self.N+1)
		for i in range(self.N):
			gnE = np.array([i,i+1])
			gnT = gnE + self.N+1
			x0 = self.x[i]
			x1 = self.x[i+1]
			jacob = (x1-x0)/2
			
			local_E    = np.dot(self.b, u[gnE])
			local_T    = np.dot(self.b, u[gnT])
			local_dEdx = np.dot(self.dbdx, u[gnE])
			local_dTdx = np.dot(self.dbdx, u[gnT])
			k_T = self.cond * local_T**2.5
			SigA = z / local_T**3
			E = sum(self.wq * local_E)/2
			gradE = np.absolute(local_dEdx)/jacob
			k_E = 1/((3*SigA)**m + (limiter*gradE/E)**m)**(1/m)
			local_F = np.zeros(self.qorder)
			aux1 = SigA * self.wq*(local_T**4-local_E)
			aux2 =  k_T * self.wq*(local_dTdx)
			aux3 =  k_E * self.wq*(local_dEdx)
			for j in range(self.qorder):
				if self.offset == 0:
					local_F[j]  = -np.dot(aux1, self.b[:,j]   )*jacob
					local_F[j] +=  np.dot(aux3, self.dbdx[:,j])/jacob
					if i==0 and j==1:
						print(i,SigA[j], local_E[j], aux1[j])
					if i==1 and j==0:
						print(i,SigA[j], local_E[j], aux1[j])
				elif self.offset == 1:
					local_F[j]  =  np.dot(aux1, self.b[:,j]   )*jacob
					local_F[j] +=  np.dot(aux2, self.dbdx[:,j])/jacob
			F[i:i+self.qorder] += local_F
		# Apply BCs
		if   self.BC.ltype == 0:
			F[0]  -= self.BC.left(t)
		elif self.BC.ltype == 1:
			F[0]  += 1/2*u[self.offset*(self.N+1)] - 2*self.BC.left(t)
		elif self.BC.ltype == 2:
			F[0]   = u[self.offset*(self.N+1)] - self.BC.left(t)
		if   self.BC.rtype == 0:
			F[-1] -= self.BC.right(t)
		elif self.BC.rtype == 1:
			F[-1] += 1/2*u[self.N+self.offset*(self.N+1)] - 2*self.BC.right(t)
		elif self.BC.rtype == 2:
			F[-1]  = u[self.N+self.offset*(self.N+1)] - self.BC.right(t)
		#print(F)
		#x = input("pause")
		return -F
	def solve_adjoint(self, uph, uth, response):
		"""
		Takes the forward solutions uph, uth and the response function R(x,t)
		"""
		raise NotImplementedError("The Adjoint is not implemented")
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
	def __init__(self, N, T, photon_system, thermal_system):
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
		self.M = spar.csr_matrix((2*(self.N+1),2*(self.N+1)))
		self.M[:N+1,:N+1] = self.ph.M
		self.M[N+1:,N+1:] = self.th.M
		#self.iM = spinv(self.M)
		self.iM = spar.csr_matrix((2*(self.N+1),2*(self.N+1)))
		self.iM[:N+1,:N+1] = self.ph.iM
		self.iM[N+1:,N+1:] = self.th.iM
	def forward_solve(self, u0, lastT=False):
		if lastT:
			u = np.copy(u0)
			for t in range(1, self.T+1):
				u = self.solve(u,u,t/self.T)
			return u
		u_for = np.empty((self.T+1,(self.N+1)*2))
		u_for[0,:] = u0
		for t in range(1, 2):#self.T+1):
			u_for[t,:] = self.solve(u_for[t-1,:],u_for[t-1,:],t/self.T)
		return u_for
	def res2(self, u,t):
		N = self.N
		F = np.zeros(2*(N+1))
		F[:N+1] = self.ph.iM.dot(self.ph.res2(u,t))
		F[N+1:] = self.th.iM.dot(self.th.res2(u,t))
		return F
	def CrankNicholson(self,u0,dt):
		u_for = np.empty((self.T+1,(self.N+1)*2))
		u_for[0,:] = u0
		print('Beginning CN')
		for tstep in range(1, self.T+1):
			t = tstep*dt
			u_for[tstep,:] = self.CrankNicholsonStep(u_for[tstep-1,:], t-dt, dt)
		return u_for
	def CrankNicholsonStep(self, u0, t0, dt):
		g = 0.5
		t1 = t0+dt

		Res0 = self.res2(u0,t0)
		G = lambda utry: (utry-u0) - 1/2 * dt * (Res0 + self.res2(utry,t1))
		
		R0 = norm(G(u0))
		tol = max(1E-15, 1E-10*R0)
		u_ans = opt.leastsq(func = G, x0 = u0, ftol=tol)
		
		R1 = norm(G(u_ans[0]))
		print('t = %.2E completed %.2E %.2E'%(t1,R0,R1))
		return u_ans[0]
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
			u_guess = np.copy(u_old)
		ans = opt.leastsq(
				func = self.newton_residual,
				x0 = u_guess,
				args = (u_old,t),
				ftol=1E-15
				)
		print(self.newton_residual(ans[0],u_old,t))
		return ans[0]

def manufacture_solution(N, T, uph_str, uth_str, k_str, s_str, is_ph):
	"""
	Constructs a diffusion_system using sympy to create the source term
	Takes k, Σ, and the solution as strings
	"""
	from sympy.parsing.sympy_parser import parse_expr
	from sympy import diff
	# Declare Sympy symbols
	t, x, up, ut, dpdx = sympy.symbols('t x uph uth dpdx')
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
	Q = diff(u,t) - diff(ksub*diff(u,x),x) - Ssub*(uth**4 - uph)
	print("Q",Q)
	print("k",diff(ksub*diff(u,x),x))
	print("S",Ssub*(uth**4 - uph))
	print("u",diff(u,t),u)
	# Create BC by solving u for x=0,1
	ufun = sympy.lambdify((x,t),u,"numpy")
	bc_left  = lambda t: ufun(0,t)
	bc_right = lambda t: ufun(1,t)
	
	# Lambdify for construction of diffusion_system
	kfun     = sympy.lambdify((x,t,up,ut,dpdx),k,"numpy")
	Dkfun    = sympy.lambdify((x,t,up,ut,dpdx),Dk,"numpy")
	Qfun     = sympy.lambdify((x,t,up,ut,dpdx),Q,"numpy")
	sigfun   = sympy.lambdify((x,t,up,ut,dpdx),Sig,"numpy")
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
	print(res)
	print(u.newton_residual(u1,u0,1,verbose=True))
	return u, res, qoifph, errph, qoifth, errth
def Tests():
	print("%16s | %12s %12s %12s %12s | %12s %12s %12s %12s"%("Problem","QOI(PH)","","L2(PH)","AbsErr(PH)","QOI(TH)","","L2(TH)","AbsErr(TH)"))
	print("-"*(124))
	R_ph = lambda x,t: 1
	R_th = lambda x,t: 1
	# Constants used: z=12
	# k = 0
	kph = "1/3/(12/uth)**3"
	sph = "(12/uth)**3"
	kth = "0"
	sth = "-(12/uth)**3"
	# Constant
	uph = "1"
	uth = "1"
	T = 2
	N = 2
	test_problem("Constant", N, T, R_ph, R_th, (uph, kph, sph), (uth, kth, sth))
	# Linear in Time
	uph = "t+1"
	uth = "t+1"
	T = 2
	N = 2
	test_problem("Linear t+1", N, T, R_ph, R_th, (uph, kph, sph), (uth, kth, sth))
	# Linear in Space
	uph = "x+1"
	uth = "x+1"
	T = 2
	N = 2
	test_problem("Linear x+1", N, T, R_ph, R_th, (uph, kph, sph), (uth, kth, sth))
	# Linear in both
	uph = "x+t+1"
	uth = "x+t+1"
	T = 2
	N = 2
	test_problem("Linear t+x+1", N, T, R_ph, R_th, (uph, kph, sph), (uth, kth, sth))
	# Quadratic in both
	uph = "t**2+1"
	uth = "x**2+1"
	T = 2
	N = 2
	test_problem("Quadradic t&x+1", N, T, R_ph, R_th, (uph, kph, sph), (uth, kth, sth))
	# Cubic in both
	uph = "t**3+1"
	uth = "x**3+1"
	T = 10
	N = 10
	test_problem("Cubic t&x+1", N, T, R_ph, R_th, (uph, kph, sph), (uth, kth, sth))
	# Sines Cosine
	uph = "x*sin(t)+1"
	uth = "x*cos(t)+1"
	T = 10
	N = 3
	test_problem("x*trig(t)+1", N, T, R_ph, R_th, (uph, kph, sph), (uth, kth, sth))
def WaveProblem():
	T = 10
	N = 10
	k = 0  # Constant the defines the material, conductivity
	z = 1  # Constant the defines the material, atomic mass number
	SIGPH = lambda x,t,uph,uth,duphdx: (1/uth)**3
	SIGTH = lambda x,t,uph,uth,duphdx: -SIGPH(x,t,uph,uth,duphdx)
	K_TH  = lambda x,t,uph,uth,duphdx: k * uth**(5/2)
	K_PH  = lambda x,t,uph,uth,duphdx: 1/(3*SIGPH(x,t,uph,uth,duphdx))# + 1/uph*abs(duphdx))
	BC_PH = (1,1,1,0)
	BC_TH = (0,0,0,0)
	
	UP = diffusion_system(N, T, K_PH, 0, 0, BC_PH, 0, SIGPH)
	UT = diffusion_system(N, T, K_TH, 0, 0, BC_TH, 1, SIGTH)
	U  = FEM_SYS(N, T, UP, UT)
	dt = 1E-8
	U0 = 1E-8 * np.ones(2*(N+1))
	U0[N+1:] = U0[N+1:]**0.25
	U_FOR = U.CrankNicholson(U0,dt)
	fig = plt.figure()
	ax = fig.add_subplot(121)
	for t in range(0,T+1,T//10):
		ax.plot(U_FOR[t,N+1:], 'o-', label='Time %.2E'%(t*dt))
	hand, labels = ax.get_legend_handles_labels()
	ax.legend(hand, labels)
	ax = fig.add_subplot(122)
	for t in range(0,T+1,T//10):
		ax.plot(U_FOR[t,:N+1], 'o-', label='Time %.2E'%(t*dt))
	hand, labels = ax.get_legend_handles_labels()
	ax.legend(hand, labels)
	plt.show()

#Tests()

## Linear in Space
#R_ph = lambda x,t: 1
#R_th = lambda x,t: 1
## Constants used: z=12
## k = 0
#kph = "1/3/(12/uth)**3"
#sph = "(12/uth)**3"
#kth = "0"
#sth = "-(12/uth)**3"
#uph = "x+1"
#uth = "x+1"
#T = 1000
#N = 2
#test_problem("Linear x+1", N, T, R_ph, R_th, (uph, kph, sph), (uth, kth, sth))
WaveProblem()

#!/usr/bin/env python3

import numpy as npy
from numpy.linalg import solve
import scipy as spy
from scipy import sparse as sp
from FEM2 import femsys

def make_mass_matrix(N):
	h = 1/(N)
	# Create the mass matrix
	d_main = 4/3 * npy.ones(N+1) # Diagonal Entries
	d_off = 1/3 * npy.ones(N+1) # Off diagonals
	M = sp.spdiags([d_off, d_main, d_off], [-1,0,1], N+1, N+1, format='csc')
	M[0,0] /= 2
	M[-1,-1] /= 2
	return M*h/2

def getAb(Dfun, dDfun, SIGfun, Qfun, bc_fun, t, u):
	#Functions solved at time t
	D = lambda uvect, x: Dfun(uvect, x, t)
	dD = lambda uvect, x: dDfun(uvect, x, t)
	sigma = lambda uvect, x: SIGfun(uvect, x, t)
	Q = lambda x: Qfun(x, t)
	bc = bc_fun(t)
	fem = femsys(N, 1, D, dD, sigma, Q, bc)
	[A_bc, b_bc, A, b] = fem.assemble(npy.array(u))
	return (A,b,bc)
def make_operators(N, T, A, b, bc):
	h = 1/(N)
	tau = 1/(T)
	M = make_mass_matrix(N)
	B1 = -M - A * tau/2
	B2 = M - A * tau/2 # Diagonal Block
	Q = npy.matrix(b).transpose()
	return (B1, B2, Q)
def apply_bc(B1, B2, Q, bc_old, bc):
	if bc[0] == 2:
		# Dirichlet
		# Off Diagonal
		Q -= bc_old[1] * B1[:,0]
		B1[0,:] = 0
		B1[:,0] = 0
		
		# Diagonal
		Q -= bc[1] * (B2[:,0])
		B2[0,:] = 0
		B2[:,0] = 0
		B2[0,0] = 1
		
		Q[0] = bc[1]
	if bc[2] == 2:
		# Dirichlet
		# Off Diagonal
		Q -= bc_old[3] * (B1[:,-1])
		B1[-1,:] = 0
		B1[:,-1] = 0
		
		# Diagonal
		Q -= bc[3] * (B2[:,-1])
		B2[-1,:] = 0
		B2[:,-1] = 0
		B2[-1,-1] = 1
		
		Q[-1] = bc[3]
	return(B1,B2,Q)
def CN(N, T, u0, Dfun, dDfun, SIGfun, Qfun, res, bc_fun, direction=1):
	h = 1/(N)
	tau = 1/(T)
	u = npy.matrix(npy.zeros((N+1,T+1)))
	QQ = npy.matrix(npy.zeros((N+1,T+1)))
	# Time end points
	time0 = 1   if direction==1 else T-1
	timef = T+1 if direction==1 else -1
	
	A,b,bc = getAb(Dfun, dDfun, SIGfun, Qfun, bc_fun, 0, u)
	B1,B2,Q = make_operators(N, T, A, b, bc)
	
	# Set initial
	if direction == 1:
		u[:,0] = npy.matrix(u0).T
		QQ[:,0] = Q
	elif direction == -1:
		(Junk, A_L, Junk) = apply_bc(B1, B2, Q, bc_fun(T-1), bc_fun(T))
		u[:,-1] = solve(A_L,Q)
		QQ[:,-1] = Q
	# Time march
	for t in range(time0, timef, direction):
		# Operators
		if direction == 1:
			(A_R, A_L, Q_cur) = apply_bc(B1, B2, Q, bc_fun(t-direction), bc_fun(t))
			
			RHS = Q_cur - npy.dot(A_R, u[:,t-direction])
			LHS = A_L
		elif direction == -1:
			(A_R, Junk, Junk) = apply_bc(B1, B2, Q, bc_fun(t), bc_fun(t+1))
			(Junk, A_L, Junk) = apply_bc(B1, B2, Q, bc_fun(t-1), bc_fun(t))
			Q_cur = Q
			# Weird adjoint modifications
			if t==T-1:
				A_R*=0.5
			if t==0:
				A_R*=2
				A_L = npy.identity(N+1)
			RHS = tau**2*(Q_cur - npy.dot(A_R, u[:,t-direction]))
			LHS = tau**2*(A_L)
		u[:,t] = solve(LHS,RHS)
		QQ[:,t] = Q_cur
	return (u,QQ)
def CNadjoint(N, T, u0, Dfun, dDfun, SIGfun, Qfun, res, bc_fun):
	
	return CN(N, T, u0, Dfun, dDfun, SIGfun, res, Qfun, bc_fun, -1)

def innerproduct(u, res, T, N):
	tau = 1/T
	[xq, wq] = npy.polynomial.legendre.leggauss(2)
	b = npy.array([(1-xq)/2, (1+xq)/2])
	qoif_t = npy.zeros(T+1)
	for t in range(T+1):
		ut = u[:,t]
		for i in range(N):
			x0 = i/N
			x1 = (i+1)/N
			jacob = (x1-x0)/2
			x = (x1+x0)/2 + xq*jacob
			qoif_t[t] += jacob * sum(wq*npy.dot(b, ut[i:i+2])*res(x,t))
	return tau/2*sum(qoif_t[1::] + qoif_t[0:-1])
def dotproduct(x, y, T, N):
	h = 1/N
	tau = 1/T
	r = 0
	for n in range(N+1):
		for t in range(T+1):
			ra = tau*h*x[n,t]*y[n,t]
			if t == 0:
				ra *= .5
			elif t == T:
				ra *= .5
			r += ra
	return r
print('u\tQoif\tError\t\tQoia\tError')


N = 15
T = 1

# u = 1
Dfun = lambda u, x, t: 1
dDfun = lambda u, x, t: 0
SIGfun = lambda u, x, t: 0
bc_fun = lambda t: (2, 1, 2, 1)
Qfun = lambda x, t: 0
res = lambda x, t: 1

u0 = npy.ones(N+1)
(u, b) = CN(N, T, u0, Dfun, dDfun, SIGfun, Qfun, res, bc_fun)
qoif = innerproduct(u, res, T, N)
us0 = npy.ones(N+1)
(us, r) = CNadjoint(N, T, us0, Dfun, dDfun, SIGfun, Qfun, res, bc_fun)
qoia = dotproduct(us, b, T, N)
print('1\t%.4f\t%.4e\t%.4f\t%.4e'%(qoif, abs(qoif-1), qoia, abs(qoif-qoia)))

# u = t
Dfun = lambda u, x, t: 1
dDfun = lambda u, x, t: 0
SIGfun = lambda u, x, t: 0
bc_fun = lambda t: (2, 0, 2, t/T)
Qfun = lambda x, t: 1
res = lambda x, t: 1

u0 = npy.zeros(N+1)
(u, b) = CN(N, T, u0, Dfun, dDfun, SIGfun, Qfun, res, bc_fun)
qoif = innerproduct(u, res, T, N)
us0 = npy.ones(N+1)
(us, r) = CNadjoint(N, T, us0, Dfun, dDfun, SIGfun, Qfun, res, bc_fun)
qoia = dotproduct(us, b, T, N)
print('t\t%.4f\t%.4e\t%.4f\t%.4e'%(qoif, abs(qoif-1/2), qoia, abs(qoif-qoia)))


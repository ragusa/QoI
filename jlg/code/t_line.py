#!/usr/bin/env python3

import numpy as np
import numpy.matlib as matrix
from numpy.linalg import inv
from numpy.linalg import solve
from FEM2 import femsys
def apply_boundary(LHS,RHS,bc):
	# Dirichlet only
	if bc[0] == 2:
		RHS -= bc[1]*LHS[:,0]
		LHS[0,:] = 0
		LHS[:,0] = 0
		LHS[0,0] = 1
		RHS[0] = bc[1]
	elif bc[0] == 0:
		RHS += bc[1]
	if bc[2] == 2:
		RHS -= bc[3]*LHS[:,-1]
		LHS[-1,:] = 0
		LHS[:,-1] = 0
		LHS[-1,-1] = 1
		RHS[-1] = bc[3]
	elif bc[2] == 0:
		RHS += bc[3]
	return (LHS,RHS)
def getAb(Dfun, dDfun, SIGfun, qf, bc_fun, t, u):
	#Functions solved at time t
	D = lambda uvect, x: Dfun(uvect, x, t)
	dD = lambda uvect, x: dDfun(uvect, x, t)
	sigma = lambda uvect, x: SIGfun(uvect, x, t)
	Q = lambda x: qf(x, t)
	bc = bc_fun(t)
	fem = femsys(N, 1, D, dD, sigma, Q, bc)
	[Anew_bc, bnew_bc, Anew, bnew] = fem.assemble(np.transpose(u))
	return (Anew,bnew,bc)
def CN(N, T, u0, Dfun, dDfun, SIGfun, Qfun, bc_fun):
	'''
	Integrates du/dt - ∇·D(u,x)∇u(x) + Σ(u,x)u(x) = Q(x)
	'''
	h = 1/N
	tau = 1/T
	u = np.zeros((N+1,T+1))
	b = np.zeros((N+1,T+1))
	u[:,0] = u0
	I = np.identity(N+1)
	M = np.identity(N+1)
	for i in range(N+1):
		M[i,i] = 2
		if i > 0:
			M[i,i-1] = 1
		if i < 0:
			M[i,i+1] = 1
		M[i,:] /= sum(M[i,:])
	M = np.identity(N+1)
	(Anew,bnew,bc) = getAb(Dfun, dDfun, SIGfun, Qfun, bc_fun, 0, u[:,0])
	for t in range(1,T+1):
		Aold = np.copy(Anew)
		bold = np.copy(bnew)
		(Anew,bnew,bc) = getAb(Dfun, dDfun, SIGfun, Qfun, bc_fun, t, u[:,t])
		QQ = tau/2 * (bold + bnew)
		RHS = (QQ + np.dot(h*M + Aold*tau/2, u[:,t-1]))
		LHS = h*M - Anew*tau/2
		(LHS, RHS) = apply_boundary(LHS,RHS,bc)
		u[:,t] = solve(LHS, RHS)
		b[:,t] = QQ
	return (u, b)

def CNadjoint(N, T, u0, Dfun, dDfun, SIGfun, Qfun, Rfun, bc_fun):
	h = 1/N
	tau = 1/T
	u = np.zeros((N+1,T+1))
	b = np.zeros((N+1,T+1))
	u[:,T] = u0
	I = np.identity(N+1)
	M = np.identity(N+1)
	for i in range(N+1):
		M[i,i] = 2/3
		if i > 0:
			M[i,i-1] = 1/3
		if i < 0:
			M[i,i+1] = 1/3
	(A,b_res,bc) = getAb(Dfun, dDfun, SIGfun, Rfun, bc_fun, T, u[:,T])
	(Anew,bnew,bc) = getAb(Dfun, dDfun, SIGfun, Qfun, bc_fun, T, u[:,T])
	Anew = Anew.transpose()
	K = tau * np.ones(T+1) # diagonal, don't bother with 2d
	K[0] *= 1/2
	K[-1] *= 1/2
	Kinv = tau * np.ones(T+1) # diagonal, don't bother with 2d
	Kinv[0] *= 2
	Kinv[-1] *= 2
	
	(LHS, RHS) = apply_boundary(A,b_res,bc)
	u[:,-1] = solve(LHS, RHS)
	for t in range(T-1,-1,-1):
		Aold = np.copy(Anew)
		bold = np.copy(bnew)
		(Anew,bnew,bc) = getAb(Dfun, dDfun, SIGfun, Qfun, bc_fun, t, u[:,t])
		Anew = Anew.transpose()
		QQ = b_res
		RHS = (QQ - np.dot(K[t+1]*Kinv[t]*(h*M + Anew*tau/2), u[:,t+1]))
		print(QQ)
		LHS = K[t]*Kinv[t]*(h*M - Anew*tau/2)
		(LHS, RHS) = apply_boundary(LHS,RHS,bc)
		u[:,t] = solve(LHS, RHS)
		b[:,t] = QQ
	return (u, b)

def innerproduct(u, res, T, N):
	tau = 1/T
	[xq, wq] = np.polynomial.legendre.leggauss(2)
	b = np.array([(1-xq)/2, (1+xq)/2])
	qoif_t = np.zeros(T+1)
	for t in range(T+1):
		ut = u[:,t]
		for i in range(N):
			x0 = i/N
			x1 = (i+1)/N
			jacob = (x1-x0)/2
			x = (x1+x0)/2 + xq*jacob
			qoif_t[t] += jacob * sum(wq*np.dot(b, ut[i:i+2])*res(x,t))
	return tau/2*sum(qoif_t[1::] + qoif_t[0:-1])
def dotproduct(u, r, T, N):
	h = 1/N
	tau = 1/T
	qoi = 0
	for t in range(T+1):
		for n in  range(N+1):
			qoi += u[n,t] * b[n,t] * h
	return qoi
# du/dt - ∇·D(u,x)∇u(x) + Σ(u,x)u(x) = Q(x)
N = 10 # Spatial Cells
T = 1 # Number of time cells

tau = 1/T
h = 1/N

print('u\tQoif\t\tError\t\tQoia')
# u = x t
Dfun = lambda u, x, t: 1
dDfun = lambda u, x, t: 0
SIGfun = lambda u, x, t: 0
bc_fun = lambda t: (2, 0, 2, t/T)
Qfun = lambda x, t: x
res = lambda x,t: 1
u0 = np.zeros(N+1)
(u, b) = CN(N, T, u0, Dfun, dDfun, SIGfun, Qfun, bc_fun)
qoif = innerproduct(u, res, T, N)
us0 = np.ones(N+1)
(us, r) = CNadjoint(N, T, us0, Dfun, dDfun, SIGfun, Qfun, res, bc_fun)
#qoia = -innerproduct(us, Qfun, T, N) + h*(sum(u[:,-1]*us[:,-1])-sum(u[:,0]*us[:,0]))
qoia = dotproduct(us, b, T, N)
print('x t\t%.4f\t%.4e\t%.4f\t%.4e'%(qoif, abs(qoif-1/4), qoia, abs(qoif-qoia)))

# u = x t^2
Dfun = lambda u, x, t: 1
dDfun = lambda u, x, t: 0
SIGfun = lambda u, x, t: 0
bc_fun = lambda t: (2, 0, 2, (t/T)**2)
abc_fun = lambda t: (2, 1, 2, 1)
Qfun = lambda x, t: x*2*t/T
res = lambda x,t: 1
u0 = np.zeros(N+1)
(u, b) = CN(N, T, u0, Dfun, dDfun, SIGfun, Qfun, bc_fun)
qoif = innerproduct(u,res, T, N)
us0 = np.ones(N+1)
(us, r) = CNadjoint(N, T, us0, Dfun, dDfun, SIGfun, Qfun, res, bc_fun)
#qoia = -innerproduct(us, Qfun, T, N) + h*(sum(u[:,-1]*us[:,-1])-sum(u[:,0]*us[:,0]))
qoia = dotproduct(us, b, T, N)
print('x t^2\t%.2f\t%.2e\t%.2f\t%.2e'%(qoif, abs(qoif-1/6), qoia, abs(qoif-qoia)))

# u = x t^3
Dfun = lambda u, x, t: 1
dDfun = lambda u, x, t: 0
SIGfun = lambda u, x, t: 0
bc_fun = lambda t: (2, 0, 2, (t/T)**3)
abc_fun = lambda t: (2, 1, 2, 1)
Qfun = lambda x, t: x*3*(t/T)**2
res = lambda x,t: 1
u0 = np.zeros(N+1)
(u, b) = CN(N, T, u0, Dfun, dDfun, SIGfun, Qfun, bc_fun)
qoif = innerproduct(u,res, T, N)
us0 = np.ones(N+1)
(us, r) = CNadjoint(N, T, us0, Dfun, dDfun, SIGfun, Qfun, res, bc_fun)
#qoia = -innerproduct(us, Qfun, T, N) + h*(sum(u[:,-1]*us[:,-1])-sum(u[:,0]*us[:,0]))
qoia = dotproduct(us, b, T, N)
print('x t\^3t%.2f\t%.2e\t%.2f\t%.2e'%(qoif, abs(qoif-1/8), qoia, abs(qoif-qoia)))
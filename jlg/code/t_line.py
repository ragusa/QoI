#!/usr/bin/env python3

import numpy as np
import numpy.matlib as matrix
from numpy.linalg import inv
from numpy.linalg import solve
from FEM2 import femsys

def CN(N, T, u0, Dfun, dDfun, SIGfun, Qfun, bc_fun):
	'''
	Integrates du/dt - ∇·D(u,x)∇u(x) + Σ(u,x)u(x) = Q(x)
	'''
	h = 1/N
	tau = 1/T

	u = np.zeros((N+1,T+1))
	b = np.zeros((N+1,T+1))
	u[:,0] = u0

	D0 = lambda uvect, x: Dfun(uvect, x, 0)
	dD0 = lambda uvect, x: dDfun(uvect, x, 0)
	sigma0 = lambda uvect, x: SIGfun(uvect, x, 0)
	Q0 = lambda x: Qfun(x, 0)
	bc0 = bc_fun(0)
	fem = femsys(N, h, D0, dD0, sigma0, Q0, bc0)
	
	[Anew, bnew, Anew_nbc, bnew_nbc] = fem.assemble(np.transpose(u[:,0]))
	
	for t in range(1,T+1):
		#Functions solved at time t
		D = lambda uvect, x: Dfun(uvect, x, t)
		dD = lambda uvect, x: dDfun(uvect, x, t)
		sigma = lambda uvect, x: SIGfun(uvect, x, t)
		Q = lambda x: Qfun(x, t)
		bc = bc_fun(t)
		#Assemble A, b
		M = np.identity(N+1)
		M *= tau
		M[0][0] *= 0.5
		M[-1][-1] *= 0.5
		Aold = np.copy(Anew)
		bold = np.copy(bnew)
		fem = femsys(N, h, D, dD, sigma, Q, bc)
		[Anew, bnew, Anew_nbc, bnew_nbc] = fem.assemble(np.transpose(u[:,t-1]))
		old = -M + Aold*tau/2
		new = M + Anew*tau/2
		QQ = tau/2 * (bold + bnew)
		u[:,t] = solve(new, (QQ - np.dot(old, u[:,t-1])))
		b[:,t] = QQ
	return (u, b)
def CNadjoint(N, T, u0, Dfun, dDfun, SIGfun, Qfun, bc_fun):
	'''
	Integrates du/dt - ∇·D(u,x)∇u(x) + Σ(u,x)u(x) = Q(x)
	'''
	h = 1/N
	tau = 1/T

	u = np.zeros((N+1,T+1))
	u[:,-1] = u0

	D0 = lambda uvect, x: Dfun(uvect, x, 0)
	dD0 = lambda uvect, x: dDfun(uvect, x, 0)
	sigma0 = lambda uvect, x: SIGfun(uvect, x, 0)
	Q0 = lambda x: Qfun(x, 0)
	bc0 = bc_fun(0)
	fem = femsys(N, h, D0, dD0, sigma0, Q0, bc0)
	
	[Anew_bc, bnew_bc, Anew, bnew] = fem.assemble(np.transpose(u[:,-1]))
	
	M = np.identity(N+1)
	M *= tau
	M[0][0] *= 0.5
	M[-1][-1] *= 0.5
	# Apply BC (Dirichlet)
	bnew -= bc_fun(0)*A[:,0]
	Anew[0,:] = 0
	Anew[:,0] = 0
	Anew[0,0] = 1
	bnew[0] = bc_fun(0)
	bnew -= bc_fun(1)*A[:,-1]
	Anew[-1, :] = 0
	Anew[ :,-1] = 0
	Anew[-1,-1] = 1
	bnew[-1] = bc_fun(1)
	for t in range(T-1,-1,-1):
		#Functions solved at time t
		D = lambda uvect, x: Dfun(uvect, x, t)
		dD = lambda uvect, x: dDfun(uvect, x, t)
		sigma = lambda uvect, x: SIGfun(uvect, x, t)
		Q = lambda x: Qfun(x, t)
		bc = bc_fun(t)
		#Assemble A, b
		bold = np.copy(bnew)
		fem = femsys(N, h, D, dD, sigma, Q, bc)
		[Anew, bnew, Anew_nbc, bnew_nbc] = fem.assemble(np.transpose(u[:,t+1]))
		old = -M + Anew*tau/2
		new = M + Anew*tau/2
		QQ = tau*h*(bnew+bold)
		u[:,t] = solve(new, (QQ - np.dot(old, u[:,t+1])))
	return u

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

# du/dt - ∇·D(u,x)∇u(x) + Σ(u,x)u(x) = Q(x)
N = 10 # Spatial Cells
T = 10 # Number of time cells

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
qoif = innerproduct(u,res, T, N)
us0 =np.zeros(N+1)
us = CNadjoint(N, T, us0, Dfun, dDfun, SIGfun, res, bc_fun)
qoia = sum(sum(us*b))*h/2
print('x t', qoif, abs(qoif-1/4), qoia, abs(qoif-qoia), sep='\t')

# u = x t^2
Dfun = lambda u, x, t: 1
dDfun = lambda u, x, t: 0
SIGfun = lambda u, x, t: 0
bc_fun = lambda t: (2, 0, 2, (t/T)**2)
Qfun = lambda x, t: x*2*t/T
res = lambda x,t: 1
u0 = np.zeros(N+1)
(u, b) = CN(N, T, u0, Dfun, dDfun, SIGfun, Qfun, bc_fun)
qoif = innerproduct(u,res, T, N)
us0 = np.zeros(N+1)
us = CNadjoint(N, T, us0, Dfun, dDfun, SIGfun, res, bc_fun)
qoia = sum(sum(us*b))*h/2
print('x t^2', qoif, abs(qoif-1/6), qoia, abs(qoif-qoia), sep='\t')

# u = x t^3
Dfun = lambda u, x, t: 1
dDfun = lambda u, x, t: 0
SIGfun = lambda u, x, t: 0
bc_fun = lambda t: (2, 0, 2, (t/T)**3)
Qfun = lambda x, t: x*3*(t/T)**2
res = lambda x,t: 1
u0 = np.zeros(N+1)
(u, b) = CN(N, T, u0, Dfun, dDfun, SIGfun, Qfun, bc_fun)
qoif = innerproduct(u,res, T, N)
us0 = np.zeros(N+1)
us = CNadjoint(N, T, us0, Dfun, dDfun, SIGfun, res, bc_fun)
qoia = sum(sum(us*b))*h/2
print('x t^3', qoif, abs(qoif-1/8), qoia, abs(qoif-qoia), sep='\t')

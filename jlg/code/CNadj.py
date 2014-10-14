#!/usr/bin/env python3

import numpy as np
from numpy.linalg import solve
from scipy import sparse as sp
from FEM2 import femsys
import matplotlib.pyplot as plt

def mass_matrix(N,T):
	'''
	Creates a sparse diagonal matrix of the mass matrix.
	'''
	h = 1/(N)
	# Create the mass matrix
	d_main = 4/3 * np.ones(N+1) # Diagonal Entries
	d_off = 1/3 * np.ones(N+1) # Off diagonals
	M = sp.spdiags([d_off, d_main, d_off], [-1,0,1], N+1, N+1, format='csc')
	M[0,0] /= 2
	M[-1,-1] /= 2
	return M*h/2

def make_operators(D, dD, Sigma, BC, Q, R, N, T, N_WIDTH, T_WIDTH, u):
	'''
	(d/dt - ∇·D∇ + Σ) u = Q
	Operators are
	Diagonal = M - tau/2 * A
	Off Diagonal = M + tau/2 * A
	'''
	tau = T_WIDTH/T
	M = mass_matrix(N, T)
	
	# Reuse FEM code for spatial portion.
	FEM = femsys(N, N_WIDTH, D, dD, Sigma, Q, BC)
	ignored_A_bc, ignored_b_bc, A, rhs = FEM.assemble(u)
	
	# Create the operators
	on_diagonal  =   M - tau/2*A
	off_diagonal = - M - tau/2*A
	rhs *= tau
	
	return (np.array(on_diagonal), np.array(off_diagonal), rhs)

def apply_bc(on_diagonal, off_diagonal, rhs, bc, t_cur, t_old):
	bc_cur = bc(t_cur)
	bc_prev = bc(t_old)
	if bc_cur [0] == 2:
		# Left Dirichlet
		rhs -= bc_prev[1] * off_diagonal[:,0]
		off_diagonal[0,:] = 0
		off_diagonal[:,0] = 0
		
		rhs -= bc_cur[1] * on_diagonal[:,0]
		on_diagonal[0,:] = 0
		on_diagonal[:,0] = 0
		
		on_diagonal[0,0] = 1
		rhs[0] = bc_cur [1]
	if bc_cur [2] == 2:
		# Right Dirichlet
		rhs -= bc_prev[3] * off_diagonal[:,-1]
		off_diagonal[-1,:] = 0
		off_diagonal[:,-1] = 0
		
		rhs -= bc_cur[3] * on_diagonal[:,-1]
		on_diagonal[-1,:] = 0
		on_diagonal[:,-1] = 0
		
		on_diagonal[-1,-1] = 1
		rhs[-1] = bc_cur [3]
	return on_diagonal, off_diagonal, rhs

def forward(D_t, dD_t, Sigma_t, BC, Q_t, R_t, N, T, N_WIDTH, T_WIDTH, u0):
	u = np.empty((T+1, N+1))
	src = np.empty((T+1, N+1))
	u[0,:] = u0
	src[0,:] = u0
	for t in range(1, T+1):
		# Create Functions for the problem, at current time
		D = lambda u, x: D_t(u, x, t)
		dD = lambda u, x: dD_t(u, x, t)
		Sigma = lambda u, x: Sigma_t(u, x, t)
		Q = lambda x: Q_t(x, t)
		R = lambda x: R_t(x, t)
		
		# Create Operators
		on_diagonal, off_diagonal, rhs = make_operators(D, dD, Sigma, BC(t), Q, R, N, T, N_WIDTH, T_WIDTH, u[t-1,:])
		on, off, rhs = apply_bc(on_diagonal, off_diagonal, rhs, BC, t, t-1)
		# Solve
		u[t,:] = solve(on, rhs - np.dot(off, u[t-1,:]))
		src[t,:] = rhs
	return u, src

def adjoint(D_t, dD_t, Sigma_t, BC, Q_t, R_t, N, T, N_WIDTH, T_WIDTH, u0):
	u = np.empty((T+1, N+1))
	A_adj = np.zeros(((T+1)*(N+1),(T+1)*(N+1)))
	for t in range(T, -1,-1):
		# Create Functions for the problem, at current time
		D = lambda u, x: D_t(u, x, t)
		dD = lambda u, x: dD_t(u, x, t)
		Sigma = lambda u, x: Sigma_t(u, x, t)
		Q = lambda x: Q_t(x, t)
		R = lambda x: R_t(x, t)
		
		# Create Operators
		on_diagonal, off_diagonal, rhs = make_operators(D, dD, Sigma, BC(t), R, Q, N, T, N_WIDTH, T_WIDTH, u[t-1,:])
		on, off, rhs = apply_bc(on_diagonal, off_diagonal, rhs, BC, t, t)
		rhs = np.ones(N+1)
		
		if t == T:
			u[t,:] = solve(on, rhs)
		elif t == T-1:
			off *= 1/2
			u[t,:] = solve(on, rhs - np.dot(off, u[t+1,:]))
		elif t == 0:
			on = np.identity(N+1)
			off *= 2
			u[t,:] = solve(on, rhs - np.dot(off, u[t+1,:]))
		else:
			u[t,:] = solve(on, rhs - np.dot(off, u[t+1,:]))
	return u
def make_K(N,T):
	K = .5*np.identity((N+1)*(T+1))
	K[0:N+1,0:N+1] *= 0.5
	K[(T)*(N+1):(N+1)*(T+1),(T)*(N+1):(N+1)*(T+1)] *= 0.5
	K *= (1/T)
	return K
def better_print(A):
	for r in A:
		for c in r:
			print('%3E '%c, end='')
		print()
	print()
def innerproduct(u, res, N, T):
	tau = 1/T
	[xq, wq] = np.polynomial.legendre.leggauss(2)
	b = np.array([(1-xq)/2, (1+xq)/2])
	qoif_t = np.zeros(T+1)
	for t in range(T+1):
		ut = u[t,:]
		for i in range(N):
			x0 = i/N
			x1 = (i+1)/N
			jacob = (x1-x0)/2
			x = (x1+x0)/2 + xq*jacob
			qoif_t[t] += jacob * sum(wq*np.dot(b, ut[i:i+2])*res(x,t))
	return tau/2*sum(qoif_t[1::] + qoif_t[0:-1])
def dotproduct(x, y, N, T):
	K = make_K(N,T)
	r_x = x.reshape((N+1)*(T+1))
	r_y = y.reshape((N+1)*(T+1))
	Ky = np.dot(K,r_y)
	xKy = r_x*Ky
	h = 1/(N+1)
	return 2*h*sum(xKy)
def plot(U, N, T):
	fig = plt.figure()
	X = np.linspace(0,1,N+1)
	ax = fig.add_subplot(111)
	ax.set_ylim(-1,2)
	for u in U:
		ax.plot(X, u)
	plt.show()

D = lambda u, x, t: 1
dD = lambda u, x, t: 0
Sigma = lambda u, x, t: 0
R = lambda x, t: 1
print("%20s\t%12s\t%12s\t%12s\t%12s"%('Function', 'QoI_F', 'ERR_F', 'QoI_A', 'ERR_FA'))
def test_problem(tag, u0, BC, Q, N, T, qoi_known):
	u_for, src = forward(D, dD, Sigma, BC, Q, R, N, T, 1, 1, u0)
	qoif = innerproduct(u_for, R, N, T)
	u_adj = adjoint(D, dD, Sigma, BC, Q, R, N, T, 1, 1, u0)
	qoia =  dotproduct(u_adj, src, N, T)
	print("%3i %3i %12s\t%5E\t%5E\t%5E\t%5E"%(N, T, tag, qoif, abs(qoif - qoi_known), qoia, abs(qoia - qoif)))
	return qoif, abs(qoif - qoi_known), qoia, abs(qoia - qoif),qoia/qoif, u_for, u_adj
fig = plt.figure()
ax = fig.add_subplot(111)
'''Nmin = 3
Nmax = 10
maxE = []
for T in range(1,20):
	E = []
	for N in range(Nmin, Nmax+1):
		# u = 1
		u0 = np.ones(N+1)
		BC = lambda t: (2, 1, 2, 1)
		Q = lambda x, t: 0
		qoif, errf, qoia, errfa, err, u, us = test_problem('u = 1',u0, BC, Q, N, T, 1)
		E.append(errfa)
	ax.plot(np.linspace(Nmin,Nmax,(Nmax-Nmin+1)), E, label=T)
	maxE.append((N,T,max(E)))
print(maxE)
legend = ax.legend()
plt.show()'''

N = 100
T = 1
# u = 1
u0 = np.ones(N+1)
BC = lambda t: (2, 1, 2, 1)
Q = lambda x, t: 0
qoif, errf, qoia, errfa, err, u, us = test_problem('u = 1',u0, BC, Q, N, T, 1)
# u = t
u0 = np.zeros(N+1)
BC = lambda t: (2, t/T, 2, t/T)
Q = lambda x, t: 1
qoif, errf, qoia, errfa, err, u, us = test_problem('u = t',u0, BC, Q, N, T, 1/2)

# u = x^2
u0 = np.array([(x/N)**2 for x in range(N+1)])
BC = lambda t: (2, 0, 2, 1)
Q = lambda x, t: 2
qoif, errf, qoia, errfa, err, u, us = test_problem('u = x²',u0, BC, Q, N, T, 1/3)
	
# u = tx
u0 = np.zeros(N+1)
BC = lambda t: (2, 0, 2, t/T)
Q = lambda x, t: x
qoif, errf, qoia, errfa, err, u, us = test_problem('u = t x',u0, BC, Q, N, T, 1/4)


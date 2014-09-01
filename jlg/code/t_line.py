#!/usr/bin/env python3

import numpy as np
import numpy.matlib as matrix
from numpy.linalg import inv
from numpy.linalg import solve
from FEM2 import femsys

def t_line(A, M, Q, N, steps, u0):
	tau = 1/steps
	u = matrix.zeros((N+1,steps+1))
	u[:,0] = u0
	for t in range(0, steps):
		MMold = -M(tau*t) + tau/2 * A(tau*t)
		MMnew = M(tau*(t+1)) + tau/2 * A(tau*(t+1))
		QQ = tau/2 * (Q(tau*t)+Q(tau*(t+1)))
		u[:,t+1] = solve(MMnew, QQ - (MMold)*u[:,t])
	return u

N = 5 # Spatial Cellts
h = 1 # Width
t_end = 1 # final time
t_steps = 100 # Number of time cells
tau = 1/t_steps


# du/dt - ∇·D(u,x)∇u(x) + Σ(u,x)u(x) = Q(x)
sigma = lambda u, x: 1
D = lambda u, x: 1
dD = lambda u, x: 0
Q = lambda x: 1

# Use FEM2 to constuct A with BC, assumes linear
fem = femsys(N, h, D, dD, sigma, Q, (0,0,0,0))
(A, b, A_nobc, b_nobc) = fem.assemble(np.zeros(N+1))
Afun = lambda t: matrix.mat(A)
bfun = lambda t: np.transpose(matrix.mat(b))
M = matrix.eye(N+1)
Mfun = lambda t: M*tau/2
u0 = np.transpose(matrix.mat(solve(A,b)))

u = t_line(Afun, Mfun, bfun, N, t_steps, matrix.zeros((N+1,1)))
print(u)

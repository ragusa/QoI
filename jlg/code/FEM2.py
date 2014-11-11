#!/usr/bin/env python3
import numpy as np
import scipy as sp
from numpy.linalg import solve
def norm(x):
	return sum([i**2 for i in x])**0.5
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
class femsys:
	'''
	FEM system to solve
	-∇·D(u,x)∇u(x) + Σ(u,x)u(x)= q(x)
	'''
	def __init__(self, n, h, D, dD, s, Q, bc, order=2):
		self.diff = D # D(u,x)
		self.dDdu = dD # dD(u,x)/du
		self.siga = s # Σ(u,x)
		self.nel = n # Number of elements (len(u) - 1)
		self.esrc = Q # q(x)
		self.x = np.linspace(0,h,n+1) # X values of nodes (nel+1)
		self.width = h # Width of problem
		# Quadrature points used for integration
		# Weights of quadrature
		self.qorder = order
		[self.xq, self.wq] = np.polynomial.legendre.leggauss(self.qorder)
		# Basis functions at quadrature points
		self.b = np.array([(1-self.xq)/2, (1+self.xq)/2])
		# Derivative of the basis functions at quadrature points
		self.dbdx = np.array([[-.5, .5],[-.5, .5]])
		# Boundary Conditions
		self.bc = boundary(bc[0],bc[1],bc[2],bc[3])
	def assemble(self, u0):
		'''
		Assemble the system Au=b, for a given value of u
		'''
		A_nobc = np.zeros((self.nel+1,self.nel+1))
		b_nobc = np.zeros((self.nel+1))
		f = np.zeros(self.qorder)
		for i in range(self.qorder):
			f[i] = sum(self.wq*self.b[:,i])
		for k in range(self.nel):
			x0 = self.x[k]
			x1 = self.x[k+1]
			jacob = (x1-x0)/2
			# Quadrature points mapped to the element
			x = (x1+x0)/2 + self.xq*jacob
			# u and du/dx values at quadrature points
			uq = np.dot(self.b, u0[k:k+self.qorder])
			dudxq = np.dot(self.dbdx, u0[k:k+self.qorder])
			# D(x,u) values at quadrature points
			Dq = self.diff(uq, x)/jacob
			dDduq = self.dDdu(uq, x)
			# Sigma values at quadrature points
			Sq = self.siga(uq, x)*jacob
			Alet = np.zeros((self.qorder,self.qorder))
			for i in range(self.qorder):
				for j in range(self.qorder):
					Alet[i,j] = sum(self.wq*self.dbdx[:,i]*Dq*self.dbdx[:,j]) + sum(self.wq*Sq*self.b[:,i]*self.b[:,j])
			A_nobc[k:k+self.qorder,k:k+self.qorder] += Alet
			b_nobc[k:k+self.qorder] += self.esrc(x)*f*jacob
		A = np.copy(A_nobc)
		b = np.copy(b_nobc)
		# Apply Boundary Conditions
		if self.bc.left == 0: # Neumann
			b[0] += self.bc.lval
		elif self.bc.left == 1: # Robin
			b[0] += 2*self.bc.lval
			A[0,0] += 1/2
		elif self.bc.left == 2: # Dirichlet
			b -= self.bc.lval*A[:,0]
			A[0,:] = 0
			A[:,0] = 0
			A[0,0] = 1
			b[0] = self.bc.lval
		if self.bc.right == 0: # Neumann
			b[-1] += self.bc.rval
		elif self.bc.right == 1: # Robin
			b[-1] += 2*self.bc.rval
			A[-1,-1] += 1/2
		elif self.bc.right == 2: # Dirichlet
			b -= self.bc.rval*A[:,-1]
			A[-1, :] = 0
			A[ :,-1] = 0
			A[-1,-1] = 1
			b[-1] = self.bc.rval
		return (A, b, A_nobc, b_nobc)
	def newton_jacobian(self, u0):
		'''
		Compute the jacobian for a given value of u0
		'''
		J = np.zeros((self.nel+1,self.nel+1))
		for k in range(self.nel):
			x0 = self.x[k]
			x1 = self.x[k+1]
			jacob = (x1-x0)/2
			# Quadrature points mapped to the element
			x = (x1+x0)/2 + self.xq*jacob
			# u and du/dx values at quadrature points
			uq = np.dot(self.b, u0[k:k+self.qorder])
			dudxq = np.dot(self.dbdx, u0[k:k+self.qorder])
			# Sigma(x,u) values at quadrature points
			Sq = self.siga(uq, x)*jacob
			# D(x,u) values at quadrature points
			Dq = self.diff(uq, x)/jacob
			dDduq = self.dDdu(uq, x)
			jlet = np.zeros((self.qorder,self.qorder))
			for i in range(self.qorder):
				for j in range(self.qorder):
					jlet[j,i] = sum(self.wq*Sq*self.wq*self.b[:,j]) + \
							sum(self.wq*self.dbdx[:,i]*Dq*self.dbdx[:,j]) + \
							sum(self.wq*self.b[j]*dDduq*dudxq)
			J[k:k+self.qorder,k:k+self.qorder] += jlet
		# Apply Boundary Conditions
		# Neumann is the addition of a const, therefore it doesn't change the residual
		if self.bc.left == 1: # Robin
        		J[0,0] += 1/2
		elif self.bc.left == 2: # Dirichlet
			J[0,:] = 0
			J[:,0] = 0
			J[0,0] = 1
		if self.bc.right == 1: # Robin
        		J[-1,-1] += 1/2
		elif self.bc.right == 2: # Dirichlet
			J[-1, :] = 0
			J[ :,-1] = 0
			J[-1,-1] = 1
		return J
	def newton_residual(self, u0):
		'''
		Compute the residual for a given value of u0
		'''
		F = np.zeros(self.nel+1)
		for i in range(self.nel):
			x0 = self.x[i]
			x1 = self.x[i+1]
			jacob = (x1-x0)/2
			# Quadrature points mapped to the element
			x = (x1+x0)/2 + self.xq*jacob
			# u and du/dx values at quadrature points
			uq = np.dot(self.b, u0[i:i+self.qorder])
			dudxq = np.dot(self.dbdx, u0[i:i+self.qorder])
			# D(x,u) values at quadrature points
			Dq = self.diff(uq, x)/jacob
			# Sigma(x,u) values at quadrature points
			Sq = self.siga(uq, x)*jacob
			# q(x) at the quadrature points
			Qq = self.esrc(x)*jacob
			local = np.zeros(self.qorder)
			for j in range(self.qorder):
				local[j] = sum(uq*Sq*self.wq*self.b[:,j]) + np.dot(Dq*self.wq*self.dbdx[:,j], dudxq) - np.dot(Qq*self.wq, self.b[:,j])
			F[i:i+2] += local
		# Apply Boundary Conditions
		if self.bc.left == 0: # Neumann
			F[0] -= self.bc.lval
		elif self.bc.left == 1: # Robin
        		F[0] += 1/2 * u0[0] - 2*self.bc.lval
		elif self.bc.left == 2: # Dirichlet
			F[0] = u0[0] - self.bc.lval
		if self.bc.right == 0: # Neumann
			F[-1] -= self.bc.rval
		elif self.bc.right == 1: # Robin
        		F[0] += 1/2 * u0[-1] - 2*self.bc.rval
		elif self.bc.right == 2: # Dirichlet
			F[-1] = u0[-1] - self.bc.rval
		
		return F
	def newton_solve(self, ustart, tol=1E-8, max_iter=100):
		'''
		Uses newton's method to solve the system.
		Return (solution, converged status, number of iterations)
		'''
		k = 0
		resi = self.newton_residual(ustart)
		norm0 = norm(resi)
		uk = ustart.copy()
		if norm0 == 0:
			return (uk, True, 0)
		while norm(resi)/norm0 > tol and k < max_iter:
			k+=1
			jaco = self.newton_jacobian(uk)
			du = solve(jaco, -resi)
			uk += du
			resi = self.newton_residual(uk)
		return (uk, norm(resi)/norm0 <= tol, k)
	def qoi(self, res):
		'''
		Returns a forward, and adjoint Quantity of interest
		'''
		(u, converged, itr) = self.newton_solve(np.ones(self.nel+1))
		return (self.qoif(u, res), self.qoia(u, res))
	def qoif(self, u, res):
		'''
		Calculates <u,r> in the forward manner
		'''
		qoif = 0
		for i in range(self.nel):
			x0 = self.x[i]
			x1 = self.x[i+1]
			jacob = (x1-x0)/2
			x = (x1+x0)/2 + self.xq*jacob 
			qoif += jacob * sum(self.wq*np.dot(self.b, u[i:i+self.qorder])*res(x))
		return qoif
	def qoia(self, u, res):
		'''
		Calculates <u,r> in an adjoint manner. I.e. <u*,q>
		'''
		adjoint_sys = femsys(self.nel, self.width, self.diff, self.dDdu, self.siga, res, self.bc.astuple())
		[A,b,A_nbc,b_nbc] = self.assemble(u)
		[As,bs,As_nbc,bs_nbc] = adjoint_sys.assemble(u)
		us = solve(A, bs_nbc)
		qoia = sum(us*b) #+ sum(us*np.dot((As-A),u))
		return qoia
def TestFun():
	sigma = lambda u, x: 0
	diff = lambda u, x: 1
	ddiff = lambda u, x: 0
	q = lambda x: 0
	res = lambda x: 1
	fem = femsys(100, 1, diff, ddiff, sigma, q, (2,0,2,1))
	(qoif, qoia) = fem.qoi(res)
	print(qoif, qoia, abs(qoif-qoia))

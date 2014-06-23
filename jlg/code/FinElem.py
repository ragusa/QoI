#!/usr/bin/env python3

import scipy as sp
import numpy as np
from numpy.polynomial.legendre import leggauss as GL_quad
from numpy.linalg import solve
from scipy import sparse
from pylab import *
import matplotlib.pyplot as plt
from scipy.optimize import minimize
def assemble_system(par, data):
	nnx = (par.porder+1)*par.nel
	n = par.nel * par.porder + 1
	A = np.zeros((n,n))#sparse.csc_matrix((n,n))
	rhs = np.zeros(n)
	for iel in range(par.nel):
		x0 = par.x[iel]
		x1 = par.x[iel+1]
		J = (x1-x0)/2 # Jacobian
		first = par.gn[iel,0]
		last = par.gn[iel,-1]+1
		A[first:last,first:last] = A[first:last,first:last] + (data.siga * par.m * J + data.diff * par.k/J)
		if callable(data.esrc):
			rhs[par.gn[iel,:]] = rhs[par.gn[iel,:]] + [data.esrc(x0), data.esrc(x1)] * par.f * J
		elif isinstance(data.esrc, (int, long, float)):
			rhs[par.gn[iel,:]] = rhs[par.gn[iel,:]] + data.esrc * par.f * J
		elif isinstance(data.esrc, (list, np.ndarray)):
			rhs[par.gn[iel,:]] = rhs[par.gn[iel,:]] + data.esrc[par.gn[iel,:]] *par. f * J
	A_nobc = np.copy(A)
	rhs_nobc = np.copy(rhs)
	Dirichlet_nodes=[]
	Dirichlet_val=[]	
	if   data.bc.left.type == 0:
		rhs[0] += data.bc.left.C
	elif data.bc.left.type == 1:
		rhs[0] += 2*data.bc.left.C
		A[0,0] += 1/2
	elif data.bc.left.type == 2:
		Dirichlet_nodes.append(0)
		Dirichlet_val.append(data.bc.left.C)
	if   data.bc.rite.type == 0:
		rhs[-1] += data.bc.rite.C
	elif data.bc.rite.type == 1:
		rhs[-1] += 2*data.bc.rite.C
		A[-1,-1] += 1/2
	elif data.bc.rite.type == 2:
		Dirichlet_nodes.append(-1)
		Dirichlet_val.append(data.bc.rite.C)
	for i in range(len(Dirichlet_nodes)):
		id = Dirichlet_nodes[i]
		bcval = Dirichlet_val[i]
		rhs = rhs - bcval*A[:,id]
		A[id,:] = np.zeros(n)
		A[:,id] = np.zeros(n)
		A[id,id] = 1
		rhs[id] = bcval
	return [A, rhs, par, A_nobc, rhs_nobc]
def feshpln(xv, p):
	shapefun = np.zeros((len(xv),p+1))
	dhdx = np.zeros((len(xv),p+1))
	shapefun = np.array([[0.788675134594813, 0.211324865405187], [0.211324865405187, 0.788675134594813]])
	dhdx = np.array([[-.5, .5],[-.5,.5]])
	return [shapefun, dhdx]

class side:
	def __init__(self, t, C):
		#0=neumann, 1=robin, 2=dirichlet
		self.type = t
		self.C = C
class boundary:
	def __init__(self, bc):
		self.left = side(bc[0], bc[1])
		self.rite = side(bc[2], bc[3])
class dat:
	def __init__(self, d, sa, src, w, bc):
		self.diff = d
		self.siga = sa
		self.esrc = src
		self.width = w
		self.bc = boundary(bc)
class npar:
	def __init__(self, n, data):
		self.nel = n
		self.x = np.linspace(0, data.width, self.nel+1)
		self.porder = 1
		self.ndofs = self.porder*self.porder+1
		self.gn = np.int_(self.gen_globalnums(self.nel, self.porder))
		q_pts = 2;
		[xq, wq] = np.polynomial.legendre.leggauss(q_pts) #Quadrature
		self.m = np.zeros((self.porder+1,self.porder+1))
		self.k = np.copy(self.m)
		self.f = np.zeros((self.porder+1,1))
		[b, dbdx] = feshpln(xq, self.porder)
		for i in range(self.porder+1):
			for j in range(self.porder+1):
				self.m[i,j] = np.dot(wq * b[:,i], b[:,j])
				self.k[i,j] = np.dot(wq * dbdx[:,i], dbdx[:,j])
			self.f[i] = np.dot(wq,b[:,i])
	def gen_globalnums(self, nel, porder):
		gn = np.zeros((nel, porder+1))
		gn[0] = np.linspace(0,porder,porder+1)
		for i in range(1,nel):
			gn[i] = [gn[i-1,-1], gn[i-1,1]+porder]
		return gn
def QoI(nel, diff, siga, src, width, bc, res):
	(par, A, A_nbc, u, b, b_nbc) = fem(nel, diff, siga, src, width, bc)
	(pars, As, As_nbc, us, bs, bs_nbc) = fem(nel, diff, siga, res, width, (0,0,0,0))

	# Forward
	QoIf = 0;
	for i in range(par.nel):
		x0 = par.x[i]
		x1 = par.x[i+1]
		QoIf += ((x1-x0)/2) * sum([res(x0), res(x1)] * np.transpose(par.f) * u[i:i+2])
	QoIa = sum(us*b) + sum(us*np.dot((As-A),u))
	return(QoIf, QoIa, par.x, u, us)
def sensitivity(nel, diff, siga, src, width, bc, res):
	(par, A, A_nbc, u, b, b_nbc) = fem(nel, diff, siga, src, width, bc)
	(pars, As, As_nbc, us, bs, bs_nbc) = fem(nel, diff, siga, res, width, (0,0,0,0))
	deltaA = 1
	deltab = 0
	sensitivity = 0
	lin_interpol = lambda u0, u1, x: (u0-u1)/2*x+(u0+u1)/2
	[pts, wq] = np.polynomial.legendre.leggauss(2)
	for k in range(nel):
		x0 = par.x[k]
		x1 = par.x[k+1]
		h = x1 - x0
		x = lambda z: h*(z+1)/2 + x0
		u_f = lambda x: lin_interpol(u[k], u[k+1], x)
		u_a = lambda x: lin_interpol(us[k], us[k+1], x)
		sensitivity += deltab*sum(wq*h/2*u_a(pts)) - deltaA*sum(wq*h/2*u_f(pts)*u_a(pts))
	return sensitivity, par, u
def sensitivity_direct(nel, diff, siga, src, width, bc, res, perturb):
	(qoi1, x,x,x,x) = QoI(nel, diff, siga, src, width, bc, res)
	(qoi2, x,x,x,x) = QoI(nel, diff, siga+perturb, src, width, bc, res)
	return (qoi2-qoi1)/perturb
def fem(nel, diff, siga, src, width, bc):
	data = dat(diff, siga, src, width, bc)
	par = npar(nel, data)
	[A, b, par, A_nob, b_nob] = assemble_system(par, data)
	u = solve(A,b)
	return (par, A, A_nob, u, b, b_nob)
def residual(nel, width, bc, D, Dderiv, q, u):
	[xq, wq] = np.polynomial.legendre.leggauss(2)
	F = np.zeros(nel+1)
	for i in range(nel):
		# X values at the ends of the interval
		x0 = width/nel * i
		x1 = width/nel * (i+1)
		# Jacobian
		jacob = (x1-x0)/2
		# Basis function (and derivative) at the quad points
		b = np.array([(1-xq)/2, (1+xq)/2])
		delb = np.array([[-.5, .5],[-.5, .5]])
		# Quad points on the interval
		x = (x1+x0)/2 + xq*jacob
		# u values at quad points
		us = np.dot(b, u[i:i+2])
		# du/dx at the quad points
		delu = np.dot(delb, u[i:i+2])
		
		ds = D(us)/jacob #D at u values (at quadrature points)
		qs = q(x)*jacob # Source
		# Construct local residual
		fi = np.zeros(2)
		for j in range(2):
			fi[j] = np.dot(ds*wq*delb[:,j], delu) - np.dot(qs*wq, b[:,j])
		
		F[i:i+2] += fi
	# Dirchlet
	if bc[0] == 2: # Left
		F[0] = u[0] - bc[1]
	if bc[2] == 2: # Right
		F[-1] = u[-1] - bc[3]
	return F
def jacobian(nel, width, bc, D, Dderiv, q, u, unmod=False):
	[xq, wq] = np.polynomial.legendre.leggauss(2)
	J = np.zeros((nel+1, nel+1))
	for k in range(nel):
		# X values at the ends of the interval
		x0 = width/nel * k
		x1 = width/nel * (k+1)
		# Jacobian
		jacob = (width/nel)/2
		# Basis function (and derivative) at the quad points
		b = np.array([(1-xq)/2, (1+xq)/2])
		delb = np.array([[-.5, .5],[-.5, .5]])
		# Quad points on the interval
		x = (x1+x0)/2 + xq*jacob
		# u values at quad points
		us = np.dot(b, u[k:k+2])
		# du/dx at the quad points
		delu = np.dot(delb, u[k:k+2])
		ds = D(us)/jacob #D at u values (at quadrature points)
		dDdu = Dderiv(us)
		jlet = np.zeros((2,2))
		for i in range(2):
			for j in range(2):
				jlet[i,j] = sum(wq*delb[:,i]*ds*delb[:,j]) + sum(wq*b[j]*dDdu*delu)
		J[k:k+2,k:k+2] += jlet
	if unmod:
		return J
	# Dirchlet
	if bc[0] == 2: # Left
		J[0,:] = 0
		J[:,0] = 0
		J[0,0] = 1
	if bc[2] == 2: # Right
		J[-1,:] = 0
		J[:,-1] = 0
		J[-1,-1] = 1
	return J
def norm(x):
	return sqrt(sum(i**2 for i in x))
def fem_nonlinear(nel, width, bc, D, Dderiv, q, ustart):
	"""
	Performs Newton's method to solve -∇·D(u)∇u = q(x)
	"""
	converged = False
	F_init = norm(residual(nel, width, bc, D, Dderiv, q, ustart))
	tol = 1E-12
	unew = ustart
	k = 0
	deltau = 0
	while not converged and k < 1000:
		uold = unew
		F = residual(nel, width, bc, D, Dderiv, q, unew)
		# Convergence check
		if norm(F)/F_init < tol:
			converged = True
			break
		J = jacobian(nel, width, bc, D, Dderiv, q, unew)
		deltau = solve(J, -F)
		unew = uold + deltau
		#print(k, "F=", norm(F))
		#print(k, "F/F=", norm(F)/F_init)
		#print(k, "u=", unew)
		#print(k, "diff", (norm(unew-uold)/norm(uold)))
		k += 1
	return [unew, converged, k]
def nonlinearAssemble(nel, width, bc, D, q, u, unmod=False):
	nullfun = lambda x: 0
	# The Jacobian used in the newton's is the A matrix with an extra term
	# nullfun kills that term
	A = jacobian(nel, width, bc, D, nullfun, q, u, unmod)
	[xq, wq] = np.polynomial.legendre.leggauss(2)
	basis = np.array([(1-xq)/2, (1+xq)/2])
	b = np.zeros(nel+1)
	f = np.zeros(2)
	for i in range(2):
		f[i] = np.dot(wq,basis[:,i])
	for i in range(nel):
		x0 = width/nel * i
		x1 = width/nel * (i+1)
		jacob = (x1-x0)/2
		b[i:i+2] += [q(x0), q(x1)] * f * jacob
	if not unmod:
		#Dirchlet hardcoded
		b -= bc[1]*A[:,0]
		b[0] = bc[1]
		b -= bc[3]*A[:,-1]
		b[-1] = bc[3]
	return (A,b)
		
def nonlinearQoI(nel, width, bc, D, Dderiv, q, ustart, res):
	[u, isconverged, numiters] = fem_nonlinear(nel, width, bc, D, Dderiv, q, ustart)
	[xq, wq] = np.polynomial.legendre.leggauss(2)
	b = np.array([(1-xq)/2, (1+xq)/2])

	QoIf = 0
	for i in range(nel):
		x0 = width/nel * i
		x1 = width/nel * (i+1)
		jacob = (x1-x0)/2
		QoIf += jacob * sum(wq*np.dot(b,u[i:i+2])*[res(x0),res(x1)])
	
	(A, b) = nonlinearAssemble(nel, width, bc, D, q, u)
	(Astar, bstar) = nonlinearAssemble(nel, width, bc, D, res, u)
	ustar = solve(Astar,bstar)
	QoIa = sum(ustar*b) + sum(ustar*np.dot((Astar-A),u))
	
	return (QoIf, QoIa)

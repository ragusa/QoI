#!/usr/bin/env python3

import scipy as sp
import numpy as np
from numpy.polynomial.legendre import leggauss as GL_quad
from numpy.linalg import solve
from scipy import sparse
from pylab import *
import matplotlib.pyplot as plt
def assemble_system(par, data):
	nnx = (par.porder+1)*par.nel
	n = par.nel * par.porder + 1
	A = np.zeros((n,n))#sparse.csc_matrix((n,n))
	rhs = np.zeros(n)
	q_pts = 2;
	[xq, wq] = np.polynomial.legendre.leggauss(q_pts) #Quadrature
	m = np.zeros((par.porder+1,par.porder+1))
	k = np.copy(m)
	f = np.zeros((par.porder+1,1))
	[b, dbdx] = feshpln(xq, par.porder)
	for i in range(par.porder+1):
		for j in range(par.porder+1):
			m[i,j] = np.dot(wq * b[:,i], b[:,j])
			k[i,j] = np.dot(wq * dbdx[:,i], dbdx[:,j])
		f[i] = np.dot(wq,b[:,i])
	par.m = np.copy(m)
	par.k = np.copy(k)
	par.f = np.copy(f)
	for iel in range(par.nel):
		x0 = par.x[iel]
		x1 = par.x[iel+1]
		J = (x1-x0)/2 # Jacobian
		first = par.gn[iel,0]
		last = par.gn[iel,-1]+1
		A[first:last,first:last] = A[first:last,first:last] + (data.siga*m*J+data.diff*k/J)
		rhs[par.gn[iel,:]] = rhs[par.gn[iel,:]]+data.esrc*f*J
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
		self.gn = np.int_(self.gen_gn(self.nel, self.porder))
	def gen_gn(self, nel, porder):
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
	f = par.f
	for i in range(par.nel):
		x0 = par.x[i]
		x1 = par.x[i+1]
		QoIf += res*((x1-x0)/2)*sum(f*u[i:i+1])

	As = np.transpose(As)
	print(us*np.dot((As-A),u) - us*(b - b_nbc))
	QoIa =     sum(us*b) + sum(us*np.dot((As-A),u))
	QoIa_nbc = sum(us*b) + sum(us*np.dot((A-A_nbc),u))
	return(QoIf, QoIa, QoIa_nbc, par.x, u, us)
	
def fem(nel, diff, siga, src, width, bc):
	data = dat(diff, siga, src, width, bc)
	par = npar(nel, data)
	[A, b, par, A_nob, b_nob] = assemble_system(par, data)
	u = solve(A,b)
	return (par, A, A_nob, u, b, b_nob)

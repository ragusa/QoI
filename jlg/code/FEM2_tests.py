#!/usr/bin/env python3
from FEM2 import femsys
import numpy as np
from numpy.linalg import solve

def testFlatSolution():
	sigma = lambda u, x: 1
	diff = lambda u, x: 1
	ddiff = lambda u, x: 0
	q = lambda x: 1
	res = lambda x: 1
	n = 100
	h = 1
	fem = femsys(n,h,diff,ddiff,sigma,q,(0,0,0,0))
	(qoif, qoia) = fem.qoi(res)
	assert abs(qoif - qoia) < 1E-12

def testMatLabProb():
	sigma = lambda u, x: 3/10
	diff = lambda u, x: 1
	ddiff = lambda u, x: 0
	q = lambda x: 3
	res = lambda x: 1
	n = 1000
	h = 10
	fem = femsys(n,h,diff,ddiff,sigma,q,(1,0,1,0))
	(qoif, qoia) = fem.qoi(res)
	print(qoif)
	print(qoia)
	assert abs(qoif - qoia) < 1E-12

def testOneCellLinear():
	sigma = lambda u,x: 0
	diff = lambda u, x: 1
	ddiff = lambda u, x: 0
	q = lambda x: 0
	res = lambda x: 1
	n = 2
	h = 1
	fem = femsys(n,h,diff,ddiff,sigma,q,(2,0,2,1))
	[A, b, A_nbc, b_nbc] = fem.assemble(np.zeros(n+1))
	u = fem.newton_solve(np.zeros(n+1))[0]
	(qoif, qoia) = fem.qoi(res)
	print(qoif)
	print(qoia)
	print(u)
testOneCellLinear()

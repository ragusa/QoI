#!/usr/bin/env python3
from multiprocessing import Pool
from FinElem import fem, QoI
import numpy as np
from numpy import exp, sqrt
from pylab import *
import matplotlib.pyplot as plt
from numpy.linalg import solve
def analytical(diff, siga, src, width, bc):
	# Solution to -d/dx D d/dx u + siga u = q
	dirchlet = lambda x, C: [exp( x*sqrt(siga/diff)),
					 exp(-x*sqrt(siga/diff)),
					 C-src/siga]
	robin    = lambda x, C: [ exp( x*sqrt(siga/diff))/4 + np.sign(x-width/2) * sqrt(siga/diff)*exp( x*sqrt(siga/diff))/2,
					  exp(-x*sqrt(siga/diff))/4 + -np.sign(x-width/2) * sqrt(siga/diff)*exp(-x*sqrt(siga/diff))/2,
					 C-src/siga/4]
	neumann  = lambda x, C: [np.sign(x-width/2) * diff*sqrt(siga/diff)*exp( x*sqrt(siga/diff)),
					 -np.sign(x-width/2) * diff*sqrt(siga/diff)*exp(-x*sqrt(siga/diff)),
					 C]
	A = np.zeros((2,2))
	b = np.zeros(2)
	if   bc[0] == 0: # Left Neumann
		r = neumann(0, bc[1])
		A[0] = r[0:-1]
		b[0] = r[-1]
	elif bc[0] == 1: # Left Robin
		r = robin(0, bc[1])
		A[0] = r[0:-1]
		b[0] = r[-1]
	elif bc[0] == 2: # Left Dirichlet
		r = dirchlet(0, bc[1])
		A[0] = r[0:-1]
		b[0] = r[-1]
	if   bc[2] == 0: # Right Neumann
		r = neumann(width, bc[3])
		A[1] = r[0:-1]
		b[1] = r[-1]
	elif bc[2] == 1: # Right Robin
		r = robin(width, bc[3])
		A[1] = r[0:-1]
		b[1] = r[-1]
	elif bc[2] == 2: # Right Dirichlet
		r = dirchlet(width, bc[3])
		A[1] = r[0:-1]
		b[1] = r[-1]
	k = solve(A,b)
	# k[0] * exp(x * sqrt(siga/D)) + k[1] * exp(x*sqrt(siga/D)) + q/siga
	u = lambda x: k[0] * exp(x * sqrt(siga/diff)) + k[1] * exp(-x*sqrt(siga/diff)) + src/siga
	intu = lambda x: k[0] * exp(x * sqrt(siga/diff)) / (sqrt(siga/diff)) - k[1] * exp(-x*sqrt(siga/diff)) / (sqrt(siga/diff)) + x*src/siga
	return u, intu
def uh_interpol(u0, u1, x):
	return (u0-u1)/2*x+(u0+u1)/2
def L2(ue, par, uh, order):
	err = 0
	[pts, w] = np.polynomial.legendre.leggauss(order)
	w /= 2
	for k in range( len(par.x) - 1):
		x0 = par.x[k]
		x1 = par.x[k+1]
		h = x1 - x0
		x = lambda z: h*(z+1)/2 + x0
		#dx = h/2 dz
		for q in range(len(w)):
			err += w[q]*(h/2)*(ue(x(pts[q])) - uh_interpol(uh[k], uh[k+1], pts[q]))**2
			#print(uh_interpol(uh[1][k], uh[1][k+1], pts[q]), ue(x(pts[q])))
	return err
def ueuh(nel, diff, siga, src, width, bc, plot=False):
	(par, A, uh, b) = fem(nel, diff, siga, src, width, bc)
	(ue, iue) = analytical(diff, siga, src, width, bc)
	l2err = L2(ue, par, uh,2)
	if plot:
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.plot(uh[0],uh[1], 'r')
		x = np.linspace(0,width,nel)
		ax.plot(x,ue(x), 'b')
		ax.axis([0,10,0,10])
		plt.show()
	return l2err
def uplot(x, u, us=None):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(x, u, 'b')
	if us != None:
		ax.plot(x, us, 'r')
	ax.axis([0,x[-1],0,1.1*max(u)])
	plt.show()
	
def QoIexact(nel, diff, siga, src, width, bc, res):
	(ue, iue) = analytical(diff, siga, src, width, bc)
#	QoI = 0
#	for x in np.linspace(0,width,num=5000,endpoint=True):
#		QoI += ue(x)*(width/5000)
#	QoI *= res
#	return QoI
	return res * (iue(width) - iue(0))
# bc types: 0=Neumann 1=Robin 2=Dirichlet
# bc = (left type, left C, right type, right C)
diff	= 1
siga	= 3/10
src	= 3
width	= 10
res = 1/width
bcs = [(0,10,0,0),
	(0,10,1,0),
	(0,10,2,0),
	(1,10,0,0),
	(1,10,1,0),
	(1,10,2,0),
	(2,10,0,0),
	(2,10,1,0),
	(2,10,2,0)
	]

ns = [1, 10, 100, 1000]
for bc in bcs:
	print(bc, '\tQe       Qf       Qa       Qa_nbc')
	for nel in ns:
		'''
		print('%i'%nel, end='\t')
		for bc in bcs:
			t = ueuh(nel, diff, siga, src, width, bc)
			print('%.3E'%t, end='\t')
		print()
		'''
		#t = ueuh(nel, diff, siga, src, width, bc, plot=True)
		(Qf, Qa, Qa_nbc, x, u, us) = QoI(nel, diff, siga, src, width, bc, res)
		Qe = QoIexact(nel, diff, siga, src, width, bc, res)
		#uplot(x, u, us)
		print('%15i\t%8.2E %8.2E %8.2E %8.2E %8.2E %8.2E %8.2E'%(nel, Qe, Qf, Qa, Qa_nbc, abs((Qf - Qe)/Qe), abs((Qa - Qe)/Qe), abs((Qa_nbc - Qe)/Qe)))
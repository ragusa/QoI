#!/usr/bin/env python3
import scipy.sparse as sp
import scipy.integrate as spint
from scipy.sparse import identity
from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import inv
import numpy as np
from numpy.linalg import solve

def ODEsolver(A, b, r, u0, steps=1000000):
	'''
	Uses the scipy.integrate.odeint function (based on ODEpack) to solve
	du/dt = A(t) u + b(t)
	from 0 to 1 using 1E6 steps
	'''
	dudt = lambda u,t: A(t) * u + b(t)
	tEnd = 1
	timepoints = np.linspace(0,tEnd,steps)
	u = spint.odeint(dudt, u0, timepoints)
	R = r(timepoints)
	QoI = sum(u[:,0]*R)/(steps+1)
	return (QoI, u)


def CN(A, b, r, u0, steps=100):
	tEnd = 1
	tau = tEnd/steps
	timepoints = np.linspace(0,tEnd,steps+1)
	u = np.zeros(steps+1)
	u[0] = u0
	A_f = sp.lil_matrix((steps+1, steps+1))
	b_f = sp.lil_matrix((steps+1,1))
	for i in range(steps):
		t = timepoints[i]#i * tEnd/steps
		tplus = timepoints[i+1]#(i+1) * tEnd/steps
		u[i+1] = (1/2*(b(t) + b(tplus)) - (-1/tau - A(t)/2)*u[i]) / (1/tau - A(tplus)/2)
		A_f[i+1,i] = (-1/tau - A(t)/2)
		A_f[i+1,i+1] = (1/tau - A(tplus)/2)
		b_f[i+1,0] = 1/2*(b(t) + b(tplus))
	A_f[0,0] = 1
	b_f[0,0] = u0
	#QoI = sum((u[1::]+u[0:-1])/2*r(timepoints)) * tau
	return (A_f, b_f, u)

def CNadjoint(A, b, r, u0, uN, steps=100):
	tEnd = 1
	tau = tEnd/steps
	timepoints = np.linspace(0,tEnd,steps+1)
	us = np.zeros(steps+1)
	us[-1] = 0
	for i in range(steps,-1,-1):
		t = timepoints[i-1]
		us[i-1] = (r(t) - (-1/tau - A(t)/2)*us[i]) / (1/tau - A(t)/2)
	#QoI = sum( (us[0:-1]+us[1:])/2 *(b(timepoints[1::])+b(timepoints[0:-1]))/2) *tau
	#QoI -= (us[-1]*uN - us[0]*u0)
	return u

A = lambda t: 2
b = lambda t: t
r = lambda t: 1
u0 = 1

#(Q, u) = ODEsolver(A, b, r, u0)
trueQ = (5*np.exp(2)-9)/8
print("True QoI", trueQ)
print("10^x steps", "QoI" , "Log10(rel. err.)", sep='\t')
print('Order\tQoIf2\tQoIf\tdiff\tQoIa\tQoIa-f\terr forward')
for order in range(1,5):
	steps = int(10**order)
	t = np.linspace(0,1,steps+1)
	tau = 1/steps
	
	(A_f, b_f, u) = CN(A, b, r, u0, steps)
	us = CNadjoint(A, b, r, u0, u[-1], steps)
	
	R = sp.csc_matrix(np.ones((steps+1,1)))
	B = b(t)
	K = identity(steps+1)
	K = sp.csc_matrix(K)
	K[0,0] = .5
	K[-1,-1] = .5
	K *= tau
	QoIf = sum(u * (K*R))
	QoIf2 = sum((u[1::]+u[0:-1])/2*r(t)) * tau
	
	A_adj = inv(K) * A_f.transpose() * K
	us = spsolve(A_adj, R)
	Kb_f = K*b_f
	QoIa = float(us*Kb_f)
	relE = (QoIf-trueQ)/trueQ
	print(order, QoIf, abs(QoIf2-QoIf), QoIa, abs(QoIf-QoIa), np.log10(relE), sep='\t')
	
	#b_vect = (b(t[1::])+b(t[0:-1]))/2
	#usDOTb = sum(us*b_vect)
	#print(usDOTb)

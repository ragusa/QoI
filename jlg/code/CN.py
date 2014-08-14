#!/usr/bin/env python3
import scipy as sp
import scipy.integrate as spint
import numpy as np
def CN(A, b, r, u0, steps=100):
	tEnd = 1
	tau = tEnd/steps
	timepoints = np.linspace(0,tEnd,steps+1)
	u = np.zeros(steps+1)
	u[0] = u0
	for i in range(steps):
		t = i * tEnd/steps
		tplus = (i+1) * tEnd/steps
		u[i+1] = (1/2*(b(t) + b(tplus)) - (-1/tau - A(t)/2)*u[i]) / (1/tau - A(tplus)/2)
	QoI = sum(u*r(timepoints))/(steps+1)
	return (QoI, u)

def ODEsolver(A, b, r, u0):
	'''
	Uses the scipy.integrate.odeint function (based on ODEpack) to solve
	du/dt = A(t) u + b(t)
	from 0 to 1 using 1E6 steps
	'''
	dudt = lambda u,t: A(t) * u + b(t)
	tEnd = 1
	steps = 1000000
	timepoints = np.linspace(0,tEnd,steps+1)
	u = spint.odeint(dudt, u0, timepoints)
	QoI = float(sum(u*r(timepoints))/(steps+1))
	return (QoI, u)


A = lambda t: 1
b = lambda t: 1
r = lambda t: 1
u0 = 0
trueQ = 0.7182818284590452353602874713526624977572470936999595

(Q, u) = ODEsolver(A, b, r, u0)
print("True QoI", trueQ)
print("10^x steps", "QoI" , "Log10(rel. err.)", sep='\t')
for order in range(1,7):
	(CNQ, CNu) = CN(A, b, r, u0, int(10**order))
	relE = abs(CNQ-trueQ)/trueQ
	print(order, CNQ , np.log10(relE), sep='\t')
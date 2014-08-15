#!/usr/bin/env python3
import scipy as sp
import scipy.integrate as spint
import numpy as np
def CN(A, b, r, u0, steps=100, adjoint=False):
	tEnd = 1
	tau = tEnd/steps
	timepoints = np.linspace(0,tEnd,steps+1)
	u = np.zeros(steps+1)
	u[0] = u0
	if adjoint:
		u[0] = 0
	for i in range(steps):
		t = timepoints[i]#i * tEnd/steps
		tplus = timepoints[i+1]#(i+1) * tEnd/steps
		if adjoint:
			t = tEnd - t
			tplus = tEnd - tplus
		u[i+1] = (1/2*(b(t) + b(tplus)) - (-1/tau - A(t)/2)*u[i]) / (1/tau - A(tplus)/2)
	if adjoint:
		u = u[::-1]
	QoI = sum(u*r(timepoints))/(steps+1)
	if adjoint:
		QoI += u0* u[0]
	return (QoI, u)

def ODEsolver(A, b, r, u0, steps=1000000):
	'''
	Uses the scipy.integrate.odeint function (based on ODEpack) to solve
	du/dt = A(t) u + b(t)
	from 0 to 1 using 1E6 steps
	'''
	dudt = lambda u,t: A(t) * u + b(t)
	tEnd = 1
	timepoints = np.linspace(0,tEnd,steps+1)
	u = spint.odeint(dudt, u0, timepoints)
	QoI = float(sum(u*r(timepoints))/(steps+1))
	return (QoI, u)

def CNadjoint(A, b, r, u0, uN, steps=100):
	tEnd = 1
	tau = tEnd/steps
	timepoints = np.linspace(0,tEnd,steps+1)
	us = np.zeros(steps+1)
	us[-1] = 0
	for i in range(-1,-steps-1,-1):
		t = timepoints[i-1]#i * tEnd/steps
		us[i-1] = (r(t) - (-1/tau - A(t)/2)*us[i]) / (1/tau - A(t)/2)
	QoI = sum(us[0:-1]*b(timepoints[0:-1])+us[1:]*b(timepoints[1:]))*tau/2
	QoI -= (us[-1]*uN - us[0]*u0)
	return (QoI, u)


A = lambda t: 1
b = lambda t: 1
r = lambda t: 1
u0 = 1

(Q, u) = ODEsolver(A, b, r, u0)
trueQ = Q#0.7182818284590452353602874713526624977572470936999595
print("True QoI", trueQ)
print("10^x steps", "QoI" , "Log10(rel. err.)", sep='\t')
for order in range(1,7):
	(CNQ, CNu) = CN(A, b, r, u0, int(10**order))
	(CNQa, CNua) = CNadjoint(A, b, r, u0, CNu[-1], int(10**order))
	relE = (CNQ-trueQ)**2/trueQ
	print(order, CNQ, CNQa, abs(CNQ-CNQa)*int(10**order), np.log10(relE), sep='\t')

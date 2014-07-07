#!/usr/bin/env python3

import numpy as np

'''
Solve for int(ru) from 0 to T, using Heun's method 
'''

# QoI = int(r·u) = int(q·u*) - u·u*|_0^T

#  du/dt = q
# -du*/dt = r


r = lambda t, u: t if t < 0.5 else 0 # negative derivative of u*
q = lambda t, u: 1 # derivative function of u

iters = 100 # Number of steps to take
t = lambda i: i/(iters-1) if i >= 0 else (iters - i)/(iters - 1) # Convert an step number to a time between 0 and 1
h = t(1)-t(0) # Step Width

u0 = 0
u = np.array([None for x in range(iters)])
u[0] = u0

usN = 0
us = np.array([None for x in range(iters)])
us[-1] = usN

# Solve for u
for i in range(1,iters):
	u_tilde = u[i-1] + h*q(t(i-1),u[i-1])
	u[i] = u[i-1] + h/2 * (q(t(i-1),u[i-1])+q(t(i), u_tilde))
# Solve for u*
for i in range(iters-2, -1, -1):
	us_tilde = us[i+1] + - h*r(t(i+1), u[i+1])
	us[i] =  us[i+1] - h/2 * (-r(t(i+1),u[i+1]) - r(t(i), u_tilde))
ri = np.array([r(t(i), u[i]) for i in range(iters)])
qi = np.array([q(t(i),u[i]) for i in range(iters)])
# Solve forward QoI
qoif = sum(h * u * ri)
# Solve adjoint QoI
qoia = sum(h * us * qi) - (qi[-1]*usN - ri[0]*u0)
# Print Results
print("QoI Forward: ", qoif)
print("QoI Adjoint: ", qoia)
print("Absolute Difference: ", qoif - qoia)

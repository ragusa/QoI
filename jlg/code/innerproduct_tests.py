#!/usr/bin/env python3
import numpy as np
from numpy.linalg import solve
from FEM2 import femsys
sigma = lambda u, x: x
diff = lambda u, x: 2*x
ddiff = lambda u, x: 2
q = lambda x: 2*x
res = lambda x: 1
n = 100
h = 1

fem = femsys(n,h,diff,ddiff,sigma,q,(2,0,2,0))
adjoint_sys = femsys(n,h,diff,ddiff,sigma, res, fem.bc.astuple())
(u, converged, itr) = fem.newton_solve(np.ones(fem.nel+1))
[A,b,A_nbc,b_nbc] = fem.assemble(u)
[As,bs,As_nbc,bs_nbc] = adjoint_sys.assemble(u)
(us, converged, itr) = adjoint_sys.newton_solve(np.ones(fem.nel+1))

(qoif, qoia) = fem.qoi(res)
uDOTR = sum(u*bs)
usDOTQ = sum(us*b_nbc)

print(qoif, qoia, uDOTR, usDOTQ)

print('QoIf', qoif)
print('QoIa', qoia)
print('U·R', uDOTR)
print('U*·Q', usDOTQ)
print('|QoIf-U·R|', abs(qoif-uDOTR))
print('|QoIa-U*·Q|', abs(qoia-usDOTQ))
print('|QoIf-QoIa|', abs(qoia-qoif))
print('|QoIf-U*·Q|', abs(usDOTQ-qoif))
print('|QoIa-U·R|', abs(uDOTR-qoif))

#!/usr/bin/env python3

import time
import datetime
import numpy as np
import scipy as sp
import scipy.optimize as opt
from collections import namedtuple
import matplotlib.pyplot as plt
import sympy as sym
from scipy import sparse as sp
def mass_matrix(N):
  '''
  Creates a sparse diagonal matrix of the mass matrix.
  '''
  h = 1/(N)
  # Create the mass matrix
  d_main = 4/3 * np.ones((N+1)) # Diagonal Entries
  d_off  = 1/3 * np.ones((N+1)) # Off diagonals
  M = sp.spdiags([d_off, d_main, d_off], [-1,0,1], (N+1), (N+1), format='csc')
  M[0,0] /= 2
  M[-1,-1] /= 2
  B = sp.csc_matrix((2*(N+1),2*(N+1)))
  B[:N+1,:N+1] = M
  B[N+1:,N+1:] = M
  return B*h/2
def parseBC(BC):
  boundary = namedtuple('boundary', ['ltype', 'lvalue', 'rtype', 'rvalue'])
  parsed = boundary(BC[0], BC[1], BC[2], BC[3])
  return parsed
def quadmaker(N, p):
  """
  Return quadrature points, weights, basis functions, and basis derivatives
  for order p
  """
  #[xq, wq] = np.polynomial.legendre.leggauss(p)
  xq = np.array([-1,1])
  wq = np.array([1,1])
  xd = np.linspace(-1,1,p)
  b    = np.zeros((p,p))
  dbdx = np.zeros((p,p))
  gn = np.zeros((N+1,p),dtype='int')
  gn[0,:] = [i for i in range(p)]
  for i in range(1,N+1):
    gn[i,:] = gn[i-1]+1
  for i in range(p):
    NUM = 1
    DEN = 1
    for j in range(p):
      if j != i:
        NUM *= xq - xd[j]
        DEN *= xd[i] - xd[j]
    b[:,i] = NUM/DEN
  for i in range(p):
    SUM = 0
    DEN = 1
    for j in range(p):
      if j != i:
        NUM = 1
        for k in range(p):
          if k!=i and k!=j:
            NUM *= xq - xd[k]
        SUM += NUM
        DEN *= xd[i]-xd[j]
    dbdx[:,i] = SUM/DEN
  return xq,wq,b,dbdx,gn
  
class fem:
  def __init__(self, kE, kT, aE, aT, bE, bT, BCE, BCT, qE, qT, N, meth='CN', mini='krylov', verbose=False):
    """
      Solves (d/dt - ∇·kE∇) E + a(E) + b(T) = qE(T,E)
             (d/dt - ∇·kT∇) T + b(E) + a(T) = qT(T,E)
    """
    self.kE = kE
    self.kT = kT
    self.aE = aE
    self.aT = aT
    self.bE = bE
    self.bT = bT
    self.qE = qE
    self.qT = qT
    #self.mass = np.array(mass_matrix(N).todense())
    self.N = N
    self.dx = 1/N
    self.x = np.linspace(0,1,N+1)

    self.BC = (parseBC(BCE), parseBC(BCT))

    self.qorder = 2
    [self.xq, self.wq, self.b, self.dbdx, self.gn] = quadmaker(self.N, self.qorder)

    if meth == 'CN':
      self.TR_residual = self.CN_TR_residual
    elif meth == 'FE':
      raise NotImplementedError('Forward Euler not implemented')
    elif meth == 'BE':
      raise NotImplementedError('Backward Euler not implemented')
    else:
      raise IOError('Unknown method "%s"'%(meth))
    
    if mini == 'leastsqs':
      self.mini = opt.leastsq
    elif mini == 'krylov':
      self.mini = lambda f, u, args: opt.root(f, u, args=args, method='krylov')
    elif mini == 'hybr':
      self.mini = lambda f, u, args: opt.root(f, u, args=args, method='hybr')
    else:
      raise IOError('Unknown minimizing routine "%s"'%(mini))

    self.verbose = verbose
  def SS_residual(self, u, t, u_for=None, adj=1):
    F = np.zeros(2*(self.N+1))
    for i in range(self.N): # Parallel?
      gnE = np.array([i,i+1])
      gnT = gnE + (self.N+1)
      x0 = self.x[i]
      x1 = self.x[i+1]
      jacob = (x1-x0)/2
      xq = self.xq + sum(self.x[i:i+2])/2
      local_E    = np.dot(self.b, u[gnE])   # if u_for==None else np.dot(self.b, u_for[gnE])   
      local_T    = np.dot(self.b, u[gnT])   # if u_for==None else np.dot(self.b, u_for[gnT])   
      local_dEdx = np.dot(self.dbdx, u[gnE])# if u_for==None else np.dot(self.dbdx, u_for[gnE])
      local_dTdx = np.dot(self.dbdx, u[gnT])# if u_for==None else np.dot(self.dbdx, u_for[gnT])
      local_aE   = self.wq*self.aE(xq, t, local_E)#*jacob
      local_aT   = self.wq*self.aT(xq, t, local_T)#*jacob
      local_bE   = self.wq*self.bE(xq, t, local_T)#*jacob
      local_bT   = self.wq*self.bT(xq, t, local_E)#*jacob
      local_qE   = self.wq*self.qE(xq, t)         #*jacob
      local_qT   = self.wq*self.qT(xq, t)         #*jacob
      #print(local_aT, local_bT, - local_qT)
      #print(self.bT(xq[1], t, local_E), xq, local_E)
      #print(self.qT(xq[1], t), -xq[1]**2+1)
      #oia=input("Pause")
      local_F = np.zeros((2,self.qorder))
      #print(gnE)
      for j in range(self.qorder):
        local_F[0,j] =  sum((adj*local_aE + adj*local_bE - local_qE)*self.b[:,j]) + adj*np.dot(self.wq*local_dEdx, self.dbdx[:,j])/jacob
        local_F[1,j] =  sum((adj*local_aT + adj*local_bT - local_qT)*self.b[:,j]) + adj*np.dot(self.wq*local_dTdx, self.dbdx[:,j])/jacob
      F[gnE] += local_F[0,:]
      F[gnT] += local_F[1,:]
    F[1:N] *= 1/2
    F[N+2:-1] *= 1/2
    for physics in [0,1]:
      start = physics * (self.N+1)
      end = start + self.N
      if self.BC[physics].ltype == 0:
        F[start] -= self.BC[physics].lvalue(t)
      elif self.BC[physics].ltype == 1:
        F[start] += 1/2 * u[start] - 2*self.BC[physics].lvalue(t)
      elif self.BC[physics].ltype == 2:
        F[start]  = u[start] - self.BC[physics].lvalue(t)
      if self.BC[physics].rtype == 0:
        F[end] -= self.BC[physics].rvalue(t)
      elif self.BC[physics].rtype == 1:
        F[end] += 1/2 * u[end] - 2*self.BC[physics].rvalue(t)
      elif self.BC[physics].rtype == 2:
        F[end]  = u[end] - self.BC[physics].rvalue(t)
    return -F

  def CN_TR_residual(self, u1, u0, t, dt, old_res=None):
    if old_res == None:
      old_res = self.SS_residual(u0, t+dt, None)
    return (u1-u0) - 1/2 * dt * (old_res+self.SS_residual(u1, t, None))
  def CNadj_TR_residual(self, u1, u0, t, dt, old_res=None, u_for=None):
    if old_res == None:
      old_res = self.SS_residual(u0, t-dt, u_for)
    return (u1-u0) - 1/2 * dt * (old_res+self.SS_residual(u1, t, u_for))
  def CN_for(self, u0, dts):
    if self.verbose:
      start_time = time.time()
      print('Beginning solve')
    u = np.empty(((self.N+1)*2, len(dts)+1))
    t_cur = 0
    u[:,0] = u0
    for i in range(1, len(dts)+1):
      old_res = self.SS_residual(u[:,i-1], t_cur)
      t_cur += dts[i-1]
      if self.verbose:
        step_time = time.time()
      TR_Residual = lambda u1: (u1 - u[:,i-1]) - 1/2 * dts[i-1] * (self.SS_residual(u1, t_cur)+old_res)
      u[:,i] = self.mini(TR_Residual, u[:,i-1], args=()).x
      if self.verbose:
        print("  Old res", self.SS_residual(u[:,i-1], t_cur-dts[i-1]))
        print("  New res", self.SS_residual(u[:,i], t_cur))
        print("  New TR ", TR_Residual(u[:,i]))
        print("  u i-1", u[:,i-1])
        print("  u i", u[:,i])
        print('Step %i done. t_cur = %f. Step Time = %f seconds. Total Elapsed Time = %f seconds'%(i, t_cur, time.time()-step_time, time.time()-start_time))
    return u
  def CN_adj(self, u0, dts, u_for):
    if self.verbose:
      start_time = time.time()
      print('Beginning solve')
    u = np.empty(((self.N+1)*2, len(dts)+1))
    t_cur = sum(dts)
    init_res = lambda u1: dts[-1]/2*self.SS_residual(u1, t_cur, u_for[:,-1])
    u[:,-1] = self.mini(init_res, u0, args=()).x
    for i in range(len(dts)-1,-1,-1):
      if self.verbose:
        step_time = time.time()
      adjmod = 1
      if i==len(dts)-1:
        adjmod *= 1/2
      if i==0:
        adjmod *=2
      old_res = self.SS_residual(u[:,i+1], t_cur, u_for[:,i+1], adj=adjmod)
      t_cur -= dts[i]
      print("ADJMOD", adjmod)
      TR_Residual = lambda u1: (u1 - adjmod*u[:,i+1]) - dts[i]/2 * (old_res + self.SS_residual(u1, t_cur, u_for[:,i]))
      u[:,i] = self.mini(TR_Residual, u[:,i+1], args=()).x
      
      if self.verbose:
        print("  Old res", self.SS_residual(u[:,i+1], t_cur-dts[i]))
        print("  New res", self.SS_residual(u[:,i], t_cur))
        print("  u i+1", u[:,i+1])
        print("  u i", u[:,i])
        print('Step %i done. t_cur = %f. Step Time = %f seconds. Total Elapsed Time = %f seconds'%(i, t_cur, time.time()-step_time, time.time()-start_time))
    return u
def innerproduct(u, res, N, T):
  tau = 1/T
  [xq, wq] = np.polynomial.legendre.leggauss(2)
  b = np.array([(1-xq)/2, (1+xq)/2])
  qoif_t = np.zeros(T+1)
  for t in range(T+1):
    ut = u[:,t]
    for i in range(N-1):
      x0 = i/(N-1)
      x1 = (i+1)/(N-1)
      jacob = (x1-x0)/2
      x = (x1+x0)/2 + xq*jacob
      qoif_t[t] += jacob*sum(wq*np.dot(b, ut[i:i+2])*res(x,t))
  return tau/2*sum(qoif_t[1::] + qoif_t[0:-1])
def TS_integral(a, b):
  N = len(a)
  T = len(a[0])
  tau = 1/(T-1)
  s = np.array([spatial_integral(a[:,i], b[:,i]) for i in range(T)])
  #print("S", s)
  return time_integral(s)
def spatial_integral(a, b):
  N = len(a)
  [xq, wq] = np.polynomial.legendre.leggauss(2)
  c = np.array([(1-xq)/2, (1+xq)/2])
  ans = 0
  for i in range(N-1):
    x0 = i/(N-1)
    x1 = (i+1)/(N-1)
    jacob = (x1-x0)/2
    x = (x1+x0)/2 + xq*jacob
    ans += jacob*sum(wq*np.dot(c, a[i:i+2])*np.dot(c, b[i:i+2]))
  return ans
def time_integral(b):
  T = len(b)
  tau = 1/(T-1)
  return tau/2*sum(b[1:] + b[0:-1])

def plotter(u):
  (N, T) = (u.shape)
  N = N//2 - 1
  fig = plt.figure()
  ax1 = fig.add_subplot(121)
  ax1.set_title=('T_r')
  ax2 = fig.add_subplot(122)
  ax2.set_title=('T_mat')
  for t in range(0,T,max(1, T//10)):
    ax1.plot(u[:N+1,t], 'o-')
    ax2.plot(u[N+1:,t], 'o-')
  plt.show()
  
kE = lambda x, t, E, T: 1
kT = lambda x, t, E, T: 1
aE = lambda x,t,E:  E   
aT = lambda x,t,T:  T 
bE = lambda x,t,T: -2*T 
bT = lambda x,t,E: -E   


from sympy.parsing.sympy_parser import parse_expr
from sympy import diff
from sympy import tanh
x, t = sym.symbols("x t")
E = x# 1/2 * (1-tanh(x-5*t))
T = 1# 1/2 * (1-tanh(x-5*t))


fE_expr = diff(E,t) - diff(diff(E,x),x) + E - 2*T
fT_expr = diff(T,t) - diff(diff(T,x),x) + T - E
E_fun = sym.lambdify((x,t), E, "numpy")
T_fun = sym.lambdify((x,t), T, "numpy")
qE = sym.lambdify((x,t), fE_expr, "numpy")
qT = sym.lambdify((x,t), fT_expr, "numpy")

gradE = sym.lambdify((x,t), diff(E, x), "numpy")
gradT = sym.lambdify((x,t), diff(T, x), "numpy")

BCE = (0, lambda t: -gradE(0,t), 0, lambda t: gradE(1,t))
BCT = (0, lambda t: -gradT(0,t), 0, lambda t: gradT(1,t))

print(gradT(1,1))
print(gradT(0,1))

print("QE = ", fE_expr)
print("QT = ", fT_expr)

rE = lambda x,t: 1
rT = lambda x,t: 1

N = 2

sys = fem(kE, kT, aE, aT, bE, bT, BCE, BCT, qE, qT, N, verbose=False)

Tmax = 1
Tsteps = 10
dts = (Tmax/Tsteps) * np.ones(Tsteps)
u0 = [E_fun(i/N,0) if i<N+1 else T_fun((i-N-1)/N, 0) for i in range(N*2+2)]

u_for = sys.CN_for(u0, dts)

qoi_f_E = innerproduct(u_for[:N+1], rE, N+1, Tsteps)
qoi_f_T = innerproduct(u_for[N+1:], rT, N+1, Tsteps)


print("QoI Forward, E ", qoi_f_E)
print("QoI Forward, T ", qoi_f_T)
#
#aBCE = (0, 0, 0, 0)
#aBCT = (0, 0, 0, 0)
#asys = fem(kE, kT, aT, aE, bT, bE, aBCE, aBCT, rE, rT, N, verbose=True)
#u_adj = asys.CN_adj(np.ones(2*N+2), dts, u_for)
#
#qEArr = np.array([[qE(i/N, Tmax*t/Tsteps) for t in range(Tsteps+1)] for i in range(N+1)])
#qTArr = np.array([[qT(i/N, Tmax*t/Tsteps) for t in range(Tsteps+1)] for i in range(N+1)])
#WE  = time_integral(u_adj[N,:]  * BCE[3] - u_adj[0,:]   * BCE[1])# Spatial Concomittants
#WT  = time_integral(u_adj[-1,:] * BCT[3] - u_adj[N+1,:] * BCT[1])# Spatial Concomittants
#TWE = (spatial_integral(u_for[:N+1,-1], u_adj[:N+1,-1]) - spatial_integral(u_for[:N+1, 0], u_adj[:N+1, 0]))        # Time concomittants
#TWT = (spatial_integral(u_for[N+1:,-1], u_adj[N+1:,-1]) - spatial_integral(u_for[N+1:, 0], u_adj[N+1:, 0]))        # Time concomittants
#TSE = TS_integral(u_adj[:N+1,:], qEArr) # Time-Space integral
#TST = TS_integral(u_adj[N+1:,:], qTArr) # Time-Space integral
##
###print("qE=\n", qEArr)
###print("qT=\n", qTArr)
###print("uA=\n", u_adj)
#print("W", WE)
#print("T", TWE)
#print("<>", TSE)
#print("W", WT)
#print("T", TWT)
#print("<>", TST)
#
#qoi_a_E = TSE - TWE + WE
#qoi_a_T = TST - TWT + WT
#
#print("QoI Adjoint, E ", qoi_a_E)
#print("QoI Adjoint, T ", qoi_a_T)
#
plotter(u_for)
#plotter(u_adj)


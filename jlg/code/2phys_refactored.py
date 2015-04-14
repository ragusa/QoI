#!/usr/bin/env python3

import time
import datetime
import numpy as np
import scipy as sp
import scipy.optimize as opt
from collections import namedtuple
import matplotlib.pyplot as plt
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
class fem_system:
  def __init__(self, N, phBC, thBC, meth='CN', mini='leastsqs', verbose='False'):
    """
    Notes:
    * Treats all flux arrays as [E, T]
    """
    self.N = N
    self.dx = 1/N
    self.x = np.linspace(0,1,N+1)
    
    self.limiter = 1 # for diffusion limiting
    self.m = 1 # for diffusion limiting

    self.z = 1 # z in sigma_a
    self.k = 0 # material thermal conductivity

    self.BC = (parseBC(phBC), parseBC(thBC))

    self.qorder = 2
    [self.xq, self.wq, self.b, self.dbdx, self.gn] = quadmaker(self.N, self.qorder)

    if meth == 'CN':
      self.tr_residual = self.CN_tr_residual
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
  def CN_tr_residual(self, u1, u0, dt, old_res):
    new_res = self.ss_residual(u1)
    #if self.verbose:
      #print(new_res)
    return (u1-u0) - 1/2 * dt * (old_res+self.ss_residual(u1))
  def ss_residual(self, u):
    """
    Calculate the steady state residual
    """
    F = np.zeros(2*(self.N+1))
    for i in range(self.N):
      gnE = np.array([i,i+1])
      gnT = gnE + (self.N+1)
      x0 = self.x[i]
      x1 = self.x[i+1]
      jacob = (x1-x0)/2

      local_E    = np.dot(self.b, u[gnE])
      local_T    = np.dot(self.b, u[gnT])
      local_dEdx = np.dot(self.dbdx, u[gnE])
      local_dTdx = np.dot(self.dbdx, u[gnT])
      DT = self.k * np.power(local_T, 2.5)
      SIGMA = self.z / np.power(local_T, 3)
      E = sum(self.wq*local_E) / sum(self.wq)
      grad_E = np.absolute(local_dEdx) / jacob
      DE = np.power(np.power(3*SIGMA, self.m) + np.power(self.limiter*grad_E/E, self.m), -1/self.m)

      local_F = np.zeros((2,self.qorder))
      aux1 = SIGMA * self.wq * (np.power(local_T, 4) - local_E)
      aux2 =    DT * self.wq * (local_dTdx)
      aux3 =    DE * self.wq * (local_dEdx)
      for j in range(self.qorder):
        local_F[0,j] = -np.dot(aux1, self.b[:,j])*jacob + np.dot(aux3, self.dbdx[:,j])/jacob
        local_F[1,j] =  np.dot(aux1, self.b[:,j])*jacob + np.dot(aux2, self.dbdx[:,j])/jacob
      F[gnE] += local_F[0,:]
      F[gnT] += local_F[1,:]
    for physics in [0,1]:
      start = physics * (self.N+1)
      end = start + self.N
      if self.BC[physics].ltype == 0:
        F[start] -= self.BC[physics].lvalue
      elif self.BC[physics].ltype == 1:
        F[start] += 1/2 * u[start] - 2*self.BC[physics].lvalue
      elif self.BC[physics].ltype == 2:
        F[start]  = u[start] - self.BC[physics].lvalue
      if self.BC[physics].rtype == 0:
        F[end] -= self.BC[physics].rvalue
      elif self.BC[physics].rtype == 1:
        F[end] += 1/2 * u[end] - 2*self.BC[physics].rvalue
      elif self.BC[physics].rtype == 2:
        F[end]  = u[end] - self.BC[physics].rvalue
    return -F
  def CN_solve(self, u0, dts):
    if self.verbose:
      start_time = time.time()
      print('Beginning solve')
    u = np.empty(((self.N+1)*2, len(dts)+1))
    u[:,0] = u0
    t_cur = 0
    for i in range(1, len(dts)+1):
      if self.verbose:
        step_time = time.time()
      old_res = self.ss_residual(u[:,i-1])
      u[:,i] = self.mini(self.tr_residual, u[:,i-1], args=(u[:,i-1], dts[i-1], old_res)).x
      t_cur += dts[i-1]
      if self.verbose:
        print('Step %i done. t_cur = %f. Step Time = %f seconds. Total Elapsed Time = %f seconds'%(i, t_cur, time.time()-step_time, time.time()-start_time))
    return u
def plotter(u):
  (N, T) = (u.shape)
  N = N//2 - 1
  fig = plt.figure()
  ax1 = fig.add_subplot(121)
  ax1.set_title=('T_r')
  ax2 = fig.add_subplot(122)
  ax2.set_title=('T_mat')
  for t in range(0,T,max(1, T//10)):
    ax1.plot(u[:N+1,t]**0.25, 'o-')
    ax2.plot(u[N+1:,t], 'o-')
  plt.show()
def test():
  N = 100
  Tmax = 3
  Tsteps = 100
  FEM = fem_system(N, (1,1,1,0), (0,0,0,0), mini='krylov', verbose=True)
  
  u_init = 1E-5 * np.ones(2*N+2)
  u_init[N+1:] = (1E-5)**(1/4)
  
  dts = (Tmax/Tsteps) * np.ones(Tsteps)
  
  u = FEM.CN_solve(u_init, dts)
  filename = 'P1 %s'%(str(datetime.datetime.now()).split(' ')[0])
  np.save(filename, u)
  plotter(u)
def load_and_plot(fn):
  u = np.load(fn)
  plotter(u)
print(":-) Hello World")
test()
#load_and_plot("2nd 2PHYSICS RUN 2015-01-13.npy")
print(":-( Goodbye World")

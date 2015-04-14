import numpy as np
from numpy import vectorize
from numpy.linalg import norm
import scipy as sp
import scipy.sparse as spar
from scipy.sparse.linalg import spsolve
import scipy.optimize as opt
from scipy.integrate import dblquad
import sympy
import matplotlib.pyplot as plt
from copy import deepcopy

class boundary_condition:
		def __init__(self, BC):
			self.ltype = BC[0]
			self.rtype = BC[2]
			self.left  = BC[1] if callable(BC[1]) else lambda t: BC[1]
			self.right = BC[3] if callable(BC[3]) else lambda t: BC[3]
class FEM:
	def __init__(self, N, dt, k, z, E_BC, T_BC):
		self. N = N
		self.dx = 1/self.N
		self.dt = dt
		self.k = k
		self.z = z
		self.E_BC = boundary_condition(E_BC)
		self.T_BC = boundary_condition(T_BC)
		self.Mass = self.generate_mass_matrix()
		self.qorder = 2
		[self.xq, self.wq, self.b, self.dbdx, self.gn] = self.generate_quadrature()
	def generate_quadrature(self):
		
	def generate_mass_matrix(self):
		
	def residual(u, t):

	def CrankNicholson(T,dt):
		

The following files can be found in the heat conduction example directory:

A. First, a certain number of test-cases with analytical solutions (no FEM discretization)
%-----------------------------------------------------------------------------------------

1. hc_analytical:

  Solves the linear HC equation analytically (the source term is uniform; the exact T is a 2nd order polynomial!). 

  Solves the adjoint HC equaton analytically. The response function r is piece-wise constant. We distinguish 3 cases
	r is uniform (1 value), r is 0 in 1 region (out of 2), r is 0 in 2 regions (out of 3); thus phi can be
	piece-wise linear/quadratic.

  The QoI (integral over space) is computed used a high-accuracy quadrature.

  We compute the sensitivities by solving the perturbed forward problem AND by using the adjoint relationships 
  and compare the 2 answers.
  The adjoint-based sensitivity is split into the various terms that make up the formula. In this simple problem, we have
  access to the perturbed forward solution, so we can assess the effects of the first-order adjoint formula and the exact one 
  (the former uses the unperturbed solution, the latter the perturbed solution).
  


2. hc_nonlinear_analytical:
  same as hc_analytical but for a nonlinear problem. the conductivity was chosen such that an analytical solution is
  available. not fully finished yet I believe


3. hc_fuel_nonlinear_analytical:
   work in progress. almost done.
  



B. Second, some 1D FEM code for heat conduction problems 
%-------------------------------------------------------

1. hc_ss_x_CG:
  Solves the linear HC equation, in steady state, ss, in x-slab geoemtry, using continuous FEM (Continuous Galerkin).
  Simulation parameters are stored in the numerical paramaters sturct, npar, and the data struct, dat.
  We can have different regions with different source functions and conductivity functions. 
  Each region is meshed (that is, 1 region .ne. one mesh cell).
  We use the inner products to compute the QoI int he forward and adjoitn manners.
  We also use the solution directly to verify the QoI and compute them "by hand".
  The same process applies to the sensitivities.
  Do not that the use of the inner products requires additional terms that come due to the Dirichlet bc, dot(r_functional_u,Tdiru)/


2. hc_ss_x_CG_nonlinear:
  Nonlinear version of hc_ss_x_CG. Noty finished yet, notably the Dirichlet bc stuff.
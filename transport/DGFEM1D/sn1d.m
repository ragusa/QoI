function sn1d
% Linear Discontinous FEM code for Sn transport in 1D
% Jean Ragusa, Ian Halvic
close all; clc; clear variables; clear global;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global npar dat snq IO_opts results
% set variable once for all
dat.pb_ID=24;
dat.forward_flux = true;
dat.adjoint_flux = ~dat.forward_flux;
do_transport_adjoint=false;
IO_opts.show_dqoi_pre_bc=false;
IO_opts.show_dqoi_from_bc=false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select angular approx (must be an even number)
sn=8;
n_moments=1; logi_galerkin=false;
[M,D,omega,weights] = gauss_legendre_1d(sn,n_moments,logi_galerkin);
snq.D = D; snq.M = M; snq.n_dir = sn;
snq.mu = omega; snq.w = weights;
%snq.mu = [-1,1];
snq.mu
snq.w
% sum of the weights
snq.sw = sum(snq.w);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select data problem
load_input(dat.pb_ID);
IO_opts.console_io = false;
dat.do_dsa = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[results.phi,results.E,results.Ebd,results.psi]=solveForwardSn;
[results.phia,results.Ea,results.Ebda,results.psia]=solveAdjointSn;
[results.phiVEF]=solveForwardVEF;
[results.phiVEFa]=solveAdjointVEF;

calcUnpertQOI
displayUnpertQOI


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Trying alternate method. VET of adjoint
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [phiVEF_alt]=solve_VEF_math_adjoint(dat.forward_flux,Ea,Ebda);
% do_plot(phiVEF_alt,'VEF-alt',0,dat.forward_flux)
% % solve forward VEF problem using IP
% [phiVEFa_alt]=solve_VEF(dat.adjoint_flux,Ea,Ebda);
% do_plot(phiVEFa_alt,'VEF-alt',100,dat.adjoint_flux)

%singleSolve
multiSolve
%multiSurfaceSolve

end
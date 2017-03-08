function sn1d
% Linear Discontinous FEM code for Sn transport in 1D
% Jean Ragusa,
close all; clc; % closes all plotting figures, clears console

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global npar dat snq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select angular approx (must be an even number)
sn=8;
n_moments=1;
logi_galerkin=false;
[M,D,omega,weights] = gauss_legendre_1d(sn,n_moments,logi_galerkin);
snq.D = D; snq.M = M; snq.n_dir = sn;
snq.mu = omega; snq.w = weights;
% sum of the weights
snq.sw = sum(snq.w);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select data problem
pb_ID=1;
load_input(pb_ID);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve forward transport problem using sweeps
forward = true;
do_dsa = true;
[phi,E,psi]=solve_transport(forward,do_dsa);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pretty plots
% do_plot(phi,E,psi)
do_plot(phi,E)

return
end
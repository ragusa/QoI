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
pb_ID=7;
load_input(pb_ID);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve forward transport problem using sweeps
forward = true;
do_dsa = true;
[phi,E,Ebd,psi]=solve_transport(forward,do_dsa);
Ebd
% pretty plots
do_plot(phi,0,E)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve forward VEF problem using IP
forward = true;
[phiVEF]=solve_VEF(forward,E,Ebd);
% pretty plots
do_plot(phiVEF,0)
do_plot((phiVEF./reshape(phi,npar.ndofs,1)-1),22)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve forward transport problem using sweeps
forward = false;
do_dsa = true;
[phia,Ea,Ebda,psia]=solve_transport(forward,do_dsa);
Ebda
% pretty plots
do_plot(phia,100,Ea)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve forward VEF problem using IP
forward = false;
[phiVEFa]=solve_VEF(forward,Ea,Ebda);
% pretty plots
do_plot(phiVEFa,100)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
forward = true;
qoi_sn_f = compute_qoi(~forward,phi);
fprintf('qoi using sn forward: \t %g \n',qoi_sn_f);
qoi_vef_f = compute_qoi(~forward,phiVEF);
fprintf('qoi using VEFforward: \t %g \n',qoi_vef_f);
forward=false;
qoi_sn_a = compute_qoi(~forward,phia);
fprintf('qoi using sn adjoint: \t %g \n',qoi_sn_a);
qoi_vef_a = compute_qoi(~forward,phiVEFa);
fprintf('qoi using VEFadjoint: \t %g \n',qoi_vef_a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve forward VEF problem using IP
forward = false;
[phiVEFmath]=solve_VEF_math_adjoint(forward,E,Ebd);
% pretty plots
do_plot(phiVEFmath,300)
qoi_vef_math_adj = compute_qoi(~forward,phiVEFmath);
fprintf('qoi using VEFmathadj: \t %g \n',qoi_vef_math_adj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


return
end
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

dat.bcincPert = 0.00;
dat.bcLPert = 0.00;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select data problem
pb_ID=14;
load_input(pb_ID);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
console_io = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve forward transport problem using sweeps
forward = true;
do_dsa = true;
[phi,E,Ebd,psi]=solve_transport(forward,do_dsa,console_io);
% pretty plots
do_plot(phi,0,forward,E)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve adjoint transport problem using sweeps
forward = false;
do_dsa = true;
[phia,Ea,Ebda,psia]=solve_transport(forward,do_dsa,console_io);
% pretty plots
do_plot(phia,100,forward,Ea)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

forward = true; sn=true;
qoi_sn_f = compute_qoi(~forward,phi,sn);
fprintf('qoi using sn forward: \t %g \n',qoi_sn_f);
%
forward=false;
qoi_sn_a = compute_qoi(~forward,phia,sn);
fprintf('qoi using sn adjoint: \t %g \n',qoi_sn_a);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % solve forward VEF problem using IP
% forward = true;
% [phiVEF]=solve_VEF(forward,E,Ebd);
% % pretty plots
% do_plot(phiVEF,0,forward)
% %do_plot((phiVEF./reshape(phi,npar.ndofs,1)-1),22)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % solve adjoint VEF problem using IP
% forward = false;
% [phiVEFa]=solve_VEF(forward,Ea,Ebda);
% % pretty plots
% do_plot(phiVEFa,100,forward)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% qoi_vef_f = compute_qoi(~forward,phiVEF,~sn);
% fprintf('qoi using VEFforward: \t %g \n',qoi_vef_f);
% %
% qoi_vef_a = compute_qoi(~forward,phiVEFa,~sn);
% fprintf('qoi using VEFadjoint: \t %g \n',qoi_vef_a);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % solve forward VEF problem using IP
% forward = false;
% [phiVEFmath]=solve_VEF_math_adjoint(forward,E,Ebd);
% % pretty plots
% do_plot(phiVEFmath,300,forward)
% qoi_vef_math_adj = compute_qoi(~forward,phiVEFmath,~sn);
% fprintf('qoi using VEFmathadj: \t %g \n',qoi_vef_math_adj);
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perturbations
dat.sourcePert = 0.0;
dat.sigsPert = +[0.1 0.1 0.1 0.1 0.1].*dat.sigs*0;
dat.sigtPert = -[0.1 0.1 0.1 0.1 0.1].*dat.sigt;
dat.sigaPert = dat.sigtPert - dat.sigsPert;
dat.psiIncPert = 0.0*dat.inc_forward;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute perturbed QoIs using adjoint method and unperturbed forward
forward=false;
qoi_sn_a = compute_perturbed_qoi_Sn(~forward,phia,phi,E,psi,psia,sn);
fprintf('perturbed qoi using sn adjoint: \t %g \n',qoi_sn_a);

error('ddd');

% qoi_vef_a = compute_perturbed_qoi_VEF(~forward,phiVEFa,phiVEF,E,Ebd);
% fprintf('perturbed qoi using VEFadjoint: \t %g \n',qoi_vef_a);
% qoi_vef_math_adj = compete_purturbed_qoi_VEF(~forward,phiVEFmath,phiVEF,E,Ebd);
% fprintf('perturbed qoi using VEFmathadj: \t %g \n',qoi_vef_math_adj);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Perturb data values for forward solve
dat.sigs = dat.sigs+dat.sigsPert;
dat.siga = dat.siga+dat.sigaPert;
dat.sigt = dat.sigt+dat.sigtPert;
dat.cdif = 1./(3*dat.sigt);
dat.qv_forward=(1+dat.sourcePert)*dat.qv_forward;
dat.inc_forward = dat.inc_forward + dat.psiIncPert ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load_input(13);
% solve forward transport problem using sweeps
forward = true;
do_dsa = true;
[phip,Ep,Ebdp,psip]=solve_transport(forward,do_dsa,console_io);
Ebdp;
% pretty plots
do_plot(phip,0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve forward VEF problem using IP
forward = true;
[phiVEFp]=solve_VEF(forward,Ep,Ebdp);
% pretty plots
%do_plot(phiVEFp,0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
forward = true;
qoi_sn_f = compute_qoi(~forward,phip);
fprintf('perturbed qoi using sn forward: \t %g \n',qoi_sn_f);
qoi_vef_f = compute_qoi(~forward,phiVEFp);
fprintf('perturbed qoi using VEFforward: \t %g \n',qoi_vef_f);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%Check 2nd order stuff
load_input(pb_ID);
forward=false;
qoi_sn_a_second = compute_perturbed_qoi_Sn(~forward,phia,phip,Ep,psip,psia);
fprintf('perturbed qoi using sn adjoint 2nd order: \t %g \n',qoi_sn_a_second);
qoi_vef_a = compute_perturbed_qoi_VEF(~forward,phiVEFa,phiVEFp,Ep,Ebd);
fprintf('perturbed qoi using VEFadjoint 2nd order: \t %g \n',qoi_vef_a);
qoi_vef_math_adj = compute_perturbed_qoi_VEF(~forward,phiVEFmath,phiVEFp,Ep,Ebd);
fprintf('perturbed qoi using VEFmathadj 2nd order: \t %g \n',qoi_vef_math_adj);
return
end
function sn1d
% Linear Discontinous FEM code for Sn transport in 1D
% Jean Ragusa,
close all; clc; clear variables;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global npar dat snq
% set variable once for all
forward_flux = true;
adjoint_flux = ~forward_flux;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select angular approx (must be an even number)
sn=8;
n_moments=1; logi_galerkin=false;
[M,D,omega,weights] = gauss_legendre_1d(sn,n_moments,logi_galerkin);
snq.D = D; snq.M = M; snq.n_dir = sn;
snq.mu = omega; snq.w = weights;
% sum of the weights
snq.sw = sum(snq.w);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select data problem
pb_ID=13;
load_input(pb_ID);
console_io = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve forward transport problem using sweeps
do_dsa = true;
[phi,E,Ebd,psi]=solve_transport(forward_flux,do_dsa,console_io);
do_plot(phi,0,forward_flux)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve adjoint transport problem using sweeps
do_dsa = true;
[phia,Ea,Ebda,psia]=solve_transport(adjoint_flux,do_dsa,console_io);
do_plot(phia,100,forward_flux)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sn=true;
qoi_sn_f = compute_qoi(forward_flux,phi,sn,[],[]);
fprintf('qoi using sn forward: \t %g \n',qoi_sn_f);
%
qoi_sn_a = compute_qoi(adjoint_flux,phia,sn,psi,psia);
fprintf('qoi using sn adjoint: \t %g \n',qoi_sn_a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve forward VEF problem using IP
[phiVEF]=solve_VEF(forward_flux,E,Ebd);
do_plot(phiVEF,0,forward_flux)
% figure(22); plot(phiVEF./reshape(phi,npar.ndofs,1)-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve adjoint VEF problem using IP
[phiVEFa]=solve_VEF(adjoint_flux,Ea,Ebda);
do_plot(phiVEFa,100,adjoint_flux)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qoi_vef_f = compute_qoi(forward_flux,phiVEF,~sn,[],[]);
fprintf('qoi using VEFforward: \t %g \n',qoi_vef_f);
qoi_vef_a = compute_qoi(adjoint_flux,phiVEFa,~sn,[],[]);
fprintf('qoi using VEFadjoint: \t %g \n',qoi_vef_a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve forward VEF problem using IP
[phiVEFmath]=solve_VEF_math_adjoint(adjoint_flux,E,Ebd);
do_plot(phiVEFmath,300,adjoint_flux)
qoi_vef_math_adj = compute_qoi(adjoint_flux,phiVEFmath,~sn,[],[]);
fprintf('qoi using VEFmathadj: \t %g \n',qoi_vef_math_adj);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perturbations
dat.sourcePert =[0 0 0 0 5];
dat.sigsPert = +[0.1 0.1 0.1 0.1 0.1].*dat.sigs*0;
dat.sigtPert = -[0.1 0.1 0.1 0.1 0.1]*0.*dat.sigt;
dat.sigaPert = dat.sigtPert - dat.sigsPert;
dat.psiIncPert = 0.0*dat.inc_forward;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perturbations
dat_saved = dat;
dat.qv_forward = dat.qv_forward + dat.sourcePert;
dat.sigs = dat.sigs + dat.sigsPert;
dat.sigt = dat.sigt + dat.sigtPert;
dat.siga = dat.siga + dat.sigaPert;
dat.inc_forward = dat.inc_forward + dat.psiIncPert;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve perturbed forward transport problem using sweeps
[phi_pert,E_pert,Ebd_pert,psi_pert]=solve_transport(forward_flux,do_dsa,console_io);
do_plot(phi_pert,10,forward_flux,E_pert)
% reset
dat = dat_saved;
qoi_sn_f_pert = compute_qoi(forward_flux,phi_pert,sn,[],[]);
fprintf('delta qoi using 2 sn forward runs: \t %g \n',qoi_sn_f_pert - qoi_sn_f);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute perturbed QoIs using adjoint method and unperturbed forward
delta_qoi_sn_a = compute_perturbed_qoi_Sn(adjoint_flux,phia,phi,psi,psia,sn);
fprintf('delta qoi using sn adjoint: \t\t %g \n',delta_qoi_sn_a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error('stopping here, not fully done with delta qoi. some of things below need to go');

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
% solve forward transport problem using sweeps
forward_flux = true;
do_dsa = true;
[phip,Ep,Ebdp,psip]=solve_transport(forward_flux,do_dsa,console_io);
Ebdp;
% pretty plots
do_plot(phip,0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve forward VEF problem using IP
forward_flux = true;
[phiVEFp]=solve_VEF(forward_flux,Ep,Ebdp);
% pretty plots
%do_plot(phiVEFp,0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
forward_flux = true;
qoi_sn_f = compute_qoi(~forward_flux,phip);
fprintf('perturbed qoi using sn forward: \t %g \n',qoi_sn_f);
qoi_vef_f = compute_qoi(~forward_flux,phiVEFp);
fprintf('perturbed qoi using VEFforward: \t %g \n',qoi_vef_f);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%Check 2nd order stuff
load_input(pb_ID);
forward_flux=false;
qoi_sn_a_second = compute_perturbed_qoi_Sn(~forward_flux,phia,phip,Ep,psip,psia);
fprintf('perturbed qoi using sn adjoint 2nd order: \t %g \n',qoi_sn_a_second);
qoi_vef_a = compute_perturbed_qoi_VEF(~forward_flux,phiVEFa,phiVEFp,Ep,Ebd);
fprintf('perturbed qoi using VEFadjoint 2nd order: \t %g \n',qoi_vef_a);
qoi_vef_math_adj = compute_perturbed_qoi_VEF(~forward_flux,phiVEFmath,phiVEFp,Ep,Ebd);
fprintf('perturbed qoi using VEFmathadj 2nd order: \t %g \n',qoi_vef_math_adj);
return
end
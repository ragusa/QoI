function sn1d
% Linear Discontinous FEM code for Sn transport in 1D
% Jean Ragusa,
close all; clc; clear variables; clear global;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global npar dat snq IO_opts
% set variable once for all
forward_flux = true;
adjoint_flux = ~forward_flux;
do_transport_adjoint=false;
IO_opts.show_dqoi_pre_bc=true;
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
pb_ID=25;
load_input(pb_ID);
console_io = false;
do_dsa = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve forward transport problem using sweeps
[phi,E,Ebd,psi]=solve_transport(forward_flux,do_dsa,console_io);
do_plot(phi,'Sn',0,forward_flux)
do_plot(E,'E',50,forward_flux,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve adjoint transport problem using sweeps
[phia,Ea,Ebda,psia]=solve_transport(adjoint_flux,do_dsa,console_io);
do_plot(phia,'Sn',100,adjoint_flux)
% do_plot(Ea,'Ea',50,forward_flux,true)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qoi_sn_f = compute_qoi(forward_flux,phi,sn,[],[]);
fprintf('qoi using sn forward: \t %g \n',qoi_sn_f);
%
qoi_sn_a = compute_qoi(adjoint_flux,phia,sn,psi,psia);
fprintf('qoi using sn adjoint: \t %g \n',qoi_sn_a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve forward VEF problem using IP
[phiVEF]=solve_VEF(forward_flux,E,Ebd);
do_plot(phiVEF,'VEF',0,forward_flux)
% figure(22); plot(phiVEF./reshape(phi,npar.ndofs,1)-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve forward VEF problem using IP
[phiVEFa_math]=solve_VEF_math_adjoint(adjoint_flux,E,Ebd);
do_plot(phiVEFa_math,'VEF',100,adjoint_flux)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if do_transport_adjoint
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % solve adjoint VEF problem using IP (Using adjoint E)
    [phiVEFa_trans_Ea]=solve_VEF(adjoint_flux,Ea,Ebda);
    do_plot(phiVEFa_trans_Ea,'VEF-transport-Ea',100,adjoint_flux)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % solve adjoint VEF problem using IP (Using forward E)
    [phiVEFa_trans_E]=solve_VEF(adjoint_flux,E,Ebda);
    do_plot(phiVEFa_trans_E,'VEF-transport-E',100,adjoint_flux)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qoi_vef_f = compute_qoi(forward_flux,phiVEF,~sn,[],[]);
fprintf('qoi using VEFforward: \t %g \n',qoi_vef_f);
%
qoi_vef_math_adj = compute_qoi(adjoint_flux,phiVEFa_math,~sn,[],[]);
fprintf('qoi using VEF math adj: \t %g \n',qoi_vef_math_adj);
%
if do_transport_adjoint 
    qoi_vef_a = compute_qoi(adjoint_flux,phiVEFa_trans_Ea,~sn,[],[]);
    fprintf('qoi using VEF transport adjoint (Ea): \t %g \n',qoi_vef_a);
    %
    qoi_vef_a_E = compute_qoi(adjoint_flux,phiVEFa_trans_E,~sn,[],[]);
    fprintf('qoi using VEF transport adjoint (E): \t %g \n',qoi_vef_a_E);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('-----BEGIN PERTURBED DATA OUTPUT----- \n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Perturbations. Used in adjoint sensitivity
dat.sigtPert = dat.sigtPertRegion.*dat.sigt*0.0+dat.sigtPertRegion.*0;
dat.sigsPert = dat.sigsPertRegion.*dat.sigs*-0.1+dat.sigsPertRegion.*0;
dat.sigaPert = dat.sigtPert - dat.sigsPert;
dat.sourcePert =dat.sourcePertRegion.*dat.qv_forward*0;
dat.psiIncPert = dat.incPertRegion.*dat.inc_forward*0+dat.incPertRegion*0.0;
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
do_plot(phi_pert,'Sn-pert',0,forward_flux)
do_plot(E_pert,'Epert',50,forward_flux,1)
do_plot(E_pert-E,'Epert-E',51,forward_flux,2)
%
qoi_sn_f_pert = compute_qoi(forward_flux,phi_pert,sn,[],[]);
fprintf('delta qoi using 2 sn forward runs: \t %g \n',qoi_sn_f_pert - qoi_sn_f);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve forward VEF problem using IP
[phiVEF_pert]=solve_VEF(forward_flux,E,Ebd);
do_plot(phiVEF_pert,'VEF-pert',0,forward_flux)
qoi_vef_f_pert = compute_qoi(forward_flux,phiVEF_pert,~sn,[],[]);
fprintf('delta qoi using 2 VEF forward runs: \t %g \n',qoi_vef_f_pert - qoi_vef_f);
% figure(22); plot(phiVEF./reshape(phi,npar.ndofs,1)-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute perturbed QoIs using Sn adjoint method and unperturbed forward
% reset
dat = dat_saved;
%
delta_qoi_sn_a = compute_perturbed_qoi_Sn(adjoint_flux,phia,phi,psi,psia,sn);
fprintf('delta qoi using sn adjoint: \t\t %g \n',delta_qoi_sn_a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute perturbed QoIs using VEF math adjoint method and unperturbed forward
delta_qoi_VEF_math_a = compute_perturbed_qoi_VEF(adjoint_flux,phiVEFa_math,phiVEF,E);
fprintf('delta qoi using VEF math adjoint: \t\t %g \n',delta_qoi_VEF_math_a);

if do_transport_adjoint
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute perturbed QoIs using VEF transport adjoint method (Ea) and unperturbed forward
    delta_qoi_VEF_trans_a = compute_perturbed_qoi_VEF(adjoint_flux,phiVEFa_trans_Ea,phiVEF,Ea);
    fprintf('delta qoi using VEF (Ea) transport adjoint: \t\t %g \n',delta_qoi_VEF_trans_a);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute perturbed QoIs using VEF transport adjoint method (E) and unperturbed forward
    delta_qoi_VEF_trans_a = compute_perturbed_qoi_VEF(adjoint_flux,phiVEFa_trans_E,phiVEF,E);
    fprintf('delta qoi using VEF (E) transport adjoint: \t\t %g \n',delta_qoi_VEF_trans_a);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

Rel_L1_diff=find_Eddington_diff(E,E_pert);
fprintf('relative L1 difference in E: \t\t %g \n',Rel_L1_diff);
%error('stopping here, not fully done with delta qoi.');

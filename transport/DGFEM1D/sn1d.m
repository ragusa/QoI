function sn1d
% Linear Discontinous FEM code for Sn transport in 1D
% Jean Ragusa, Ian Halvic
close all; clc; clear variables; clear global;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global npar dat snq IO_opts
% set variable once for all
forward_flux = true;
adjoint_flux = ~forward_flux;
do_transport_adjoint=false;
IO_opts.show_dqoi_pre_bc=false;
IO_opts.show_dqoi_from_bc=false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select angular approx (must be an even number)
sn=2;
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
pb_ID=57;
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
% solve forward VEF problem using IP
[phiVEF]=solve_VEF(forward_flux,E,Ebd);
do_plot(phiVEF,'VEF',0,forward_flux)
% figure(22); plot(phiVEF./reshape(phi,npar.ndofs,1)-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve forward VEF problem using IP
[phiVEFa]=solve_VEF_math_adjoint(adjoint_flux,E,Ebd);
do_plot(phiVEFa,'VEF',100,adjoint_flux)
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n-----BEGIN UNPERTURBED QOI DATA OUTPUT----- \n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qoi_sn_f = compute_qoi(forward_flux,phi,sn,psi,psia);
fprintf('qoi using sn forward: \t %g \n',qoi_sn_f);
%
qoi_sn_a = compute_qoi(adjoint_flux,phia,sn,psi,psia);
fprintf('qoi using sn adjoint: \t %g \n',qoi_sn_a);
%
qoi_vef_f = compute_qoi(forward_flux,phiVEF,~sn,[],[]);
fprintf('qoi using VEFforward: \t %g \n',qoi_vef_f);
%
qoi_vef_a = compute_qoi(adjoint_flux,phiVEFa,~sn,[],[]);
fprintf('qoi using VEF math adj: \t %g \n',qoi_vef_a);
%
if do_transport_adjoint 
    qoi_vef_a_trans = compute_qoi(adjoint_flux,phiVEFa_trans_Ea,~sn,[],[],0);
    fprintf('qoi using VEF transport adjoint (Ea): \t %g \n',qoi_vef_a_trans);
    %
    qoi_vef_a_trans_E = compute_qoi(adjoint_flux,phiVEFa_trans_E,~sn,[],[],0);
    fprintf('qoi using VEF transport adjoint (E): \t %g \n',qoi_vef_a_trans_E);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Trying alternate method. VET of adjoint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % solve forward transport problem using sweeps
% [phi_alt,E_alt,Ebd_alt,psi_alt]=solve_transport(forward_flux,do_dsa,console_io,1);
% do_plot(phi_alt,'Sn-alt',0,forward_flux)
% %do_plot(E_alt,'E',50,forward_flux,1)
% % solve adjoint transport problem using sweeps
% [phia_alt,Ea_alt,Ebda_alt,psia_alt]=solve_transport(adjoint_flux,do_dsa,console_io,1);
% do_plot(phia_alt,'Sn-alt',100,adjoint_flux)
% % do_plot(Ea,'Ea',50,forward_flux,true)
% solve forward VEF problem using IP
[phiVEF_alt]=solve_VEF_math_adjoint(forward_flux,Ea,Ebda);
do_plot(phiVEF_alt,'VEF-alt',0,forward_flux)
% solve forward VEF problem using IP
[phiVEFa_alt]=solve_VEF(adjoint_flux,Ea,Ebda);
do_plot(phiVEFa_alt,'VEF-alt',100,adjoint_flux)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n-----BEGIN PERTURBATION SENSITIVITY DATA OUTPUT ALT METHOD----- \n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% qoi_sn_f_alt = compute_qoi(forward_flux,phi_alt,sn,psi_alt,psia_alt,1);
% fprintf('qoi using sn forward: \t %g \n',qoi_sn_f_alt);
% %
% qoi_sn_a_alt = compute_qoi(adjoint_flux,phia_alt,sn,psi_alt,psia_alt,1);
% fprintf('qoi using sn adjoint-alt: \t %g \n',qoi_sn_a_alt);
% %
qoi_vef_f_alt = compute_qoi(forward_flux,phiVEF_alt,~sn,[],[]);
fprintf('qoi using VEFforward: \t %g \n',qoi_vef_f_alt);
%
qoi_vef_a_alt = compute_qoi(adjoint_flux,phiVEFa_alt,~sn,[],[]);
fprintf('qoi using VEF math adj: \t %g \n',qoi_vef_a_alt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n-----BEGIN PERTURBATION SENSITIVITY DATA OUTPUT----- \n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Perturbations. Used in adjoint sensitivity
dat.sigaPert = dat.sigaPertRegion.*dat.siga*0.0+dat.sigaPertRegion.*0;
dat.sigsPert = dat.sigsPertRegion.*dat.sigs*0.1+dat.sigsPertRegion.*0;
dat.sigtPert = dat.sigaPert + dat.sigsPert;
dat.sourcePert =dat.sourcePertRegion.*dat.qv_forward*0.0;
dat.psiIncPert = dat.incPertRegion.*dat.inc_forward*0.0+dat.incPertRegion*0.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perturbations
dat_saved = dat;
dat.qv_forward = dat.qv_forward + dat.sourcePert;
dat.sigs = dat.sigs + dat.sigsPert;
dat.sigt = dat.sigt + dat.sigtPert;
dat.siga = dat.siga + dat.sigaPert;
dat.cdif = 1./(3*dat.sigt);
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
fprintf('delta qoi using 2 sn forward runs: \t \t %g \n',qoi_sn_f_pert - qoi_sn_f);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve forward VEF problem using IP
[phiVEF_pert]=solve_VEF(forward_flux,E_pert,Ebd_pert);
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
fprintf('delta qoi using sn adjoint: \t\t\t %g \n',delta_qoi_sn_a);
%
delta_qoi_sn_a_pert = compute_perturbed_qoi_Sn(adjoint_flux,phia,phi_pert,psi_pert,psia,sn);
fprintf('delta qoi using sn adjoint (pert phi): \t\t\t %g \n',delta_qoi_sn_a_pert);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute perturbed QoIs using VEF math adjoint method and unperturbed forward
delta_qoi_VEF_a = compute_perturbed_qoi_VEF(adjoint_flux,phiVEFa,phiVEF,E);
fprintf('delta qoi using VEF math adjoint: \t\t %g \n',delta_qoi_VEF_a);
% Compute perturbed QoIs using VEF math adjoint method and unperturbed forward
delta_qoi_VEF_a_pert = compute_perturbed_qoi_VEF(adjoint_flux,phiVEFa,phiVEF_pert,E_pert);
fprintf('delta qoi using VEF math adjoint (phi pert): \t\t %g \n',delta_qoi_VEF_a_pert);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute perturbed QoIs using VEF math adjoint method and unperturbed
% forward from the Sn solve (which would be done to find E)
delta_qoi_VEF_SNphi = compute_perturbed_qoi_VEF(adjoint_flux,phiVEFa,phi,E);
fprintf('delta qoi using VEF adj w Sn fwd: \t\t %g \n',delta_qoi_VEF_SNphi);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Alternate Method
% Compute perturbed QoIs using VEF math adjoint method and unperturbed
% forward, alternate method
delta_qoi_VEF_a_alt = compute_perturbed_qoi_VEF(adjoint_flux,phiVEFa_alt,phiVEF_alt,Ea);
fprintf('delta qoi using VEF math adjoint alt: \t\t %g \n',delta_qoi_VEF_a_alt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n-----BEGIN DATA ON SENSITIVITY ERROR----- \n')
fprintf('This section contains data on the error in the sensitivity calculations. \n')
fprintf('In general, the "two forward solves" method will be used as the "true" value. \n')
fprintf('This means that there is a "true" sensitivity for two Sn solves and another \n')
fprintf('From two forward VEF solves. \n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('VEF adjoint error from SN forward: \t\t %g \n',(qoi_sn_f_pert - qoi_sn_f)-delta_qoi_VEF_a);
fprintf('VEF adjoint error from VEF forward: \t %g \n',(qoi_vef_f_pert - qoi_vef_f)-delta_qoi_VEF_a);
fprintf('VEF adjoint error from SN adjoint: \t\t %g \n',delta_qoi_sn_a-delta_qoi_VEF_a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n-----BEGIN DATA ON EDDINGTON PERTURBATION----- \n')
fprintf('This section contains data on the terms resulting from the \n')
fprintf('perturbed Eddington and boundary Eddington. \n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[deltaE_term,deltaB_L,deltaB_R,volValues]=compute_deltaE_QoI_term(phiVEFa,phiVEF,phiVEF_pert,E,E_pert,Ebd,Ebd_pert);
fprintf('Total delta E+B term: \t %g \n',deltaE_term+deltaB_L+deltaB_R);
fprintf('-deltaE term: \t\t\t %g \n',deltaE_term);
fprintf('-deltaB term: \t\t\t %g \n',deltaB_L+deltaB_R);
fprintf('--deltaB_L term: \t\t %g \n',deltaB_L);
fprintf('--deltaB_R term: \t\t %g \n',deltaB_R);
%Plots the contribution of the volumetric (delta E) term at each element.
%Crude plot, needs refined.
figure(70)
plot(volValues)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n-----BEGIN DATA ON OTHER/under development METRICS---- \n')
fprintf('This section contains data on other metrics. \n')
fprintf('A bit of a holding pen for exploratory values/development. \n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve forward VEF problem using IP
[phiVEFa_Epert]=solve_VEF_math_adjoint(adjoint_flux,E_pert,Ebd_pert);
do_plot(phiVEFa_Epert,'VEF Using E_{pert}',100,adjoint_flux)
do_plot(phiVEFa_Epert-phiVEFa,'VEF_{diff}',120,adjoint_flux)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rel_L1_diff=find_Eddington_diff(E,E_pert);
delta_qoi_VEF_math_a_Epert = compute_perturbed_qoi_VEF(adjoint_flux,phiVEFa_Epert,phiVEF_pert,E_pert);
fprintf('delta qoi using VEF math adjoint and Epert: \t\t %g \n',delta_qoi_VEF_math_a_Epert);
fprintf('Diff from using Epert: \t\t %g \n',delta_qoi_VEF_a-delta_qoi_VEF_math_a_Epert);
fprintf('relative L1 difference in E: \t\t %g \n',Rel_L1_diff);
%adjDiff = compute_qoi(adjoint_flux,phiVEFa_Epert,~sn,[],[])-compute_qoi(adjoint_flux,phiVEFa,~sn,[],[]);
%fprintf('NEED TO USE PERT FIX SOURCE USELESS. Diff from using two VEF adjoints: \t\t %g \n',adjDiff);


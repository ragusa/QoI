function sn1d_iterator
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
pb_ID=51;
load_input(pb_ID);
console_io = false;
do_dsa = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Specify Output File Path
outputMatrix=[];
today=datestr(now,'mmddyyHHMM');
%[pathstr,name,ext] = fileparts(pwd); 
file=['DGFEM1D_prob',int2str(pb_ID),'_',today,'.csv'];
filename=fullfile('C:\Users\Ian\checkout','output',file);
%filename=fullfile('E:\Storage\git_checkout','output',file);
% outputMatrix=['pid' 'QoISNf' 'QoISNa' 'QoIVEFf' 'QoIVEFa'];
% outputMatrix=[outputMatrix '%sigtPert' '%sigaPert' '%sourcePert' '%incPert'];
% outputMatrix=[outputMatrix 'dQoISNf' 'dQoIVEFf' 'dQoISNa' 'dQoIVEFa' 'Ediff'];
% Iterators for perturbation factors
sigtPertFactor=[0];
sigsPertFactor=linspace(-0.1,0.1,21);
sourcePertFactor=[0];
incPertFactor=[0];
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
% do_plot(Ea,'Ea',50,forward_flux,1)
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
[phiVEFa]=solve_VEF_math_adjoint(adjoint_flux,E,Ebd);
do_plot(phiVEFa,'VEF',100,adjoint_flux)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qoi_vef_f = compute_qoi(forward_flux,phiVEF,~sn,[],[]);
fprintf('qoi using VEFforward: \t %g \n',qoi_vef_f);
%
qoi_vef_a = compute_qoi(adjoint_flux,phiVEFa,~sn,[],[]);
fprintf('qoi using VEF math adj: \t %g \n',qoi_vef_a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('-----BEGIN PERTURBED DATA OUTPUT----- \n')
for ii=1:numel(sigtPertFactor)
    for jj=1:numel(sigsPertFactor)
        for kk=1:numel(sourcePertFactor)
            for ll=1:numel(incPertFactor)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Load Perturbations. Used in adjoint sensitivity
                    dat.sigtPert = dat.sigtPertRegion.*dat.sigt*sigtPertFactor(ii);
                    dat.sigsPert = dat.sigsPertRegion.*dat.sigs*sigsPertFactor(jj);
                    dat.sigaPert = dat.sigtPert - dat.sigsPert;
                    dat.sourcePert =dat.sourcePertRegion.*dat.qv_forward*sourcePertFactor(kk);
                    dat.psiIncPert = dat.incPertRegion.*dat.inc_forward*incPertFactor(ll);
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
                    %do_plot(phi_pert,'Sn-pert',0,forward_flux)
                    %do_plot(E_pert,'Epert',50,forward_flux,1)
                    %do_plot(E_pert-E,'Epert-E',51,forward_flux,2)
                    %
                    qoi_sn_f_pert = compute_qoi(forward_flux,phi_pert,sn,[],[]);
                    delta_qoi_sn_f=qoi_sn_f_pert - qoi_sn_f;
                    %fprintf('delta qoi using 2 sn forward runs: \t %g \n',delta_qoi_sn_f);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % solve forward VEF problem using IP
                    [phiVEF_pert]=solve_VEF(forward_flux,E_pert,Ebd_pert);
                    %do_plot(phiVEF_pert,'VEF-pert',0,forward_flux)
                    qoi_vef_f_pert = compute_qoi(forward_flux,phiVEF_pert,~sn,[],[]);
                    delta_qoi_VEF_f=qoi_vef_f_pert - qoi_vef_f;
                    %fprintf('delta qoi using 2 VEF forward runs: \t %g \n',delta_qoi_VEF_f);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Compute perturbed QoIs using Sn adjoint method and unperturbed forward
                    % reset
                    dat = dat_saved;
                    %
                    delta_qoi_sn_a = compute_perturbed_qoi_Sn(adjoint_flux,phia,phi,psi,psia,sn);
                    %fprintf('delta qoi using sn adjoint: \t\t %g \n',delta_qoi_sn_a);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Compute perturbed QoIs using VEF math adjoint method and unperturbed forward
                    delta_qoi_VEF_math_a = compute_perturbed_qoi_VEF(adjoint_flux,phiVEFa,phiVEF,E);
                    %fprintf('delta qoi using VEF math adjoint: \t\t %g \n',delta_qoi_VEF_math_a);

                    Rel_L1_diff=find_Eddington_diff(E,E_pert);
                    %fprintf('relative L1 difference in E: \t\t %g \n',Rel_L1_diff);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    [deltaE_term,deltaB_L,deltaB_R,volValues]=compute_deltaE_QoI_term(phiVEFa,phiVEF,phiVEF_pert,E,E_pert,Ebd,Ebd_pert);
                    total_deltaE_term=deltaE_term+deltaB_L+deltaB_R;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    outputLine=[pb_ID qoi_sn_f qoi_sn_a qoi_vef_f qoi_vef_a];
                    outputLine=[outputLine sigtPertFactor(ii) sigsPertFactor(jj) sourcePertFactor(kk) incPertFactor(ll) ];
                    outputLine=[outputLine  delta_qoi_sn_f  delta_qoi_VEF_f  delta_qoi_sn_a  delta_qoi_VEF_math_a Rel_L1_diff];
                    outputLine=[outputLine  total_deltaE_term];
                    outputMatrix = [outputMatrix; outputLine];
            end
        end
    end
end

csvwrite(filename,outputMatrix)
disp(['output stored in ',filename]);
%error('stopping here, not fully done with delta qoi.');

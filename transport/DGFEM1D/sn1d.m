function sn1d
% Linear Discontinous FEM code for Sn transport in 1D
% Jean Ragusa, Ian Halvic
close all; clc; clear variables; clear global;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global npar dat snq IO_opts results
% set variable once for all
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
dat.pb_ID=108;
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

singleSolve
%multiSolve

end

function singleSolve
global npar dat snq IO_opts results
% Load Perturbations. Used in adjoint sensitivity
dat.sigaPert = dat.sigaPertRegion.*dat.siga*-0.0+dat.sigaPertRegion.*0;
dat.sigsPert = dat.sigsPertRegion.*dat.sigs*0.0+dat.sigsPertRegion.*0;
dat.sigtPert = dat.sigaPert + dat.sigsPert;
dat.sourcePert =dat.sourcePertRegion.*dat.qv_forward*0.1;
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
% solve perturbed forward transport problem using sweeps
[results.phi_pert,results.E_pert,results.Ebd_pert,results.psi_pert]=solveForwardSn;
[results.phiVEF_pert]=solveForwardVEF;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
calcForwardSensitivity
dat = dat_saved;
calcAdjointSensitivity
displaySensitivity
generatePlots
displaySensitivityError
end

function multiSolve
global npar dat snq IO_opts results
fprintf('\n-----BEGIN MULTIPLE SENSITIVITY SOLVE----- \n')
%Specify Output File Path
outputMatrix=[];
today=datestr(now,'mmddyyHHMM');
%[pathstr,name,ext] = fileparts(pwd); 
file=['DGFEM1D_prob',int2str(dat.pb_ID),'_',today,'.csv'];
filename=fullfile('C:\Users\Ian\checkout','output',file);
figurePath='C:\Users\Ian\checkout\QoI\docs\IanProposal\figures2';
figureFolder='figures2';
%filename=fullfile('E:\Storage\git_checkout','output',file);
% outputMatrix=['pid' 'QoISNf' 'QoISNa' 'QoIVEFf' 'QoIVEFa'];
% outputMatrix=[outputMatrix '%sigtPert' '%sigaPert' '%sourcePert' '%incPert'];
% outputMatrix=[outputMatrix 'dQoISNf' 'dQoIVEFf' 'dQoISNa' 'dQoIVEFa' 'Ediff'];
% Iterators for perturbation factors
% linspace(-0.1,0.1,21);
sigaPertFactor=linspace(-0.1,0.1,21);
sigsPertFactor=linspace(-0.1,0.1,21);
sourcePertFactor=linspace(-0.1,0.1,21);
incPertFactor=linspace(-0.1,0.1,21);
%Begin iterations over perturbations.
dat.sigaPert = dat.sigaPertRegion.*0;
dat.sigsPert = dat.sigsPertRegion.*0;
dat.sigtPert = dat.sigaPert + dat.sigsPert;
dat.sourcePert =dat.sourcePertRegion.*0;
dat.psiIncPert = dat.incPertRegion.*0;
sigaPertValues=[];
sigaSensValues=[];
sigsPertValues=[];
sigsSensValues=[];
qPertValues=[];
qSensValues=[];
incPertValues=[];
incSensValues=[];
datStart=dat;
for ii=1:numel(sigaPertFactor)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load Perturbations. Used in adjoint sensitivity
    dat.sigaPert = dat.sigaPertRegion.*dat.siga*sigaPertFactor(ii);
    dat.sigtPert = dat.sigaPert;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perturbations
    dat_saved = dat;
    dat.sigt = dat.sigt + dat.sigtPert;
    dat.siga = dat.siga + dat.sigaPert;
    dat.cdif = 1./(3*dat.sigt);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % solve perturbed forward transport problem using sweeps
    [results.phi_pert,results.E_pert,results.Ebd_pert,results.psi_pert]=solveForwardSn;
    [results.phiVEF_pert]=solveForwardVEF;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    calcForwardSensitivity
    dat = dat_saved;
    calcAdjointSensitivity
    Rel_L1_diff=find_Eddington_diff(results.E,results.E_pert);
    sigaSensValues=[sigaSensValues; [results.delta_qoi_sn_f  results.delta_qoi_VEF_f  results.delta_qoi_sn_a  results.delta_qoi_VEF_a]];
end
dat=datStart;
for ii=1:numel(sigsPertFactor)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load Perturbations. Used in adjoint sensitivity
    dat.sigsPert = dat.sigsPertRegion.*dat.sigs*sigsPertFactor(ii);
    dat.sigtPert = dat.sigsPert;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perturbations
    dat_saved = dat;
    dat.sigt = dat.sigt + dat.sigtPert;
    dat.sigs = dat.sigs + dat.sigsPert;
    dat.cdif = 1./(3*dat.sigt);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % solve perturbed forward transport problem using sweeps
    [results.phi_pert,results.E_pert,results.Ebd_pert,results.psi_pert]=solveForwardSn;
    [results.phiVEF_pert]=solveForwardVEF;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    calcForwardSensitivity
    dat = dat_saved;
    calcAdjointSensitivity
    Rel_L1_diff=find_Eddington_diff(results.E,results.E_pert);
    sigsSensValues=[sigsSensValues; [results.delta_qoi_sn_f  results.delta_qoi_VEF_f  results.delta_qoi_sn_a  results.delta_qoi_VEF_a]];
end
dat=datStart;
for ii=1:numel(sourcePertFactor)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load Perturbations. Used in adjoint sensitivity
    dat.sourcePert =dat.sourcePertRegion.*dat.qv_forward*sourcePertFactor(ii);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perturbations
    dat_saved = dat;
    dat.qv_forward = dat.qv_forward + dat.sourcePert;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % solve perturbed forward transport problem using sweeps
    [results.phi_pert,results.E_pert,results.Ebd_pert,results.psi_pert]=solveForwardSn;
    [results.phiVEF_pert]=solveForwardVEF;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    calcForwardSensitivity
    dat = dat_saved;
    calcAdjointSensitivity
    Rel_L1_diff=find_Eddington_diff(results.E,results.E_pert);
    qSensValues=[qSensValues; [results.delta_qoi_sn_f  results.delta_qoi_VEF_f  results.delta_qoi_sn_a  results.delta_qoi_VEF_a]];
end
dat=datStart;
for ii=1:numel(incPertFactor)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load Perturbations. Used in adjoint sensitivity
    dat.psiIncPert = dat.incPertRegion.*dat.inc_forward*incPertFactor(ii);
    dat.psiIncPert = dat.psiIncPert+dat.incPertRegion*0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perturbations
    dat_saved = dat;
    dat.inc_forward = dat.inc_forward + dat.psiIncPert;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % solve perturbed forward transport problem using sweeps
    [results.phi_pert,results.E_pert,results.Ebd_pert,results.psi_pert]=solveForwardSn;
    [results.phiVEF_pert]=solveForwardVEF;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    calcForwardSensitivity
    dat = dat_saved;
    calcAdjointSensitivity
    Rel_L1_diff=find_Eddington_diff(results.E,results.E_pert);
    incSensValues=[incSensValues; [results.delta_qoi_sn_f  results.delta_qoi_VEF_f  results.delta_qoi_sn_a  results.delta_qoi_VEF_a]];
end

figure(600)
hold on
xlabel('\sigma_a % change')
ylabel('QoI % response')
plot(sigaPertFactor,sigaSensValues(:,1)./(results.qoi_sn_f),'-+r')
plot(sigaPertFactor,sigaSensValues(:,2)./(results.qoi_sn_f),'--+g')
plot(sigaPertFactor,sigaSensValues(:,3)./(results.qoi_sn_f),'-+m')
plot(sigaPertFactor,sigaSensValues(:,4)./(results.qoi_sn_f),'--+b')
legend({'sn forward','VEF forward','sn adjoint','VEF adjoint'},'Position',[0.5 0.80 0.01 0.01])
figureFile=[int2str(dat.pb_ID),'sigaSens'];
dataFile=['data\',int2str(dat.pb_ID),'siga.csv'];
fullPath=fullfile(figurePath,{figureFile,dataFile});
outputMatrix = [sigaPertFactor'  sigaSensValues];
csvwrite(fullPath{2},outputMatrix)
title(dataFile,'interpreter','none')
print(fullPath{1},'-dpng');

figure(601)
hold on
xlabel('\sigma_s % change')
ylabel('QoI % response')
plot(sigsPertFactor,sigsSensValues(:,1)./(results.qoi_sn_f),'-+r')
plot(sigsPertFactor,sigsSensValues(:,2)./(results.qoi_sn_f),'--+g')
plot(sigsPertFactor,sigsSensValues(:,3)./(results.qoi_sn_f),'-+m')
plot(sigsPertFactor,sigsSensValues(:,4)./(results.qoi_sn_f),'--+b')
legend({'sn forward','VEF forward','sn adjoint','VEF adjoint'},'Position',[0.5 0.80 0.01 0.01])
figureFile=[int2str(dat.pb_ID),'sigsSens'];
dataFile=['data\',int2str(dat.pb_ID),'sigs.csv'];
fullPath=fullfile(figurePath,{figureFile,dataFile});
outputMatrix = [sigsPertFactor'  sigsSensValues];
csvwrite(fullPath{2},outputMatrix)
title(dataFile,'interpreter','none')
print(fullPath{1},'-dpng');

figure(602)
hold on
xlabel('q % change')
ylabel('QoI % response')
plot(sourcePertFactor,qSensValues(:,1)./(results.qoi_sn_f),'-+r')
plot(sourcePertFactor,qSensValues(:,2)./(results.qoi_sn_f),'--+g')
plot(sourcePertFactor,qSensValues(:,3)./(results.qoi_sn_f),'-+m')
plot(sourcePertFactor,qSensValues(:,4)./(results.qoi_sn_f),'--+b')
legend({'sn forward','VEF forward','sn adjoint','VEF adjoint'},'Position',[0.5 0.80 0.01 0.01])
figureFile=[int2str(dat.pb_ID),'qSens'];
dataFile=['data\',int2str(dat.pb_ID),'q.csv'];
fullPath=fullfile(figurePath,{figureFile,dataFile});
outputMatrix = [sourcePertFactor'  qSensValues];
csvwrite(fullPath{2},outputMatrix)
title(dataFile,'interpreter','none')
print(fullPath{1},'-dpng');

figure(603)
hold on
xlabel('\Psi^- % change')
ylabel('QoI % response')
plot(incPertFactor,incSensValues(:,1)./(results.qoi_sn_f),'-+r')
plot(incPertFactor,incSensValues(:,2)./(results.qoi_sn_f),'--+g')
plot(incPertFactor,incSensValues(:,3)./(results.qoi_sn_f),'-+m')
plot(incPertFactor,incSensValues(:,4)./(results.qoi_sn_f),'--+b')
legend({'sn forward','VEF forward','sn adjoint','VEF adjoint'},'Position',[0.5 0.80 0.01 0.01])
figureFile=[int2str(dat.pb_ID),'incSens'];
dataFile=['data\',int2str(dat.pb_ID),'inc.csv'];
fullPath=fullfile(figurePath,{figureFile,dataFile});
outputMatrix = [incPertFactor'  incSensValues];
csvwrite(fullPath{2},outputMatrix)
title(dataFile,'interpreter','none')
print(fullPath{1},'-dpng');

do_plot(results.phi,'Sn',610,dat.forward_flux)
do_plot(results.phiVEF,'VEF',610,dat.forward_flux)
figure(611)
figureFile=[int2str(dat.pb_ID),'phi'];
figureName=fullfile(figurePath,figureFile);
print(figureName,'-dpng');

do_plot(results.phia,'Sn',611,dat.adjoint_flux)
do_plot(results.phiVEFa,'VEF',611,dat.adjoint_flux)
figure(612)
figureFile=[int2str(dat.pb_ID),'phia'];
figureName=fullfile(figurePath,figureFile);
print(figureName,'-dpng');

end

function [phi,E,Ebd,psi]=solveForwardSn
global npar dat snq IO_opts 
[phi,E,Ebd,psi]=solve_transport(dat.forward_flux,dat.do_dsa,IO_opts.console_io);
return
end

function [phia,Ea,Ebda,psia]=solveAdjointSn
global npar dat snq IO_opts 
[phia,Ea,Ebda,psia]=solve_transport(dat.adjoint_flux,dat.do_dsa,IO_opts.console_io);
return
end

function phiVEF=solveForwardVEF
global npar dat snq IO_opts results
[phiVEF]=solve_VEF(dat.forward_flux,results.E,results.Ebd);
return
end

function phiVEFa=solveAdjointVEF
global npar dat snq IO_opts results
[phiVEFa]=solve_VEF_math_adjoint(dat.adjoint_flux,results.E,results.Ebd);
return
end

function calcUnpertQOI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global npar dat snq IO_opts results
sn=snq.n_dir;
results.qoi_sn_f = compute_qoi(dat.forward_flux,results.phi,sn,results.psi,results.psia);
results.qoi_sn_a = compute_qoi(dat.adjoint_flux,results.phia,sn,results.psi,results.psia);
results.qoi_VEF_f = compute_qoi(dat.forward_flux,results.phiVEF,~sn,[],[]);
results.qoi_VEF_a = compute_qoi(dat.adjoint_flux,results.phiVEFa,~sn,[],[]);
end

function calcForwardSensitivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global npar dat snq IO_opts results
sn=snq.n_dir;
%Forward (direct) sensitivities
results.qoi_sn_f_pert = compute_qoi(dat.forward_flux,results.phi_pert,sn,[],[]);
results.qoi_VEF_f_pert = compute_qoi(dat.forward_flux,results.phiVEF_pert,~sn,[],[]);
results.delta_qoi_sn_f=results.qoi_sn_f_pert-results.qoi_sn_f;
results.delta_qoi_VEF_f=results.qoi_VEF_f_pert-results.qoi_VEF_f;
end

function calcAdjointSensitivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global npar dat snq IO_opts results
sn=snq.n_dir;
%Basic Adjoint sensitivities
results.delta_qoi_sn_a = compute_perturbed_qoi_Sn(dat.adjoint_flux,results.phia,results.phi,results.psi,results.psia,sn);
results.delta_qoi_VEF_a = compute_perturbed_qoi_VEF(dat.adjoint_flux,results.phiVEFa,results.phiVEF,results.E);
results.delta_qoi_VEF_SNphi = compute_perturbed_qoi_VEF(dat.adjoint_flux,results.phiVEFa,results.phi,results.E);
%Some other unrealistic methods (using values we wouldnt expect to have)
results.delta_qoi_VEF_a_pert = compute_perturbed_qoi_VEF(dat.adjoint_flux,results.phiVEFa,results.phiVEF_pert,results.E_pert);
results.delta_qoi_sn_a_pert = compute_perturbed_qoi_Sn(dat.adjoint_flux,results.phia,results.phi_pert,results.psi_pert,results.psia,sn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Alternate Method
% Compute perturbed QoIs using VEF math adjoint method and unperturbed
% forward, alternate method
%delta_qoi_VEF_a_alt = compute_perturbed_qoi_VEF(dat.adjoint_flux,phiVEFa_alt,phiVEF_alt,Ea);
end

function displayUnpertQOI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n-----BEGIN UNPERTURBED QOI DATA OUTPUT----- \n')
global npar dat snq IO_opts results
fprintf('qoi using sn forward: \t %g \n',results.qoi_sn_f);
fprintf('qoi using sn adjoint: \t %g \n',results.qoi_sn_a);
fprintf('qoi using VEF forward: \t %g \n',results.qoi_VEF_f);
fprintf('qoi using VEF adjoint: \t %g \n',results.qoi_VEF_a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function displaySensitivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global npar dat snq IO_opts results
fprintf('\n-----BEGIN PERTURBATION SENSITIVITY DATA OUTPUT----- \n')
fprintf('delta qoi using 2 sn forward runs: \t \t %g \n',results.delta_qoi_sn_f);
fprintf('delta qoi using 2 VEF forward runs: \t %g \n',results.delta_qoi_VEF_f);
fprintf('delta qoi using sn adjoint: \t\t\t %g \n',results.delta_qoi_sn_a);
fprintf('delta qoi using sn adjoint (pert phi): \t\t\t %g \n',results.delta_qoi_sn_a_pert);
fprintf('delta qoi using VEF math adjoint: \t\t %g \n',results.delta_qoi_VEF_a);
fprintf('delta qoi using VEF math adjoint (phi pert): \t\t %g \n',results.delta_qoi_VEF_a_pert);
fprintf('delta qoi using VEF adj w Sn fwd: \t\t %g \n',results.delta_qoi_VEF_SNphi);
%fprintf('delta qoi using VEF math adjoint alt: \t\t %g \n',results.delta_qoi_VEF_a_alt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function generatePlots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global npar dat snq IO_opts results
do_plot(results.phi,'Sn',0,dat.forward_flux)
do_plot(results.phiVEF,'VEF',0,dat.forward_flux)
do_plot(results.phi_pert,'Sn-pert',0,dat.forward_flux)
do_plot(results.phiVEF_pert,'VEF-pert',0,dat.forward_flux)

do_plot(results.phia,'Sn',100,dat.adjoint_flux)
do_plot(results.phiVEFa,'VEF',100,dat.adjoint_flux)

do_plot(results.E,'E',50,dat.forward_flux,1)
do_plot(results.E_pert,'Epert',50,dat.forward_flux,1)
do_plot(results.E_pert-results.E,'Epert-E',51,dat.forward_flux,2)
end


function displaySensitivityError
global npar dat snq IO_opts results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n-----BEGIN DATA ON SENSITIVITY ERROR----- \n')
fprintf('This section contains data on the error in the sensitivity calculations. \n')
fprintf('In general, the "two forward solves" method will be used as the "true" value. \n')
fprintf('This means that there is a "true" sensitivity for two Sn solves and another \n')
fprintf('From two forward VEF solves. \n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('VEF adjoint error from SN forward: \t\t %g \n',(results.qoi_sn_f_pert - results.qoi_sn_f)-results.delta_qoi_VEF_a);
fprintf('VEF adjoint error from VEF forward: \t %g \n',(results.qoi_VEF_f_pert - results.qoi_VEF_f)-results.delta_qoi_VEF_a);
fprintf('VEF adjoint error from SN adjoint: \t\t %g \n',results.delta_qoi_sn_a-results.delta_qoi_VEF_a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function otherStuff
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
[phiVEFa_Epert]=solve_VEF_math_adjoint(dat.adjoint_flux,E_pert,Ebd_pert);
do_plot(phiVEFa_Epert,'VEF Using E_{pert}',100,dat.adjoint_flux)
do_plot(phiVEFa_Epert-phiVEFa,'VEF_{diff}',120,dat.adjoint_flux)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rel_L1_diff=find_Eddington_diff(E,E_pert);
delta_qoi_VEF_math_a_Epert = compute_perturbed_qoi_VEF(dat.adjoint_flux,phiVEFa_Epert,phiVEF_pert,E_pert);
fprintf('delta qoi using VEF math adjoint and Epert: \t\t %g \n',delta_qoi_VEF_math_a_Epert);
fprintf('Diff from using Epert: \t\t %g \n',delta_qoi_VEF_a-delta_qoi_VEF_math_a_Epert);
fprintf('relative L1 difference in E: \t\t %g \n',Rel_L1_diff);
%adjDiff = compute_qoi(dat.adjoint_flux,phiVEFa_Epert,~sn,[],[])-compute_qoi(dat.adjoint_flux,phiVEFa,~sn,[],[]);
%fprintf('NEED TO USE PERT FIX SOURCE USELESS. Diff from using two VEF adjoints: \t\t %g \n',adjDiff);
end
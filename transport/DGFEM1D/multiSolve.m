function multiSolve
global npar dat snq IO_opts results
fprintf('\n-----BEGIN MULTIPLE SENSITIVITY SOLVE----- \n')

%dEdq=0;
%dEdsa=0;
%dEdss=0;
[dEdq, dEdsa, dEdss,dEdinc, dBdq, dBdsa, dBdss,dBdinc]=getEslopes;

%Specify Output File Path
outputMatrix=[];
today=datestr(now,'mmddyyHHMM');
%[pathstr,name,ext] = fileparts(pwd); 
file=['DGFEM1D_prob',int2str(dat.pb_ID),'_',today,'.csv'];
%filename=fullfile('C:\Users\Ian\checkout','output',file);
%figurePath='C:\Users\Ian\checkout\QoI\docs\IanProposal\figures2';
%figurePath='C:\Users\Ian\Projects\QoI\docs\IanProposal\figures2';
figurePath='C:\Users\Ian\Projects\QoI\docs\Thesis\figures2';
figureFolder='figures2';
%filename=fullfile('E:\Storage\git_checkout','output',file);
% outputMatrix=['pid' 'QoISNf' 'QoISNa' 'QoIVEFf' 'QoIVEFa'];
% outputMatrix=[outputMatrix '%sigtPert' '%sigaPert' '%sourcePert' '%incPert'];
% outputMatrix=[outputMatrix 'dQoISNf' 'dQoIVEFf' 'dQoISNa' 'dQoIVEFa' 'Ediff'];
% Iterators for perturbation factors
% linspace(-0.1,0.1,21);
sigaPertFactor=linspace(-0.1,0.1,11);
sigsPertFactor=linspace(-0.1,0.1,11);
sourcePertFactor=linspace(-0.1,0.1,11);
incPertFactor=linspace(-0.1,0.1,11);
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
    E_interp=results.E+(dEdsa.*sigaPertFactor(ii));
    B_interp=results.Ebd+(dBdsa.*sigaPertFactor(ii));
    [results.deltaE_interp,deltaB_L,deltaB_R,volValues]=compute_deltaE_QoI_term(results.phiVEFa,results.phiVEF,results.phiVEF_pert,results.E,E_interp,results.Ebd,B_interp);
    delta_qoi_VEF_a_Eint=results.delta_qoi_VEF_a+results.deltaE_interp+deltaB_L+deltaB_R;
    sigaSensValues=[sigaSensValues; [results.delta_qoi_sn_f  results.delta_qoi_VEF_f  results.delta_qoi_sn_a  results.delta_qoi_VEF_a results.delta_qoi_blend_a delta_qoi_VEF_a_Eint results.delta_qoi_AltVEF_a]];   
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
    E_interp=results.E+(dEdss.*sigsPertFactor(ii)); 
    B_interp=results.Ebd+(dBdss.*sigsPertFactor(ii));
    [results.deltaE_interp,deltaB_L,deltaB_R,volValues]=compute_deltaE_QoI_term(results.phiVEFa,results.phiVEF,results.phiVEF_pert,results.E,E_interp,results.Ebd,B_interp);
    delta_qoi_VEF_a_Eint=results.delta_qoi_VEF_a+results.deltaE_interp+deltaB_L+deltaB_R;
    sigsSensValues=[sigsSensValues; [results.delta_qoi_sn_f  results.delta_qoi_VEF_f  results.delta_qoi_sn_a  results.delta_qoi_VEF_a results.delta_qoi_blend_a delta_qoi_VEF_a_Eint results.delta_qoi_AltVEF_a]];
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
    E_interp=results.E+(dEdq.*sourcePertFactor(ii));  
    B_interp=results.Ebd+(dBdq.*sourcePertFactor(ii));
    [results.deltaE_interp,deltaB_L,deltaB_R,volValues]=compute_deltaE_QoI_term(results.phiVEFa,results.phiVEF,results.phiVEF_pert,results.E,E_interp,results.Ebd,B_interp);
    delta_qoi_VEF_a_Eint=results.delta_qoi_VEF_a+results.deltaE_interp+deltaB_L+deltaB_R;
    qSensValues=[qSensValues; [results.delta_qoi_sn_f  results.delta_qoi_VEF_f  results.delta_qoi_sn_a  results.delta_qoi_VEF_a results.delta_qoi_blend_a delta_qoi_VEF_a_Eint results.delta_qoi_AltVEF_a]];
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
    E_interp=results.E+(dEdinc.*incPertFactor(ii));
    B_interp=results.Ebd+(dBdinc.*incPertFactor(ii));
    [results.deltaE_interp,deltaB_L,deltaB_R,volValues]=compute_deltaE_QoI_term(results.phiVEFa,results.phiVEF,results.phiVEF_pert,results.E,E_interp,results.Ebd,B_interp);
    delta_qoi_VEF_a_Eint=results.delta_qoi_VEF_a+results.deltaE_interp+deltaB_L+deltaB_R;
    incSensValues=[incSensValues; [results.delta_qoi_sn_f  results.delta_qoi_VEF_f  results.delta_qoi_sn_a  results.delta_qoi_VEF_a results.delta_qoi_blend_a delta_qoi_VEF_a_Eint results.delta_qoi_AltVEF_a]];
end

% %%%%%%%%%Try 3D graph
% dat=datStart;
% multiSensValues=[];
% xPert=[];
% yPert=[];
% for ii=1:numel(sourcePertFactor)
%     for jj=1:numel(sigaPertFactor)
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Load Perturbations. Used in adjoint sensitivity
%         dat.sourcePert =dat.sourcePertRegion.*dat.qv_forward*sourcePertFactor(ii);
%         dat.sigaPert = dat.sigaPertRegion.*dat.siga*sigaPertFactor(jj);
%         dat.sigtPert = dat.sigaPert;
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % perturbations
%         dat_saved = dat;
%         dat.qv_forward = dat.qv_forward + dat.sourcePert;
%         dat.sigt = dat.sigt + dat.sigtPert;
%         dat.siga = dat.siga + dat.sigaPert;
%         dat.cdif = 1./(3*dat.sigt);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % solve perturbed forward transport problem using sweeps
%         [results.phi_pert,results.E_pert,results.Ebd_pert,results.psi_pert]=solveForwardSn;
%         [results.phiVEF_pert]=solveForwardVEF;
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         calcForwardSensitivity
%         dat = dat_saved;
%         calcAdjointSensitivity
%         Rel_L1_diff=find_Eddington_diff(results.E,results.E_pert);
%         E_interp=results.E+(dEdq.*sourcePertFactor(ii))+(dEdsa.*sigaPertFactor(jj));
%         B_interp=results.Ebd+(dBdq.*sourcePertFactor(ii))+(dBdsa.*sigaPertFactor(jj));
%         [results.deltaE_interp,deltaB_L,deltaB_R,volValues]=compute_deltaE_QoI_term(results.phiVEFa,results.phiVEF,results.phiVEF_pert,results.E,E_interp,results.Ebd,B_interp);
%         delta_qoi_VEF_a_Eint=results.delta_qoi_VEF_a+results.deltaE_interp;
%         multiSensValues=[multiSensValues; [results.delta_qoi_sn_f  results.delta_qoi_VEF_f  results.delta_qoi_sn_a  results.delta_qoi_VEF_a results.delta_qoi_blend_a delta_qoi_VEF_a_Eint results.delta_qoi_AltVEF_a]];
%         xPert=[xPert sourcePertFactor(ii)];
%         yPert=[yPert sigaPertFactor(jj)];
%     end
% end

% figure(500)
% hold on
% xlabel('q % change')
% ylabel('\sigma_a % change')
% zlabel('QoI % response')
% trisurf(delaunay(xPert,yPert),xPert,yPert,multiSensValues(:,1)./(results.qoi_sn_f),[],'FaceColor',[0 1 0])
% trisurf(delaunay(xPert,yPert),xPert,yPert,multiSensValues(:,2)./(results.qoi_sn_f),[],'FaceColor',[1 0 0])
% trisurf(delaunay(xPert,yPert),xPert,yPert,multiSensValues(:,3)./(results.qoi_sn_f),[],'FaceColor',[0 0 1])
% trisurf(delaunay(xPert,yPert),xPert,yPert,multiSensValues(:,4)./(results.qoi_sn_f),[],'FaceColor',[1 1 1])
% trisurf(delaunay(xPert,yPert),xPert,yPert,multiSensValues(:,5)./(results.qoi_sn_f),[],'FaceColor',[0 1 1])
% trisurf(delaunay(xPert,yPert),xPert,yPert,multiSensValues(:,6)./(results.qoi_sn_f),[],'FaceColor',[1 1 0])
% % plot3(xPert,yPert,multiSensValues(:,1)./(results.qoi_sn_f),'+g')
% % plot3(xPert,yPert,multiSensValues(:,2)./(results.qoi_sn_f),'or')
% % plot3(xPert,yPert,multiSensValues(:,3)./(results.qoi_sn_f),'db')
% % plot3(xPert,yPert,multiSensValues(:,4)./(results.qoi_sn_f),'+k')
% % plot3(xPert,yPert,multiSensValues(:,5)./(results.qoi_sn_f),'oc')
% % plot3(xPert,yPert,multiSensValues(:,6)./(results.qoi_sn_f),'dm')
% legend({'sn forward','VEF forward','sn adjoint','VEF adjoint','Blend adjoint','VEF adjoint E_{appx}'},'Position',[0.5 0.80 0.01 0.01])
% figureFile=[int2str(dat.pb_ID),'plot3d'];
% dataFile=['data\',int2str(dat.pb_ID),'plot3d.csv'];
% fullPath=fullfile(figurePath,{figureFile,dataFile});
% outputMatrix = [sigaPertFactor'  sigaSensValues];
% csvwrite(fullPath{2},outputMatrix)
% title(dataFile,'interpreter','none')
% print(fullPath{1},'-dpng');


%%%%%%End 3D graph
%%%%%%Start multi Pert
dat=datStart;
twoSensValues=[];
for ii=1:numel(sourcePertFactor)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load Perturbations. Used in adjoint sensitivity
    dat.sourcePert =dat.sourcePertRegion.*dat.qv_forward*1*sourcePertFactor(ii);
    dat.sigaPert = dat.sigaPertRegion.*dat.siga*-1*sigaPertFactor(ii);
    dat.sigtPert = dat.sigaPert;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perturbations
    dat_saved = dat;
    dat.sigt = dat.sigt + dat.sigtPert;
    dat.siga = dat.siga + dat.sigaPert;
    dat.cdif = 1./(3*dat.sigt);
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
    E_interp=results.E+(dEdq.*sourcePertFactor(ii))+(-1*dEdsa.*sourcePertFactor(ii));  
    B_interp=results.Ebd+(dBdq.*sourcePertFactor(ii))+(-1*dBdsa.*sourcePertFactor(ii));
    [results.deltaE_interp,deltaB_L,deltaB_R,volValues]=compute_deltaE_QoI_term(results.phiVEFa,results.phiVEF,results.phiVEF_pert,results.E,E_interp,results.Ebd,B_interp);
    delta_qoi_VEF_a_Eint=results.delta_qoi_VEF_a+results.deltaE_interp+deltaB_L+deltaB_R;
    twoSensValues=[twoSensValues; [results.delta_qoi_sn_f  results.delta_qoi_VEF_f  results.delta_qoi_sn_a  results.delta_qoi_VEF_a results.delta_qoi_blend_a delta_qoi_VEF_a_Eint results.delta_qoi_AltVEF_a]];
end
figure(502)
hold on
xlabel('q % change, - \sigma_a % change')
ylabel('QoI % response')
plot(sourcePertFactor,twoSensValues(:,1)./(results.qoi_sn_f),'-+g')
%plot(sourcePertFactor,twoSensValues(:,2)./(results.qoi_sn_f),'-or')
plot(sourcePertFactor,twoSensValues(:,3)./(results.qoi_sn_f),'-db')
plot(sourcePertFactor,twoSensValues(:,4)./(results.qoi_sn_f),'--+k')
plot(sourcePertFactor,twoSensValues(:,5)./(results.qoi_sn_f),'--oc')
plot(sourcePertFactor,twoSensValues(:,6)./(results.qoi_sn_f),'--dm')
plot(sourcePertFactor,twoSensValues(:,7)./(results.qoi_sn_f),'--+r')
legend({'sn forward','sn adjoint','VET adjoint','Blend adjoint','VET adjoint E_{appx}','aVET'},'Position',[0.5 0.78 0.01 0.01])
figureFile=[int2str(dat.pb_ID),'qsigaSens'];
dataFile=['data\',int2str(dat.pb_ID),'qsiga.csv'];
fullPath=fullfile(figurePath,{figureFile,dataFile});
outputMatrix = [sourcePertFactor'  twoSensValues];
csvwrite(fullPath{2},outputMatrix)
title(dataFile,'interpreter','none')
print(fullPath{1},'-dpng');



%%%%%%Start multi Pert
dat=datStart;
twoSensValues=[];
for ii=1:numel(incPertFactor)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load Perturbations. Used in adjoint sensitivity
    dat.psiIncPert = dat.incPertRegion.*dat.inc_forward*incPertFactor(ii);
    dat.sigaPert = dat.sigaPertRegion.*dat.siga*-1*sigaPertFactor(ii);
    dat.sigtPert = dat.sigaPert;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perturbations
    dat_saved = dat;
    dat.sigt = dat.sigt + dat.sigtPert;
    dat.siga = dat.siga + dat.sigaPert;
    dat.cdif = 1./(3*dat.sigt);
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
    E_interp=results.E+(dEdinc.*incPertFactor(ii))+(-1*dEdsa.*incPertFactor(ii));  
    B_interp=results.Ebd+(dBdinc.*incPertFactor(ii))+(-1*dBdsa.*incPertFactor(ii));
    [results.deltaE_interp,deltaB_L,deltaB_R,volValues]=compute_deltaE_QoI_term(results.phiVEFa,results.phiVEF,results.phiVEF_pert,results.E,E_interp,results.Ebd,B_interp);
    delta_qoi_VEF_a_Eint=results.delta_qoi_VEF_a+results.deltaE_interp+deltaB_L+deltaB_R;
    twoSensValues=[twoSensValues; [results.delta_qoi_sn_f  results.delta_qoi_VEF_f  results.delta_qoi_sn_a  results.delta_qoi_VEF_a results.delta_qoi_blend_a delta_qoi_VEF_a_Eint results.delta_qoi_AltVEF_a]];
end
figure(503)
hold on
xlabel('\Psi^- % change, - \sigma_a % change')
ylabel('QoI % response')
plot(incPertFactor,twoSensValues(:,1)./(results.qoi_sn_f),'-+g')
%plot(incPertFactor,twoSensValues(:,2)./(results.qoi_sn_f),'-or')
plot(incPertFactor,twoSensValues(:,3)./(results.qoi_sn_f),'-db')
plot(incPertFactor,twoSensValues(:,4)./(results.qoi_sn_f),'--+k')
plot(incPertFactor,twoSensValues(:,5)./(results.qoi_sn_f),'--oc')
plot(incPertFactor,twoSensValues(:,6)./(results.qoi_sn_f),'--dm')
plot(incPertFactor,twoSensValues(:,7)./(results.qoi_sn_f),'--+r')
legend({'sn forward','sn adjoint','VET adjoint','Blend adjoint','VET adjoint E_{appx}','aVET'},'Position',[0.5 0.78 0.01 0.01])
figureFile=[int2str(dat.pb_ID),'incsigaSens'];
dataFile=['data\',int2str(dat.pb_ID),'incsiga.csv'];
fullPath=fullfile(figurePath,{figureFile,dataFile});
outputMatrix = [sourcePertFactor'  twoSensValues];
csvwrite(fullPath{2},outputMatrix)
title(dataFile,'interpreter','none')
print(fullPath{1},'-dpng');

%%%%%%End multi Pert

figure(600)
hold on
xlabel('\sigma_a % change')
ylabel('QoI % response')
plot(sigaPertFactor,sigaSensValues(:,1)./(results.qoi_sn_f),'-+g')
%plot(sigaPertFactor,sigaSensValues(:,2)./(results.qoi_sn_f),'-or')
plot(sigaPertFactor,sigaSensValues(:,3)./(results.qoi_sn_f),'-db')
plot(sigaPertFactor,sigaSensValues(:,4)./(results.qoi_sn_f),'--+k')
plot(sigaPertFactor,sigaSensValues(:,5)./(results.qoi_sn_f),'--oc')
plot(sigaPertFactor,sigaSensValues(:,6)./(results.qoi_sn_f),'--dm')
plot(sigaPertFactor,sigaSensValues(:,7)./(results.qoi_sn_f),'--+r')
legend({'sn forward','sn adjoint','VET adjoint','Blend adjoint','VET adjoint E_{appx}','aVET'},'Position',[0.5 0.78 0.01 0.01])
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
plot(sigsPertFactor,sigsSensValues(:,1)./(results.qoi_sn_f),'-+g')
%plot(sigsPertFactor,sigsSensValues(:,2)./(results.qoi_sn_f),'-or')
plot(sigsPertFactor,sigsSensValues(:,3)./(results.qoi_sn_f),'-db')
plot(sigsPertFactor,sigsSensValues(:,4)./(results.qoi_sn_f),'--+k')
plot(sigsPertFactor,sigsSensValues(:,5)./(results.qoi_sn_f),'--oc')
plot(sigsPertFactor,sigsSensValues(:,6)./(results.qoi_sn_f),'--dm')
plot(sigsPertFactor,sigsSensValues(:,7)./(results.qoi_sn_f),'--+r')
legend({'sn forward','sn adjoint','VET adjoint','Blend adjoint','VET adjoint E_{appx}','aVET'},'Position',[0.5 0.78 0.01 0.01])
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
plot(sourcePertFactor,qSensValues(:,1)./(results.qoi_sn_f),'-+g')
%plot(sourcePertFactor,qSensValues(:,2)./(results.qoi_sn_f),'-or')
plot(sourcePertFactor,qSensValues(:,3)./(results.qoi_sn_f),'-db')
plot(sourcePertFactor,qSensValues(:,4)./(results.qoi_sn_f),'--+k')
plot(sourcePertFactor,qSensValues(:,5)./(results.qoi_sn_f),'--oc')
plot(sourcePertFactor,qSensValues(:,6)./(results.qoi_sn_f),'--dm')
plot(sourcePertFactor,qSensValues(:,7)./(results.qoi_sn_f),'--+r')
legend({'sn forward','sn adjoint','VET adjoint','Blend adjoint','VET adjoint E_{appx}','aVET'},'Position',[0.5 0.78 0.01 0.01])
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
plot(incPertFactor,incSensValues(:,1)./(results.qoi_sn_f),'-+g')
%plot(incPertFactor,incSensValues(:,2)./(results.qoi_sn_f),'-or')
plot(incPertFactor,incSensValues(:,3)./(results.qoi_sn_f),'-db')
plot(incPertFactor,incSensValues(:,4)./(results.qoi_sn_f),'--+k')
plot(incPertFactor,incSensValues(:,5)./(results.qoi_sn_f),'--oc')
plot(incPertFactor,incSensValues(:,6)./(results.qoi_sn_f),'--dm')
plot(incPertFactor,incSensValues(:,7)./(results.qoi_sn_f),'--+r')
legend({'sn forward','sn adjoint','VET adjoint','Blend adjoint','VET adjoint E_{appx}','aVET'},'Position',[0.5 0.78 0.01 0.01])
figureFile=[int2str(dat.pb_ID),'incSens'];
dataFile=['data\',int2str(dat.pb_ID),'inc.csv'];
fullPath=fullfile(figurePath,{figureFile,dataFile});
outputMatrix = [incPertFactor'  incSensValues];
csvwrite(fullPath{2},outputMatrix)
title(dataFile,'interpreter','none')
print(fullPath{1},'-dpng');

do_plot(results.phi,'$$\phi$$ Sn',610,dat.forward_flux)
do_plot(results.phiVEF,'$$\phi$$ VET',610,dat.forward_flux)
do_plot(results.phiAltVEF,'$$\varphi$$ aVET',610,dat.forward_flux)
figure(611)
figureFile=[int2str(dat.pb_ID),'phi'];
figureName=fullfile(figurePath,figureFile);
print(figureName,'-dpng');

do_plot(results.phia,'$$\phi^\dag$$ Sn',611,dat.adjoint_flux)
do_plot(results.phiVEFa,'$$\varphi^\dag$$ VET',611,dat.adjoint_flux)
do_plot(results.phiAltVEFa,'$$\phi^\dag$$ aVET',611,dat.adjoint_flux)
figure(612)
figureFile=[int2str(dat.pb_ID),'phia'];
figureName=fullfile(figurePath,figureFile);
print(figureName,'-dpng');

do_plot(results.phi,'$$\phi$$ Sn',612,dat.forward_flux)
do_plot(results.phiVEF,'$$\phi$$ VET',612,dat.forward_flux)
do_plot(results.phiAltVEF,'$$\varphi$$ aVET',612,dat.forward_flux)
do_plot(2*results.phiAltVEF,'$$2\varphi$$ aVET',612,dat.forward_flux)
figure(613)
figureFile=[int2str(dat.pb_ID),'phi2'];
figureName=fullfile(figurePath,figureFile);
print(figureName,'-dpng');

do_plot(results.phia,'$$\phi^\dag$$ Sn',613,dat.adjoint_flux)
do_plot(results.phiVEFa,'$$\varphi^\dag$$ VET',613,dat.adjoint_flux)
do_plot(results.phiAltVEFa,'$$\phi^\dag$$ aVET',613,dat.adjoint_flux)
do_plot(2*results.phiVEFa,'$$2\varphi^\dag$$ VET',613,dat.adjoint_flux)
figure(614)
figureFile=[int2str(dat.pb_ID),'phia2'];
figureName=fullfile(figurePath,figureFile);
print(figureName,'-dpng');
end
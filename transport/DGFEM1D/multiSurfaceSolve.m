function multiSurfaceSolve
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
%%%%%%%%%Try 3D graph
dat=datStart;
multiSensValues=[];
xPert=[];
yPert=[];
for ii=1:numel(sourcePertFactor)
    for jj=1:numel(sigaPertFactor)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Load Perturbations. Used in adjoint sensitivity
        dat.sourcePert =dat.sourcePertRegion.*dat.qv_forward*sourcePertFactor(ii);
        dat.sigaPert = dat.sigaPertRegion.*dat.siga*sigaPertFactor(jj);
        dat.sigtPert = dat.sigaPert;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % perturbations
        dat_saved = dat;
        dat.qv_forward = dat.qv_forward + dat.sourcePert;
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
        E_interp=results.E+(dEdq.*sourcePertFactor(ii))+(dEdsa.*sigaPertFactor(jj));
        B_interp=results.Ebd+(dBdq.*sourcePertFactor(ii))+(dBdsa.*sigaPertFactor(jj));
        [results.deltaE_interp,deltaB_L,deltaB_R,volValues]=compute_deltaE_QoI_term(results.phiVEFa,results.phiVEF,results.phiVEF_pert,results.E,E_interp,results.Ebd,B_interp);
        delta_qoi_VEF_a_Eint=results.delta_qoi_VEF_a+results.deltaE_interp;
        multiSensValues=[multiSensValues; [results.delta_qoi_sn_f  results.delta_qoi_VEF_f  results.delta_qoi_sn_a  results.delta_qoi_VEF_a results.delta_qoi_blend_a delta_qoi_VEF_a_Eint]];
        xPert=[xPert sourcePertFactor(ii)];
        yPert=[yPert sigaPertFactor(jj)];
    end
end

figure(500)
hold on
xlabel('q % change')
ylabel('\sigma_a % change')
zlabel('QoI % response')
trisurf(delaunay(xPert,yPert),xPert,yPert,multiSensValues(:,1)./(results.qoi_sn_f),[],'FaceColor',[0 1 0])
trisurf(delaunay(xPert,yPert),xPert,yPert,multiSensValues(:,2)./(results.qoi_sn_f),[],'FaceColor',[1 0 0])
trisurf(delaunay(xPert,yPert),xPert,yPert,multiSensValues(:,3)./(results.qoi_sn_f),[],'FaceColor',[0 0 1])
trisurf(delaunay(xPert,yPert),xPert,yPert,multiSensValues(:,4)./(results.qoi_sn_f),[],'FaceColor',[1 1 1])
trisurf(delaunay(xPert,yPert),xPert,yPert,multiSensValues(:,5)./(results.qoi_sn_f),[],'FaceColor',[0 1 1])
trisurf(delaunay(xPert,yPert),xPert,yPert,multiSensValues(:,6)./(results.qoi_sn_f),[],'FaceColor',[1 1 0])
% plot3(xPert,yPert,multiSensValues(:,1)./(results.qoi_sn_f),'+g')
% plot3(xPert,yPert,multiSensValues(:,2)./(results.qoi_sn_f),'or')
% plot3(xPert,yPert,multiSensValues(:,3)./(results.qoi_sn_f),'db')
% plot3(xPert,yPert,multiSensValues(:,4)./(results.qoi_sn_f),'+k')
% plot3(xPert,yPert,multiSensValues(:,5)./(results.qoi_sn_f),'oc')
% plot3(xPert,yPert,multiSensValues(:,6)./(results.qoi_sn_f),'dm')
legend({'sn forward','VEF forward','sn adjoint','VEF adjoint','Blend adjoint','VEF adjoint E_{appx}'},'Position',[0.5 0.80 0.01 0.01])
figureFile=[int2str(dat.pb_ID),'plot3d'];
dataFile=['data\',int2str(dat.pb_ID),'plot3d.csv'];
fullPath=fullfile(figurePath,{figureFile,dataFile});
outputMatrix = [sigaPertFactor'  sigaSensValues];
csvwrite(fullPath{2},outputMatrix)
title(dataFile,'interpreter','none')
print(fullPath{1},'-dpng');


%%%%%%3D graph
end
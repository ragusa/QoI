function [dEdq, dEdsa, dEdss,dEdinc,dBdq, dBdsa, dBdss,dBdinc]=getEslopes
global dat results
dat_saved=dat;
factor=0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Try to estimate Ep
dat.sourcePert =dat.sourcePertRegion.*dat.qv_forward*factor;
dat.sigaPert = dat.sigaPertRegion.*dat.siga*0.0+dat.sigaPertRegion.*0;
dat.sigsPert = dat.sigsPertRegion.*dat.sigs*0.0+dat.sigsPertRegion.*0;
dat.sigtPert = dat.sigaPert + dat.sigsPert;
dat.qv_forward = dat.qv_forward + dat.sourcePert;
dat.sigs = dat.sigs + dat.sigsPert;
dat.sigt = dat.sigt + dat.sigtPert;
dat.siga = dat.siga + dat.sigaPert;
% solve perturbed forward transport problem using sweeps
[Dummyphi_pert,E_pert1,Ebd_pert1,Dummypsi_pert]=solveForwardSn;
dEdq=((E_pert1-results.E)./factor);
dBdq=((Ebd_pert1-results.Ebd)./factor);
dat=dat_saved;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Try to estimate Ep
dat.sourcePert =dat.sourcePertRegion.*dat.qv_forward*0.0;
dat.sigaPert = dat.sigaPertRegion.*dat.siga*factor+dat.sigaPertRegion.*0;
dat.sigsPert = dat.sigsPertRegion.*dat.sigs*0.0+dat.sigsPertRegion.*0;
dat.sigtPert = dat.sigaPert + dat.sigsPert;
dat.qv_forward = dat.qv_forward + dat.sourcePert;
dat.sigs = dat.sigs + dat.sigsPert;
dat.sigt = dat.sigt + dat.sigtPert;
dat.siga = dat.siga + dat.sigaPert;
% solve perturbed forward transport problem using sweeps
[Dummyphi_pert,E_pert1,Ebd_pert1,Dummypsi_pert]=solveForwardSn;
dEdsa=((E_pert1-results.E)./factor);
dBdsa=((Ebd_pert1-results.Ebd)./factor);
dat=dat_saved;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Try to estimate Ep
dat.sourcePert =dat.sourcePertRegion.*dat.qv_forward*0.0;
dat.sigaPert = dat.sigaPertRegion.*dat.siga*0.0+dat.sigaPertRegion.*0;
dat.sigsPert = dat.sigsPertRegion.*dat.sigs*factor+dat.sigsPertRegion.*0;
dat.sigtPert = dat.sigaPert + dat.sigsPert;
dat.qv_forward = dat.qv_forward + dat.sourcePert;
dat.sigs = dat.sigs + dat.sigsPert;
dat.sigt = dat.sigt + dat.sigtPert;
dat.siga = dat.siga + dat.sigaPert;
% solve perturbed forward transport problem using sweeps
[Dummyphi_pert,E_pert1,Ebd_pert1,Dummypsi_pert]=solveForwardSn;
dEdss=((E_pert1-results.E)./factor);
dBdss=((Ebd_pert1-results.Ebd)./factor);
dat=dat_saved;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Try to estimate Ep
dat.sourcePert =dat.sourcePertRegion.*dat.qv_forward*0.0;
dat.sigaPert = dat.sigaPertRegion.*dat.siga*0.0+dat.sigaPertRegion.*0;
dat.sigsPert = dat.sigsPertRegion.*dat.sigs*0.0+dat.sigsPertRegion.*0;
dat.sigtPert = dat.sigaPert + dat.sigsPert;
dat.psiIncPert = dat.incPertRegion.*dat.inc_forward*factor;
dat.qv_forward = dat.qv_forward + dat.sourcePert;
dat.sigs = dat.sigs + dat.sigsPert;
dat.sigt = dat.sigt + dat.sigtPert;
dat.siga = dat.siga + dat.sigaPert;
dat.inc_forward = dat.inc_forward + dat.psiIncPert;
% solve perturbed forward transport problem using sweeps
[Dummyphi_pert,E_pert1,Ebd_pert1,Dummypsi_pert]=solveForwardSn;
dEdinc=((E_pert1-results.E)./factor);
dBdinc=((Ebd_pert1-results.Ebd)./factor);
dat=dat_saved;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return 
end
function singleSolve
global npar dat snq IO_opts results

[dEdq, dEdsa, dEdss, dEdinc]=getEslopes;

qPertFac=0.0;
sigaPertFac=0.0;
sigsPertFac=0.1;
incPertFac=0.0;
% Load Perturbations. Used in adjoint sensitivity
dat.sigaPert = dat.sigaPertRegion.*dat.siga*sigaPertFac+dat.sigaPertRegion.*0;
dat.sigsPert = dat.sigsPertRegion.*dat.sigs*sigsPertFac+dat.sigsPertRegion.*0.0;
dat.sigtPert = dat.sigaPert + dat.sigsPert;
dat.sourcePert =dat.sourcePertRegion.*dat.qv_forward*qPertFac;
dat.psiIncPert = dat.incPertRegion.*dat.inc_forward*incPertFac+dat.incPertRegion*0.0;
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
calcForwardSensitivity
dat=dat_saved;
E_interp=results.E+(dEdq.*qPertFac+dEdsa.*sigaPertFac+dEdss.*sigsPertFac+dEdinc*incPertFac);
calcAdjointSensitivity

% delta_qoi_VEF_interp = compute_perturbed_qoi_VEF(dat.adjoint_flux,results.phiVEFa,results.phiVEF,E_interp);
% fprintf('delta qoi using VEF math adjoint E Interp: \t\t %g \n',delta_qoi_VEF_interp);
% delta_qoi_VEF_Epert = compute_perturbed_qoi_VEF(dat.adjoint_flux,results.phiVEFa,results.phiVEF,results.E_pert);
% fprintf('delta qoi using VEF math adjoint E Interp: \t\t %g \n',delta_qoi_VEF_Epert);
[results.deltaE_interp,results.deltaB_L,results.deltaB_R,volValues]=compute_deltaE_QoI_term(results.phiVEFa,results.phiVEF,results.phiVEF_pert,results.E,E_interp,results.Ebd,results.Ebd_pert);

fprintf('Delta source term: \t\t %g \n',results.dqoiparts_VEF_a(1));
fprintf('Delta abs term: \t\t %g \n',results.dqoiparts_VEF_a(2));
fprintf('Delta scat term: \t\t %g \n',results.dqoiparts_VEF_a(3));
fprintf('Delta inc term: \t\t %g \n',results.dqoiparts_VEF_a(4));
fprintf('E Interp Correction: \t\t %g \n',results.deltaE_interp);
fprintf('B Interp Correction: \t\t %g \n',results.deltaB_L+results.deltaB_R);

displaySensitivity
generatePlots
displaySensitivityError
eddingtonTerms
%do_plot(E_interp,'Einterp',150,dat.forward_flux,1)
do_plot(E_interp-results.E,'$$\delta$$ E estimate',252,dat.forward_flux,1)
do_plot(results.E_pert-results.E,'$$\delta$$ E real',252,dat.forward_flux,1)

do_plot(results.phi,'$$\phi$$',253,dat.forward_flux)
do_plot(results.phi_pert,'$$\phi$$ pert',253,dat.forward_flux)

do_plot(results.phi,'$$\phi$$ Sn',300,dat.forward_flux)
do_plot(results.phiVEF,'$$\phi$$ VET',300,dat.forward_flux)

do_plot(results.phia,'$$\phi^\dagger$$ Sn',301,dat.adjoint_flux)
do_plot(results.phiVEFa,'$$\varphi^\dagger$$ VET',301,dat.adjoint_flux)


results.Ebd
results.Ebd_pert

results.dqoiparts_VEF_a
end
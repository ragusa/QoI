function generatePlots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global npar dat snq IO_opts results
do_plot(results.phi,'Sn',0,dat.forward_flux)
do_plot(results.phiVEF,'VET',0,dat.forward_flux)
do_plot(results.phi_pert,'Sn-pert',0,dat.forward_flux)
do_plot(results.phiVEF_pert,'VET-pert',0,dat.forward_flux)
%do_plot(results.phiAltVEF,'AltVET',0,dat.forward_flux)

do_plot(results.phia,'Sn',100,dat.adjoint_flux)
do_plot(results.phiVEFa,'VET',100,dat.adjoint_flux)
%do_plot(2*results.phiVEFa,'2*VET',100,dat.adjoint_flux)
%do_plot(results.phiAltVEFa,'AltVET',100,dat.adjoint_flux)

do_plot(results.E,'E',149,dat.forward_flux,1)
do_plot(results.E,'E',150,dat.forward_flux,1)
do_plot(results.E_pert,'Epert',150,dat.forward_flux,1)
do_plot(results.E_pert-results.E,'$$\delta E$$',151,dat.forward_flux,2)

do_plot(100*(results.E_pert-results.E)./results.E,'$$\delta E$$ ($$\%$$ change)',155,dat.forward_flux,2)
end
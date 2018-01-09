function generatePlots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global npar dat snq IO_opts results
do_plot(results.phi,'Sn',0,dat.forward_flux)
do_plot(results.phiVEF,'VEF',0,dat.forward_flux)
do_plot(results.phi_pert,'Sn-pert',0,dat.forward_flux)
do_plot(results.phiVEF_pert,'VEF-pert',0,dat.forward_flux)

do_plot(results.phia,'Sn',100,dat.adjoint_flux)
do_plot(results.phiVEFa,'VEF',100,dat.adjoint_flux)
do_plot(2*results.phiVEFa,'2*VEF',100,dat.adjoint_flux)

do_plot(results.E,'E',150,dat.forward_flux,1)
do_plot(results.E_pert,'Epert',150,dat.forward_flux,1)
do_plot(results.E_pert-results.E,'Epert-E',151,dat.forward_flux,2)
end
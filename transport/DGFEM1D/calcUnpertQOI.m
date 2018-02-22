function calcUnpertQOI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global npar dat snq IO_opts results
sn=snq.n_dir;
results.qoi_sn_f = compute_qoi(dat.forward_flux,results.phi,sn,results.psi,results.psia);
results.qoi_sn_a = compute_qoi(dat.adjoint_flux,results.phia,sn,results.psi,results.psia);
results.qoi_VEF_f = compute_qoi(dat.forward_flux,results.phiVEF,~sn,[],[]);
results.qoi_VEF_a = compute_qoi(dat.adjoint_flux,results.phiVEFa,~sn,[],[]);
results.qoi_AltVEF_f = compute_qoi_Alt(dat.forward_flux,results.phiAltVEF,~sn,results.psi,results.psia,results.phia);
results.qoi_AltVEF_a = compute_qoi_Alt(dat.adjoint_flux,results.phiAltVEFa,~sn,results.psi,results.psia,results.phia);
end
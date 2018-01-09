function calcUnpertQOI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global npar dat snq IO_opts results
sn=snq.n_dir;
results.qoi_sn_f = compute_qoi(dat.forward_flux,results.phi,sn,results.psi,results.psia);
results.qoi_sn_a = compute_qoi(dat.adjoint_flux,results.phia,sn,results.psi,results.psia);
results.qoi_VEF_f = compute_qoi(dat.forward_flux,results.phiVEF,~sn,[],[]);
results.qoi_VEF_a = compute_qoi(dat.adjoint_flux,results.phiVEFa,~sn,[],[]);
end
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
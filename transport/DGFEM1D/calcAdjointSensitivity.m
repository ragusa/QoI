function calcAdjointSensitivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global npar dat snq IO_opts results
sn=snq.n_dir;
%Basic Adjoint sensitivities
results.delta_qoi_sn_a = compute_perturbed_qoi_Sn(dat.adjoint_flux,results.phia,results.phi,results.psi,results.psia,sn);
results.delta_qoi_VEF_a = compute_perturbed_qoi_VEF(dat.adjoint_flux,results.phiVEFa,results.phiVEF,results.E);
results.delta_qoi_blend_a = compute_perturbed_qoi_VEF(dat.adjoint_flux,results.phiVEFa,results.phiVEF,results.E,results.phia,results.psia);
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
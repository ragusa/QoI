function displaySensitivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global npar dat snq IO_opts results
fprintf('\n-----BEGIN PERTURBATION SENSITIVITY DATA OUTPUT----- \n')
fprintf('delta qoi using 2 sn forward runs: \t \t %g \n',results.delta_qoi_sn_f);
fprintf('delta qoi using 2 VEF forward runs: \t %g \n',results.delta_qoi_VEF_f);
fprintf('delta qoi using sn adjoint: \t\t\t %g \n',results.delta_qoi_sn_a);
fprintf('delta qoi using VEF adjoint: \t\t %g \n',results.delta_qoi_VEF_a);
fprintf('delta qoi using VEF adjoint (E interp): \t\t %g \n',results.delta_qoi_VEF_a+results.deltaE_interp);
fprintf('delta qoi using VEF-blend adjoint: \t\t %g \n',results.delta_qoi_blend_a);
fprintf('delta qoi using VEF-blend adjoint (E interp): \t\t %g \n',results.delta_qoi_blend_a+results.deltaE_interp);
fprintf('\n-----BEGIN PERTURBATION SENSITIVITY DATA OUTPUT (2nd set)----- \n')
fprintf('delta qoi using sn adjoint (pert phi): \t\t\t %g \n',results.delta_qoi_sn_a_pert);
fprintf('delta qoi using VEF adjoint (phi pert): \t\t %g \n',results.delta_qoi_VEF_a_pert);
fprintf('delta qoi using VEF adj w Sn fwd: \t\t %g \n',results.delta_qoi_VEF_SNphi);
%fprintf('delta qoi using VEF math adjoint alt: \t\t %g \n',results.delta_qoi_VEF_a_alt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
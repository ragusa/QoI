function displaySensitivityError
global npar dat snq IO_opts results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n-----BEGIN DATA ON SENSITIVITY ERROR----- \n')
fprintf('This section contains data on the error in the sensitivity calculations. \n')
fprintf('In general, the "two forward solves" method will be used as the "true" value. \n')
fprintf('This means that there is a "true" sensitivity for two Sn solves and another \n')
fprintf('From two forward VEF solves. \n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('VEF adjoint error from SN forward: \t\t %g \n',(results.qoi_sn_f_pert - results.qoi_sn_f)-results.delta_qoi_VEF_a);
fprintf('VEF adjoint error from VEF forward: \t %g \n',(results.qoi_VEF_f_pert - results.qoi_VEF_f)-results.delta_qoi_VEF_a);
fprintf('VEF adjoint error from SN adjoint: \t\t %g \n',results.delta_qoi_sn_a-results.delta_qoi_VEF_a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
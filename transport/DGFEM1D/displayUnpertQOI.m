function displayUnpertQOI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n-----BEGIN UNPERTURBED QOI DATA OUTPUT----- \n')
global npar dat snq IO_opts results
fprintf('qoi using sn forward: \t %g \n',results.qoi_sn_f);
fprintf('qoi using sn adjoint: \t %g \n',results.qoi_sn_a);
fprintf('qoi using VEF forward: \t %g \n',results.qoi_VEF_f);
fprintf('qoi using VEF adjoint: \t %g \n',results.qoi_VEF_a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
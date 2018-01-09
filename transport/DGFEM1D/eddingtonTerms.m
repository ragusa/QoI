function eddingtonTerms
global npar dat snq IO_opts results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n-----BEGIN DATA ON EDDINGTON PERTURBATION----- \n')
fprintf('This section contains data on the terms resulting from the \n')
fprintf('perturbed Eddington and boundary Eddington. \n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[deltaE_term,deltaB_L,deltaB_R,volValues]=compute_deltaE_QoI_term(results.phiVEFa,results.phiVEF,results.phiVEF_pert,results.E,results.E_pert,results.Ebd,results.Ebd_pert);
fprintf('Total delta E+B term: \t %g \n',deltaE_term+deltaB_L+deltaB_R);
fprintf('-deltaE term: \t\t\t %g \n',deltaE_term);
fprintf('-deltaB term: \t\t\t %g \n',deltaB_L+deltaB_R);
fprintf('--deltaB_L term: \t\t %g \n',deltaB_L);
fprintf('--deltaB_R term: \t\t %g \n',deltaB_R);
%Plots the contribution of the volumetric (delta E) term at each element.
%Crude plot, needs refined.
figure(70)
plot(volValues)
end
function otherStuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n-----BEGIN DATA ON OTHER/under development METRICS---- \n')
fprintf('This section contains data on other metrics. \n')
fprintf('A bit of a holding pen for exploratory values/development. \n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve forward VEF problem using IP
[phiVEFa_Epert]=solve_VEF_math_adjoint(dat.adjoint_flux,E_pert,Ebd_pert);
do_plot(phiVEFa_Epert,'VEF Using E_{pert}',100,dat.adjoint_flux)
do_plot(phiVEFa_Epert-phiVEFa,'VEF_{diff}',120,dat.adjoint_flux)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rel_L1_diff=find_Eddington_diff(E,E_pert);
delta_qoi_VEF_math_a_Epert = compute_perturbed_qoi_VEF(dat.adjoint_flux,phiVEFa_Epert,phiVEF_pert,E_pert);
fprintf('delta qoi using VEF math adjoint and Epert: \t\t %g \n',delta_qoi_VEF_math_a_Epert);
fprintf('Diff from using Epert: \t\t %g \n',delta_qoi_VEF_a-delta_qoi_VEF_math_a_Epert);
fprintf('relative L1 difference in E: \t\t %g \n',Rel_L1_diff);
%adjDiff = compute_qoi(dat.adjoint_flux,phiVEFa_Epert,~sn,[],[])-compute_qoi(dat.adjoint_flux,phiVEFa,~sn,[],[]);
%fprintf('NEED TO USE PERT FIX SOURCE USELESS. Diff from using two VEF adjoints: \t\t %g \n',adjDiff);
end
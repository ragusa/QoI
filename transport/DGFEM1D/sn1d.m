function sn1d
% Linear Discontinous FEM code for Sn transport in 1D
% Jean Ragusa,
close all; clc; % closes all plotting figures, clears console

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global npar dat snq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select angular approx (must be an even number)
sn=8;
n_moments=1;
logi_galerkin=false;
[M,D,omega,weights] = gauss_legendre_1d(sn,n_moments,logi_galerkin);
snq.D = D; snq.M = M; snq.n_dir = sn;
snq.mu = omega; snq.w = weights;
% sum of the weights
snq.sw = sum(snq.w);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select data problem
pb_ID=14;
load_input(pb_ID);
console_io = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No Perturbations
dat.sourcePert = 0.0;
dat.sigsPert = dat.sigs*0;
dat.sigtPert = dat.sigt*0;
dat.sigaPert = dat.sigtPert - dat.sigsPert;
dat.inc_forward_pert = 0*dat.inc_forward;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve forward transport problem using sweeps
forward = true;     sn=true;    do_dsa = true;
[phi,E,Ebd,psi]=solve_transport(forward,do_dsa,console_io);
qoi_sn_f = compute_qoi_forward(phi,sn);
do_plot(phi,0,forward,E)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve adjoint transport problem using sweeps
forward = false;    sn=true;    do_dsa = true;
[phia,Ea,Ebda,psia]=solve_transport(forward,do_dsa,console_io);
qoi_sn_a = compute_qoi_adjoint(phia,sn,0,phi,E,psi,psia);
do_plot(phia,10,forward,Ea)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve forward VEF problem using IP
forward = true;     sn=false;
[phiVEF]=solve_VEF(forward,E,Ebd);
qoi_vef_f = compute_qoi_forward(phiVEF,sn);
do_plot(phiVEF,20,forward)
%do_plot((phiVEF./reshape(phi,npar.ndofs,1)-1),22)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve adjoint VEF problem using IP
forward = false;    sn=false;
[phiVEFa]=solve_VEF(forward,Ea,Ebda);
qoi_vef_a = compute_qoi_adjoint(phiVEFa,sn,0,phiVEF,E,psi,psia);
do_plot(phiVEFa,30,forward)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Perturbations. Used in adjoint sensitivity
dat.sourcePert =[1 1 1 1 1]*0;
dat.sigsPert = +[0.1 0.1 0.1 0.1 0.1].*dat.sigs*0;
dat.sigtPert = -[0.1 0.1 0.1 0.1 0.1].*dat.sigt*0;
dat.sigaPert = dat.sigtPert - dat.sigsPert;
dat.inc_forward_pert = 0.2*dat.inc_forward;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
forward=false; sn=true;
qoi_sn_a_pert = compute_qoi_adjoint(phia,sn,1,phi,E,psi,psia);
qoi_sn_a_delta = compute_qoi_adjoint(phia,sn,2,phi,E,psi,psia);
sn=false;
qoi_vef_a_pert = compute_qoi_adjoint(phiVEFa,sn,1,phiVEF,E,psi,psia);
qoi_vef_a_delta = compute_qoi_adjoint(phiVEFa,sn,2,phiVEF,E,psi,psia);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Apply perturbation to input values.
dat_saved = dat;
dat.qv_forward = dat.qv_forward + dat.sourcePert;
dat.sigs = dat.sigs + dat.sigsPert;
dat.sigt = dat.sigt + dat.sigtPert;
dat.siga = dat.siga + dat.sigaPert;
dat.inc_forward = dat.inc_forward + dat.inc_forward_pert;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve perturbed forward transport problem using sweeps
forward = true;     sn=true;    do_dsa = true;
[phi_pert,E_pert,Ebd_pert,psi_pert]=solve_transport(forward,do_dsa,console_io);
qoi_sn_f_pert = compute_qoi_forward(phi_pert,sn);
do_plot(phi_pert,40,forward,E)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve perturbed forward VEF problem using IP
forward = true;     sn=false;
[phiVEF_pert]=solve_VEF(forward,E_pert,Ebd_pert);
qoi_vef_f_pert = compute_qoi_forward(phiVEF_pert,sn);
do_plot(phiVEF_pert,50,forward)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write QoI Values
fprintf('---UNPERTURBED QOIS--- \n')
fprintf('qoi using sn forward: \t %g \n',qoi_sn_f);
fprintf('qoi using sn adjoint: \t %g \n',qoi_sn_a);
fprintf('qoi using VEFforward: \t %g \n',qoi_vef_f);
fprintf('qoi using VEFadjoint: \t %g \n',qoi_vef_a);
fprintf('---PERTURBED QOIS USING ADJOINT--- \n')
fprintf('perturbed qoi using sn adjoint: \t %g \n',qoi_sn_a_pert);
fprintf('perturbed qoi using VEFadjoint: \t %g \n',qoi_vef_a_pert);
fprintf('---DELTA QOI USING ADJOINT--- \n')
fprintf('delta qoi using sn adjoint: \t %g \n',qoi_sn_a_delta);
fprintf('delta qoi using VEFadjoint: \t %g \n',qoi_vef_a_delta);
fprintf('---PERTURBED QOIS USING FORWARD SOLVE--- \n')
fprintf('perturbed qoi using sn forward: \t %g \n',qoi_sn_f_pert);
fprintf('perturbed qoi using VEFforward: \t %g \n',qoi_vef_f_pert);
fprintf('---TEST LINE END OF OUTPUT--- \n')

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% END OF ROUTINE - below are code scraps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % forward = true; sn=true;
% % qoi_sn_f = compute_qoi(~forward,phi,sn);
% % fprintf('qoi using sn forward: \t %g \n',qoi_sn_f);
% % %
% % forward=false;
% % qoi_sn_a = compute_qoi(~forward,phia,sn);
% % fprintf('qoi using sn adjoint: \t %g \n',qoi_sn_a);
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % sn=false;
% % forward=true;
% % qoi_vef_f = compute_qoi(~forward,phiVEF,sn);
% % fprintf('qoi using VEFforward: \t %g \n',qoi_vef_f);
% % %
% % forward=false;
% % qoi_vef_a = compute_qoi(~forward,phiVEFa,sn);
% % fprintf('qoi using VEFadjoint: \t %g \n',qoi_vef_a);
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % solve forward VEF problem using IP
% % forward = false;
% % [phiVEFmath]=solve_VEF_math_adjoint(forward,E,Ebd);
% % % pretty plots
% % do_plot(phiVEFmath,300,forward)
% % qoi_vef_math_adj = compute_qoi(~forward,phiVEFmath,~sn);
% % fprintf('qoi using VEFmathadj: \t %g \n',qoi_vef_math_adj);
% % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % solve forward transport problem using sweeps
% forward = true;
% do_dsa = true;
% dat_saved = dat;
% dat.qv_forward = dat.qv_forward + dat.sourcePert;
% dat.sigs = dat.sigs + dat.sigsPert;
% dat.sigt = dat.sigt + dat.sigtPert;
% dat.siga = dat.siga + dat.sigaPert;
% dat.inc_forward = dat.inc_forward + dat.psiIncPert;
% [phi_pert,E_pert,Ebd_pert,psi_pert]=solve_transport(forward,do_dsa,console_io);
% % pretty plots
% do_plot(phi_pert,10,forward,E_pert)
% % reset
% dat = dat_saved;
% sn=true;
% qoi_sn_f_pert = compute_qoi(~forward,phi_pert,sn,0);
% fprintf('delta qoi using 2 sn forward runs: \t %g \n',qoi_sn_f_pert-qoi_sn_f);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% %Compute perturbed QoIs using adjoint method and unperturbed forward
% forward=false;
% delta_qoi_sn_a = compute_perturbed_qoi_Sn(~forward,phia,phi,E,psi,psia,sn);
% fprintf('delta qoi using sn adjoint: \t\t %g \n',delta_qoi_sn_a);
% 
% %error('ddd');
% 
% % qoi_vef_a = compute_perturbed_qoi_VEF(~forward,phiVEFa,phiVEF,E,Ebd);
% % fprintf('perturbed qoi using VEFadjoint: \t %g \n',qoi_vef_a);
% % qoi_vef_math_adj = compete_purturbed_qoi_VEF(~forward,phiVEFmath,phiVEF,E,Ebd);
% % fprintf('perturbed qoi using VEFmathadj: \t %g \n',qoi_vef_math_adj);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %Perturb data values for forward solve
% % dat.sigs = dat.sigs+dat.sigsPert;
% % dat.siga = dat.siga+dat.sigaPert;
% % dat.sigt = dat.sigt+dat.sigtPert;
% % dat.cdif = 1./(3*dat.sigt);
% % dat.qv_forward=dat.qv_forward+dat.sourcePert;
% % dat.inc_forward = dat.inc_forward + dat.psiIncPert ;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % solve forward transport problem using sweeps
% forward = true;
% do_dsa = true;
% [phip,Ep,Ebdp,psip]=solve_transport(forward,do_dsa,console_io);
% Ebdp;
% % pretty plots
% do_plot(phip,0,forward)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % solve forward VEF problem using IP
% forward = true;
% [phiVEFp]=solve_VEF(forward,Ep,Ebdp);
% % pretty plots
% %do_plot(phiVEFp,0,forward)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% forward = true; sn=true;
% qoi_sn_f = compute_qoi(~forward,phip,sn,0);
% fprintf('perturbed qoi using sn forward: \t %g \n',qoi_sn_f);
% sn=false;
% qoi_vef_f = compute_qoi(~forward,phiVEFp,sn,0);
% fprintf('perturbed qoi using VEFforward: \t %g \n',qoi_vef_f);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% return
% %%%%%%%%%%%Check 2nd order stuff
% load_input(pb_ID);
% forward=false;
% qoi_sn_a_second = compute_perturbed_qoi_Sn(~forward,phia,phip,Ep,psip,psia);
% fprintf('perturbed qoi using sn adjoint 2nd order: \t %g \n',qoi_sn_a_second);
% qoi_vef_a = compute_perturbed_qoi_VEF(~forward,phiVEFa,phiVEFp,Ep,Ebd);
% fprintf('perturbed qoi using VEFadjoint 2nd order: \t %g \n',qoi_vef_a);
% qoi_vef_math_adj = compute_perturbed_qoi_VEF(~forward,phiVEFmath,phiVEFp,Ep,Ebd);
% fprintf('perturbed qoi using VEFmathadj 2nd order: \t %g \n',qoi_vef_math_adj);
% return
% end

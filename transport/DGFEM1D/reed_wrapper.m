function [Y]=reed_wrapper(X)
global npar dat snq IO_opts results
nominal_input=dat.nominal_input;
load_reed(nominal_input);
[results.phi,results.E,results.Ebd,results.psi]=solveForwardSn;
nom_phi=results.phi;
nom_E=results.E;
nom_Ebd=results.Ebd;
nom_psi=results.psi;
%do_plot(nom_phi,'$$\phi_{nom}$$',1,dat.forward_flux)
%qoi = compute_qoi(dat.forward_flux,results.phi,snq.n_dir,results.psi,results.psi);
%fprintf('Nominal: %d \n',qoi);


nSamp=length(X);
Y =zeros(nSamp,5+800);
for ii=1:nSamp
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load Reed Problem
    if dat.use_reduced_input
        test_input=nominal_input;
        test_input(4)=X(ii,1);
        test_input(7)=X(ii,2);
    else
        test_input=X(ii,:);
    end
    load_reed(test_input);
    [results.phi_p,E_p,results.Ebd_p,results.psi_p]=solveForwardSn;
    [results.phiVEF_p]=solveForwardVEF;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Sore values of interest
    dat.qv_adjoint=[0 0 0 0 1];
    qoi5_sn = compute_qoi(dat.forward_flux,results.phi_p,snq.n_dir,results.psi_p,results.psi_p);
    qoi5_VET = compute_qoi(dat.forward_flux,results.phiVEF_p,~snq.n_dir,[],[]);
    dat.qv_adjoint=[0 0 1 0 0];
    qoi3_sn = compute_qoi(dat.forward_flux,results.phi_p,snq.n_dir,results.psi_p,results.psi_p);
    qoi3_VET = compute_qoi(dat.forward_flux,results.phiVEF_p,~snq.n_dir,[],[]);
    dE=nom_E-E_p;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Progress: %d of %d:  \t qoi_sn= %g  \n qoi_VEF= %g  \n',ii,nSamp,qoi5_sn,qoi5_VET); 
    Y(ii,1) =  qoi5_sn;
    Y(ii,2) =  qoi5_VET;
    Y(ii,3) =  qoi3_sn;
    Y(ii,4) =  qoi3_VET;
    Y(ii,5) =  norm(dE);
    Y(ii,6:805)=reshape(E_p,800,1);
end
end


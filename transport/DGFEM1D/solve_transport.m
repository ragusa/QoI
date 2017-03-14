function [varargout]=solve_transport(forward,do_dsa)

global dat npar snq

if do_dsa
    dsa_method = 'MIP';
    dsa_bc_type = 'Robin';
    alpha_sn = false;
    % build DSA matrix
    switch dsa_bc_type
        case 'Neumann'
            error('Neumann bc not accepted for DSA study')
        case 'Robin'
            dat.bc.left.type = 1; % 1=Robin
            dat.bc.rite.type = 1; % 1=Robin
        case 'Dirichlet'
            dat.bc.left.type = 2; % 2=Dirichlet
            dat.bc.rite.type = 2; % 2=Dirichlet
        otherwise
            error('Unknown bc type %s in %s',bc_type,mfilename);
    end
    if alpha_sn
        npar.alpha = dot(weights(sn/2+1:end),omega(sn/2+1:end)') / snq.sw;
    else
        npar.alpha = 0.25;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    npar.lump = false;
    npar.Ckappa = 10;
    npar.Ckappa_bd=2*npar.Ckappa;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [A,Sd] = build_dsa_matrix(dsa_method);
end

% iteration parameters
npar.max_SI_iter  = 5000;
npar.SI_tolerance = 1e-8;
keep_angular_flux = false;

% initializa scalar flux
phi_new = zeros(npar.porder+1,npar.nel);
% compute total src
q = ext_scatter_src(phi_new,forward);


%  begin Source Iteration Loop
for iter=1:npar.max_SI_iter,
    % save previous scalar flux
    phi_old = phi_new;
    % perform transport sweep
    phi_new = sweep1D_LD(q, forward, keep_angular_flux);
    % perform dsa if requested
    if do_dsa
        rhs = Sd*reshape(phi_new-phi_old,npar.ndofs,1);
        dphi = A\rhs;
        phi_new = phi_new + reshape(dphi,npar.porder+1,npar.nel);
    end
    % compute error norm
    new_norm = norm(reshape(phi_new-phi_old,npar.ndofs,1),2);
    % check for convergence, otherwise prepare next SI
    if new_norm<npar.SI_tolerance
        break
    else
        % update RHS term
        q = ext_scatter_src(phi_new,forward);
        if(iter==npar.max_SI_iter), warning('about to exit w/o convergence'); end
    end
    % console printouts
    fprintf('SI iteration %4i, error %7e',iter,new_norm);
    if iter>1
        fprintf(', NSR %7e',new_norm/old_norm);
    end
    old_norm = new_norm;
    fprintf('\n');
end


% output arguments
varargout{1} = phi_new;
if nargout==4
    keep_angular_flux=true;
    [phi_new,psi] = sweep1D_LD(q, forward, keep_angular_flux);
    varargout{4} = psi;
    [E,Ebd] = compute_Eddington_tensor(phi_new,psi);
    varargout{2} = E;
    varargout{3} = Ebd;
end


end

function [varargout]=solve_transport(forward,do_dsa)

global dat npar snq

if do_dsa
    dsa_method = 'MIP';
    dsa_bc_type = 'Robin';
    % build DSA matrix
    switch dsa_bc_type
        case 'Neumann'
            error('Neumann bc non accepted for DSA study')
        case 'Robin'
            dat.bc.left.type = 1; % 1=Robin
            dat.bc.rite.type = 1; % 1=Robin
        case 'Dirichlet'
            dat.bc.left.type = 2; % 2=Dirichlet
            dat.bc.rite.type = 2; % 2=Dirichlet
        otherwise
            error('Unknown bc type %s in %s',bc_type,mfilename);
    end
    [A,Sd] = build_dsa_matrix(dsa_method);
end

% iteration parameters
npar.max_SI_iter  = 500;
npar.SI_tolerance = 1e-8;
keep_angular_flux = false;

% forward or adjoint
forward=true;
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
    % compute error norm
    err = norm(phi_new-phi_old,inf);
    fprintf('SI iteration %4i, error %7e \n',iter,err);
    % check for convergence, otherwise prepare next SI
    if err<npar.SI_tolerance
        break
    else
        % update RHS term
        q = ext_scatter_src(phi_new,forward); 
        if(iter==npar.max_SI_iter), warning('about to exit w/o convergence'); end
    end
end


% output arguments
varargout{1} = phi_new;
if nargout==2
    keep_angular_flux=true;
    [phi_new,psi] = sweep1D_LD(q, forward, keep_angular_flux);
    varargout{2} = psi;
end


end
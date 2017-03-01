function sn1d
% Jean Ragusa, 2012, Prague
close all; clc; % closes all plotting figures, clears console
spa_disc=input('Spatial disc: Enter 0 (SD), 1(DD), 2(LD) ');
switch spa_disc
    case{0,1}
        nukn_cell=1;
    case{2}
        nukn_cell=2;
    otherwise
        error('Wrong spa_disc %i has been entered. Only 0,1,2 are valid',spa_disc);
end
% load material, domain and source properties
[x,pb_data]=load_input(input('material zone refinement (>0) '));
% load the angular quadrature
snq=load_quadrature(input('Enter quadrature (2,4,6,8) '));

% initialize angular flux (ncell, ndir, nunknowns)
ncells=length(x)-1;
psi=zeros(ncells,snq.sn,nukn_cell); % IMPORTANT: psi NEED NOT be kept (only for plotting in this demo)
phi=zeros(ncells,nukn_cell);    % initialize scalar flux
q  =zeros(ncells,nukn_cell);    % initialize total src (q = external src + scatter src)
% compute initial q
[q] = ext_scatter_src(q, phi, ncells, pb_data);

% iteration parameters
maxiter = 500;
tol     = 1e-8;

%  begin Source Iteration Loop
for iter=1:maxiter,
    phi_old = phi;   % save previous scalar flux
    % perform transport sweep, AGAIN, psi NEED NOT be kept
    [psi] = sweep1D(psi, snq, pb_data, q, x, spa_disc);
    % update scalar flux (USUALLY done witin sweeps. Do you understand why?)
    [phi] = scalar_flx(psi, phi, snq, ncells);
    % compute error norm
    err = norm(phi-phi_old,inf);
    disp(sprintf('iteration %4i, error %7e',iter,err));
    % check for convergence, otherwise prepare next SI
    if(err<tol)
        break
    else
        % update RHS term
        [q] = ext_scatter_src(q, phi, ncells, pb_data); % compute q
        if(iter==maxiter), error('about to exit w/o convergence'); end
    end
end
% pretty plots
do_plot(phi,psi,x,spa_disc)

return
end
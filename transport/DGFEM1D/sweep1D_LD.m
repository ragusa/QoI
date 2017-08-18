function [varargout]=sweep1D_LD(q,forward,keep_angular_flux)
% perform LD sweep
% in: total source (scattering+external)
%     whether to do farward or adjoint transport solve
%     whether to keep the angular flux or not

global snq npar dat

% sanity check
if keep_angular_flux && nargout~=2
    error('To keep the angular flux, you need to have 2 output arguments in %s',mfilename);
end
% storage for the new scalar flux
phi = zeros(npar.porder+1,npar.nel);
% storage to keep angular fluxes if requested
if keep_angular_flux
    psi = zeros(npar.porder+1,npar.nel,snq.n_dir);
end

% source data (incoming) 
if forward
    inc = dat.inc_forward;
else
    inc = dat.inc_adjoint;
end
% local mass/grad/edge matrices
mass = [2 1;1 2]/6;
grad = [1 1;-1 -1]/2;
edgp = [0 0;0 1];
edgm = [1 0;0 0];

% loop over directions
for idir=1:snq.n_dir;
    
    % choose the sweeping order (L to R or R to L)
    if(snq.mu(idir)>0)
        ibeg=1;  iend=npar.nel;  incr=+1;
    else
        ibeg=npar.nel;  iend=1;  incr=-1;
    end
    % retrieve boundary condition for that direction
    psi_in = inc(idir); 
    % shortcut
    mu = snq.mu(idir);
    gradient = grad*mu;
    
    % loop over spatial cells
    for i=ibeg:incr:iend
        my_zone=npar.iel2zon(i);
        sig = dat.sigt(my_zone);     %IWH look at this angular vs non-angular sigma. Add /snq.sw ?;
        
        dx   = npar.dx(i);
        % compute L/R psi for given cell i
        rhs = mass * dx * q(:,i);
        L = gradient + mass * dx * sig;
        if (mu>0)
            rhs(1) = rhs(1) + psi_in*mu;
            L      = L      + edgp*mu;
        else
            rhs(2) = rhs(2) - psi_in*mu;
            L      = L      - edgm*mu;
        end
        % psi_out gives array of psi_R and psi_L (in order) of cell i
        psi_local = L\rhs;
        % prepare for next cell solve
        if mu>0
            psi_in = psi_local(2);
        else
            psi_in = psi_local(1);
        end
        % compute the new scalar fluxes
        phi(:,i) = phi(:,i) + snq.w(idir)*psi_local;
        % storage angular flux
        if keep_angular_flux 
            psi(:,i,idir) = psi_local;
        end
    end
end

% output arguments
varargout{1} = phi;
if keep_angular_flux
   varargout{2} = psi;
end

end

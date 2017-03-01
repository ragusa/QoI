function [psi]=sweep1D_LD( psi, mu, xs, q, z2med, x);
% perform LD sweep
ndir=length(mu);
% loop over directions
for idir=1:ndir
    % choose the sweeping order (L to R or R to L)
    if(mu(idir)>0)
        ibeg=1;  iend=length(x)-1;  incr=+1;
    else
        ibeg=length(x)-1;  iend=1;  incr=-1;
    end    
    % retrieve boundary condition for that direction
    psi_in=0. ; % imposed vaccum BC, can be improved the generality later   
    m=mu(idir); % shortcut 
    % loop over spatial cells
    for i=ibeg:incr:iend
        imed = z2med(i);
        sig  = xs(imed); 
        dx=x(i+1)-x(i);
        % compute L/R psi for given cell i
        % the source q is constant in each zone QR=QL
        if (mu(idir)>0)
            RHS=[(psi_in*m+(q(i,1)/3+q(i,2)/6)*dx);((q(i,1)/6+q(i,2)/3)*dx)];
            L=[(m/2+sig*dx/3) (m/2+sig*dx/6); (-m/2+sig*dx/6) (-m/2+m+sig*dx/3)];
        else
            RHS=[(dx*(q(i,1)/3+q(i,2)/6));(dx*(q(i,1)/6+q(i,2)/3)-m*psi_in)];
            L=[(m/2-m+sig*dx/3) (m/2+sig*dx/6); (-m/2+sig*dx/6) (-m/2+sig*dx/3)];
        end
        % psi_out gives array of psi_R and psi_L (in order) of cell i
        psi_sol=L\RHS;
        % save right and left psi for cell i for a given direction (only for plotting!!!)
        psi(i,idir,:)=psi_sol;
        % prepare for next cell
        if (mu(idir)>0)
            psi_in=psi_sol(2);
        else
            psi_in=psi_sol(1);
        end
    end
end
return
end

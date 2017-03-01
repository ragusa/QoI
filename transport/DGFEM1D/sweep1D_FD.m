function [psi]=sweep1D_FD( psi, mu, xs, q, z2med, x, a);
% perform FD sweep
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
    m=abs(mu(idir)); % shortcut 
    % loop over spatial cells
    for i=ibeg:incr:iend        
        imed = z2med(i);
        sig  = xs(imed); 
        dx = x(i+1)-x(i);       
        % compute exiting psi
        psi_out = ( q(i) + ( m/dx - (1-a)*sig ) * psi_in ) / ( m/dx + a*sig ) ;        
        % store cell-average psi (only for plotting!!!)
        psi(i,idir) = (1-a)*psi_in + a*psi_out;      
        % prepare for next cell
        psi_in=psi_out;
    end
end
return
end
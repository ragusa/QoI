function [psi]=sweep1D(psi, snq, pb, q, mesh, spa_disc);
% select Finite difference (SD or DD) sweeps or LD sweep
if(spa_disc~=2)
    a=1; % SC
    if spa_disc==1 , a=0.5; end % DD
    [psi]=sweep1D_FD( psi, snq.mu, pb.XS.tot, q, pb.ID.mat, mesh, a);
else
    [psi]=sweep1D_LD( psi, snq.mu, pb.XS.tot, q, pb.ID.mat, mesh );
end
return 
end

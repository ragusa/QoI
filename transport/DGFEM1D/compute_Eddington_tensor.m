function [E,Ebd] = compute_Eddington_tensor(phi,psi)

global npar snq

E = zeros(npar.porder+1,npar.nel);

for iel=1:npar.nel
    
    for idir=1:snq.n_dir
        E(:,iel) = E(:,iel) + snq.w(idir)*snq.mu(idir)^2*psi(:,iel,idir);
    end
    E(:,iel) = E(:,iel)./phi(:,iel);
    
end

Ebd = zeros(2,1);

for idir=1:snq.n_dir
    Ebd(1,1) = Ebd(1,1) + snq.w(idir)*abs(snq.mu(idir))*psi(1,1,idir);
    Ebd(2,1) = Ebd(2,1) + snq.w(idir)*abs(snq.mu(idir))*psi(end,end,idir);
end
Ebd(1,1) = Ebd(1,1) / phi(1,1);
Ebd(2,1) = Ebd(2,1) / phi(2,end);
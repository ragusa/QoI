function E = compute_Eddington_tensor(phi,psi)

global npar snq

E = zeros(npar.porder+1,npar.nel);

for iel=1:npar.nel
    
    for idir=1:snq.n_dir
        E(:,iel) = E(:,iel) + snq.w(idir)*snq.mu(idir)^2*psi(:,iel,idir);
    end
    E(:,iel) = E(:,iel)./phi(:,iel);
    
end
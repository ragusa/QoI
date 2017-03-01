function [phi]=scalarflx(psi, phi, SNQ, ncells);
% compute scalar flux
phi(:,:)=0.; % reset scalar flux
for i=1:ncells,
    for idir=1:SNQ.sn,
        phi(i,:) = phi(i,:)+ SNQ.w(idir) * (shiftdim(psi(i,idir,:),1));
    end
end
return
end
        
function [q]=ext_scatter_src(phi,forward)
% compute the total source (ext+scattering) 

global dat npar snq

% source data (volumetric) 
if forward
    qv  = dat.qv_forward;
else
    qv  = dat.qv_adjoint;
end

% initialize
q = zeros(npar.porder+1,npar.nel);

% loop over each spatial cell
for iel=1:npar.nel

    my_zone=npar.iel2zon(iel);
    sigs = dat.sigs(my_zone);
    qext = qv(my_zone);

    q(:,iel) = qext + sigs * phi(:,iel);
end

% isotropic source and scattering
q = q/snq.sw; 

end
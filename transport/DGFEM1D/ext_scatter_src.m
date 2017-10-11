function [q]=ext_scatter_src(phi,forward,altMethod)
% compute the total source (ext+scattering) 

global dat npar snq

% source data (volumetric) 
if forward
    % isotropic source. we need the angular dependent SRD in the forward
    % problem
    qv  = dat.qv_forward / snq.sw;
else
    % we want q^\dagger(x,mu) = response function(x) given in load_input so
    % no change here
    qv   = dat.qv_adjoint;
end

% initialize
q = zeros(npar.porder+1,npar.nel);

% loop over each spatial cell
for iel=1:npar.nel

    my_zone=npar.iel2zon(iel);
    % isotropic scattering
    sigs = dat.sigs(my_zone)/snq.sw;
    qext = qv(my_zone);

    q(:,iel) = qext + sigs * phi(:,iel);
end

end
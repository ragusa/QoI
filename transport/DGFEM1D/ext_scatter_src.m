function [q]=ext_scatter_src(q, phi, ncells, pb);
% compute the total source (ext+scattering) 
q(:,:)=0;
% loop over each spatial cell
for i=1:ncells,
    imed = pb.ID.mat(i); % retrieve material ID
    xs   = pb.XS.sca(imed); % retrieve scattering XS value
    isrc = pb.ID.src(i);    % retrieve ext src ID
    sext=0;
    if(isrc>0), sext = pb.src(isrc); end
    q(i,:) = (sext + xs * phi(i,:) )/ 2;
end
return
end
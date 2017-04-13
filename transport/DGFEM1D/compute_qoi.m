function qoi = compute_qoi(forward,phi,is_sn)

global npar dat snq

% source data (volumetric)
if forward
    qv  = dat.qv_forward;
    % qv is the space-dependent src rate density -SRD- [part/cm^3-s]
    % in the Sn equation, it is  made into an (isotropic) angular-dependent
    % SRD [part/cm^3-s-unit_cosine] by dividing by the sum of the angular
    % weights (2)
    %
    % when one applies the inner product, for Sn, it is the angle-dependent
    % terms that are in the definition:
    % QoI = \int dmu qv(mu) psi(mu) = \int dmu qv/2 psi(mu) = qv/2 phi
    % (note there is everywhere an integral of space dx that I omitted for
    % brevity. it plays no role in the demonstration
    %
    % hence below, we need to do the division by 2 (sw)
    if is_sn
        qv = qv / snq.sw;
    end
else
    % this is q^\dagger(x,mu) = function given as input in load input. no need to tweak
    qv  = dat.qv_adjoint;
end


% type of solution data passed in
[n1,n2]=size(phi);
if n2==1 % if DFEM solution store as vector, reshape it into a 2d array (local_dofs, elem)
    phi=reshape(phi,npar.porder+1,npar.nel);
end

% shortcuts
porder= npar.porder;
% quadrature
[xq,wq] = GLNodeWt(porder+1);
% store shapeset
[b,~] =feshpln(xq,porder);
% compute elementary mass matrix
m=zeros(porder+1,porder+1);
for i=1:porder+1
    for j=1:porder+1
        m(i,j)= dot(wq.*b(:,i), b(:,j));
    end
end

qoi=0;
% loop over elements
for iel=1:npar.nel
    % element extremities
    % x0 = npar.x(iel);
    % x1 = npar.x(iel+1);
    % get x values in the interval
    % x = (x1+x0)/2+xq*(x1-x0)/2;
    % jacobian of the transformation to the ref. element
    Jac = npar.dx(iel)/2;
    my_zone=npar.iel2zon(iel);
    qext = qv(my_zone);
    % compute local matrices + load vector
    % assemble
    qoi = qoi + Jac*qext*dot(ones(porder+1,1),m*phi(:,iel));
end

% if ~forward
%     % Add boundary terms if using the adjoitn to compute the QoI
%     dir_index_rite = 1:snq.n_dir/2; % the first half of the directions are <0
%     dir_index_left = snq.n_dir/2+1:snq.n_dir; % the second half of the directions are >0
%     
%     % Get angular flux on right
%     psi_rite = shiftdim(psi(npar.porder+1,npar.nel,:),1);
%     % overwrite with BC values
%     psi_rite(index_rite) = dat.inc_forward(dir_index_rite);
% 
%     % Get angular flux on left
%     psi_left = shiftdim(psi(1,1,:),1);
%     % overwrite with BC values
%     psi_left(dir_index_left) = dat.inc_forward(dir_indexdir_index_left_rite);
% 
%     % the reason for using the "other" direction for the adjoint flux is that
%     % psia(mu) = psi(-mu), or equivalently psi(-mu)=psi(mu) and remember we
%     % faked a "forward" solve for the adjoint
%     
%     % Incident flux on left
%     iel=1;
%     % total incident
%     incident_left = dat.inc_forward(dir_index_left);
% %     adjoint_flx_left = shiftdim(psia(1,iel,dir_index_rite),1);
%     %
%     % qoi = qoi - dot(snq.w.*[incident_left incident_rite]',[adjoint_flx_left adjoint_flx_rite]);
% end

function qoi = compute_qoi(forward,phi)

global npar dat

% source data (volumetric) 
if forward
    qv  = dat.qv_forward;
else
    qv  = dat.qv_adjoint;
end

% type of solution data passed in
[n1,n2]=size(phi);
if n2==1
    phi=reshape(phi,npar.porder+1,npar.nel);
end

% shortcuts
porder= npar.porder;
% quadrature
[xq,wq] = GLNodeWt(porder+1);
% store shapeset
[b,dbdx] =feshpln(xq,porder);

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
    siga = dat.siga(my_zone);
    qext = qv(my_zone);
    % compute local matrices + load vector
    for i=1:porder+1
        f(i)= dot(qext.*wq, b(:,i));
    end
    % assemble
    qoi = qoi + Jac*dot(f,phi(:,iel));
end

if ~forward
    % Add boundary terms if using the adjoint to compute the QoI
    neg_dir = 1:snq.n_dir/2; % the first half of the directions are <0
    pos_dir = snq.n_dir/2+1:snq.n_dir; % the second half of the directions are >0
    
    % Get angular flux at right extremity
    psi_rite = shiftdim(psi(npar.porder+1,npar.nel,:),1);
    % overwrite with BC values
    psi_rite(neg_dir) = dat.inc_forward(neg_dir);

    % Get angular flux at left extremity
    psi_left = shiftdim(psi(1,1,:),1);
    % overwrite with BC values
    psi_left(pos_dir) = dat.inc_forward(pos_dir);

    % the reason for using the "other" direction for the adjoint flux is that
    % psia(mu) = psi(-mu), or equivalently psia(-mu)=psi(mu) and remember we
    % faked a "forward" solve for the adjoint
    reverse_dir = snq.n_dir:-1:1
    psia_rite = shiftdim(psi(npar.porder+1,npar.nel,reverse_dir),1);
    psia_left = shiftdim(psi(1,1,reverse_dir),1);
    % overwrite with BC values
    psia_bc = dat.inc_adjoint(reverse_dir);
    psia_left(neg_dir) = dat.inc_adjoint(neg_dir);
    psia_rite(pos_dir) = dat.inc_adjoint(pos_dir);
    % right extremity: vo.vn = vo.ex
    qoi = qoi - dot(snq.w.*snq.mu.*psi_rite,psia_rite);
    % left extremity: vo.vn = -vo.ex
    qoi = qoi + dot(snq.w.*snq.mu.*psi_left,psia_left);
end

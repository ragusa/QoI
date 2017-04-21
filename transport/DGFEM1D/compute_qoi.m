function qoi = compute_qoi(use_forward_flux,phi,psi,psia,is_sn)

global npar dat snq

% source data (volumetric)
if use_forward_flux
    % this is q^\dagger(x,mu) = function given as input in load input. no need to tweak
    qv  = dat.qv_adjoint;
else
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

fprintf('qoi before bc %g (forward=%g) \n',qoi,use_forward_flux);

if ~use_forward_flux
    %
    warning('here, we use both forward and adjoint angular fluxes for BC terms. If 0-adjoitn BC are used, we should be able to get away with only pisa as a given');
    
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
    reverse_dir = snq.n_dir:-1:1;
    psia_rite = shiftdim(psia(npar.porder+1,npar.nel,reverse_dir),1);
    psia_left = shiftdim(psia(1,1,reverse_dir),1);
    % overwrite with BC values
    %%% psia_bc = dat.inc_adjoint(reverse_dir);
    psia_left(neg_dir) = dat.inc_adjoint(neg_dir);
    psia_rite(pos_dir) = dat.inc_adjoint(pos_dir);
    % right extremity: vo.vn = vo.ex
    qoi = qoi - dot(snq.w'.*snq.mu.*psi_rite,psia_rite);
    [psi_rite ; psia_rite]';
    % left extremity: vo.vn = -vo.ex
    qoi = qoi + dot(snq.w'.*snq.mu.*psi_left,psia_left);
    [psi_left;psia_left]';
    %
    fprintf('qoi adjoint bc left = %g \n',-dot(snq.w'.*snq.mu.*psi_rite,psia_rite));
    fprintf('qoi adjoint bc rite = %g \n',dot(snq.w'.*snq.mu.*psi_left,psia_left));
end

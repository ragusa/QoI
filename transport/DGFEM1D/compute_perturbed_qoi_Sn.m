function d_qoi = compute_perturbed_qoi_Sn(use_forward_flux,phia,phi_unpert,psi_unpert,psia,is_sn)

global npar dat snq IO_opts

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
[n1,n2]=size(phia);
if n2==1
    phia=reshape(phia,npar.porder+1,npar.nel);
end
[n1,n2]=size(phi_unpert);
if n2==1
    phi_unpert=reshape(phi_unpert,npar.porder+1,npar.nel);
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
d_qoi=0;

% loop over elements
for iel=1:npar.nel
    my_zone=npar.iel2zon(iel);
    delta_sigs= dat.sigsPert(my_zone);
    delta_sigt= dat.sigtPert(my_zone);
    delta_qext = dat.sourcePert(my_zone);
    if is_sn,~use_forward_flux;
        delta_qext = delta_qext / snq.sw;
    end
    Jac   = npar.dx(iel)/2;
    % assemble
    d_qoi = d_qoi + Jac*dot(m*ones(2,1)*delta_qext,phia(:,iel)) ;
    d_qoi = d_qoi + Jac*delta_sigs/snq.sw*dot(m*phi_unpert(:,iel),phia(:,iel));
    % Anglular integration of psi psia product
    for idir=1:snq.n_dir
        backDir = snq.n_dir+1-idir; 
        d_qoi = d_qoi - delta_sigt*Jac* snq.w(idir)* dot(m*psi_unpert(:,iel,idir), psia(:,iel,backDir));
    end
end

if IO_opts.show_dqoi_pre_bc
    fprintf('dqoi before bc %g (forward=%g) \n',d_qoi,use_forward_flux);
end

%%%% d_qoi=d_qoi - dot(snq.w.*[incident_left incident_rite]',[adjoint_flx_left adjoint_flx_rite]);

% Add boundary terms if using the adjoint to compute the QoI
neg_dir = 1:snq.n_dir/2; % the first half of the directions are <0
pos_dir = snq.n_dir/2+1:snq.n_dir; % the second half of the directions are >0

% Set ang. flux at right extremity
% in the dQoI, we do not want to have to know the perturbed outgoing
% forward flux, so we need to ensure that the adjoitn boundary src is 0
psi_rite = zeros(1,snq.n_dir); %%% shiftdim(psi(npar.porder+1,npar.nel,:),1);
% overwrite with BC values
psi_rite(neg_dir) =  dat.psiIncPert(neg_dir);

% Set angular flux at left extremity
psi_left = zeros(1,snq.n_dir); %%% shiftdim(psi(1,1,:),1);
% overwrite with BC values
psi_left(pos_dir) =  dat.psiIncPert(pos_dir);

% the reason for using the "opposite" direction for the adjoint flux is that
% psia(mu) = psi(-mu), or equivalently psia(-mu)=psi(mu) and remember we
% faked a "forward" solve for the adjoint
reverse_dir = snq.n_dir:-1:1;
psia_rite = shiftdim(psia(npar.porder+1,npar.nel,reverse_dir),1);
psia_left = shiftdim(psia(1,1,reverse_dir),1);
% overwrite with BC values
psia_left(neg_dir) = dat.inc_adjoint(neg_dir);
psia_rite(pos_dir) = dat.inc_adjoint(pos_dir);
if sum(abs(psia_left(neg_dir)))>eps || sum(abs(psia_rite(pos_dir)))>0
    error('in %s, angular adjoint boundary src must be zero',mfilename);
end
% right extremity: vo.vn = vo.ex
d_qoi = d_qoi - dot(snq.w'.*snq.mu.*psi_rite,psia_rite);
% left extremity: vo.vn = -vo.ex
d_qoi = d_qoi + dot(snq.w'.*snq.mu.*psi_left,psia_left);

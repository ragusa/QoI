function [d_qoi ,parts] = compute_perturbed_qoi_VEF(use_forward_flux,phia,phi_unpert,E,snphia,snpsia)

global npar dat snq IO_opts

blended=0;
if nargin>4
    blended=1;
end
if nargin>5
    blended=2;
end

parts=zeros(4,1);

% source data (volumetric)
if use_forward_flux
    % this is q^\dagger(x,mu) = function given as input in load input. no need to tweak
    qv  = dat.qv_adjoint;
else
    qv  = dat.qv_forward;
    % qv is the space-dependent src rate density -SRD- [part/cm^3-s]
end

% type of solution data passed in
[n1,n2]=size(phia);
if n2==1
    phia=reshape(phia,npar.porder+1,npar.nel);
    if blended>0
        snphia=reshape(snphia,npar.porder+1,npar.nel);
    end
    if blended==2
        snphia=reshape(snphia,npar.porder+1,npar.nel);
    end
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
[b,dbdx] =feshpln(xq,porder);
% compute elementary mass matrix
m=zeros(porder+1,porder+1);
md=zeros(porder+1,porder+1);
mdd=zeros(porder+1,porder+1);
for i=1:porder+1
    for j=1:porder+1
        m(i,j)= dot(wq.*b(:,i), b(:,j));
        md(i,j)= dot(wq.*dbdx(:,i), b(:,j));
        mdd(i,j)= dot(wq.*dbdx(:,i), dbdx(:,j));
    end
end
d_qoi=0;
val=0;

% loop over elements
for iel=1:npar.nel
    my_zone=npar.iel2zon(iel);
    delta_siga= dat.sigaPert(my_zone);
    isigtr = 3*dat.cdif(my_zone); % inverse of sigma_tr 
    isigtrp= 3*1./(3*(dat.sigt(my_zone)+dat.sigtPert(my_zone)));
    delta_isigtr=isigtrp-isigtr;
    E0=E(1,iel);
    E1=E(2,iel);
    Eloc = (E1+E0)/2+xq*(E1-E0)/2;
    dEdx = (E1-E0)/2;
    delta_qext = dat.sourcePert(my_zone);
    Jac   = npar.dx(iel)/2;
    % assemble
    %Source Perturbation terms
    if blended>0
        val=Jac*dot(snphia(:,iel), m*ones(2,1)*delta_qext/ snq.sw);
    else
        val=Jac*dot(phia(:,iel), m*ones(2,1)*delta_qext); 
    end
    parts(1)=parts(1)+val;
    d_qoi = d_qoi + val;
    %Absorption dQoI  Terms
    val = - Jac*delta_siga*dot(phia(:,iel), m*phi_unpert(:,iel));
    parts(2)=parts(2)+val;
    d_qoi = d_qoi +val;
    %%Think I need to redo the below, once I figure out the
    for i=1:porder+1
        for j=1:porder+1
            k(i,j)= dot(wq.*dbdx(:,i) , Eloc.*dbdx(:,j)+dEdx*b(:,j));
        end
    end
    val= - delta_isigtr/Jac*dot(phia(:,iel), k*phi_unpert(:,iel));
    parts(3)=parts(3)+val;
    d_qoi = d_qoi +val;
end

if IO_opts.show_dqoi_pre_bc
    fprintf('dqoi before bc %g (forward=%g) \n',d_qoi,use_forward_flux);
end

if blended==2
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
    psia_rite = shiftdim(snpsia(npar.porder+1,npar.nel,reverse_dir),1);
    psia_left = shiftdim(snpsia(1,1,reverse_dir),1);
    % overwrite with BC values
    psia_left(neg_dir) = dat.inc_adjoint(neg_dir);
    psia_rite(pos_dir) = dat.inc_adjoint(pos_dir);
    if sum(abs(psia_left(neg_dir)))>eps || sum(abs(psia_rite(pos_dir)))>0
        error('in %s, angular adjoint boundary src must be zero',mfilename);
    end
    % right extremity: vo.vn = vo.ex
    parts(4)= parts(4) - dot(snq.w'.*snq.mu.*psi_rite,psia_rite);
    parts(4)=parts(4) + dot(snq.w'.*snq.mu.*psi_left,psia_left);
    d_qoi = d_qoi - dot(snq.w'.*snq.mu.*psi_rite,psia_rite);
    % left extremity: vo.vn = -vo.ex
    d_qoi = d_qoi + dot(snq.w'.*snq.mu.*psi_left,psia_left);
else
    %%% Boundary Condition Stuff
    % Build forward J values (mainly copied from establish_bc_for_VEF)
    inc=dat.psiIncPert;
    ndir = snq.n_dir;
    ii = find (inc(1:ndir/2)~=0);
    if isempty(ii) % vacuum
        JForwardRite = 0;
    else
        Jinc = dot (snq.w(1:ndir/2)'.*snq.mu(1:ndir/2), inc(1:ndir/2) );
        JForwardRite = -Jinc;
    end
    ii = find (inc(ndir/2+1:end)~=0);
    if isempty(ii) % vacuum
        JForwardLeft = 0;
    else
        Jinc = dot (snq.w(ndir/2+1:end)'.*snq.mu(ndir/2+1:end), inc(ndir/2+1:end) );
        JForwardLeft = Jinc;
    end
    BCqoiRite=2*JForwardRite*phia(npar.porder+1,npar.nel);%+2*JAdjointRite*phi_unpert(npar.porder+1,npar.nel);
    BCqoiLeft=2*JForwardLeft*phia(1,1);%-2*JAdjointLeft*phi_unpert(1,1);
    if IO_opts.show_dqoi_from_bc
        fprintf('BCqoiRite value %g  \n',BCqoiRite);
        fprintf('BCqoiLeft value %g  \n',BCqoiLeft);
    end
    parts(4)= parts(4) + BCqoiLeft;
    parts(4)=parts(4) + BCqoiRite;
    d_qoi = d_qoi + BCqoiLeft + BCqoiRite;
end

function d_qoi = compute_perturbed_qoi_VEF(use_forward_flux,phia,phi_unpert,E)

global npar dat snq IO_opts

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
    d_qoi = d_qoi + Jac*dot(phia(:,iel), m*ones(2,1)*delta_qext) ;
    d_qoi = d_qoi - Jac*delta_siga*dot(phia(:,iel), m*phi_unpert(:,iel));
    %%Think I need to redo the below, once I figure out the
    for i=1:porder+1
        for j=1:porder+1
            k(i,j)= dot(wq.*dbdx(:,i) , Eloc.*dbdx(:,j)+dEdx*b(:,j));
        end
    end
    d_qoi = d_qoi - delta_isigtr/Jac*dot(phia(:,iel), k*phi_unpert(:,iel));
end

if IO_opts.show_dqoi_pre_bc
    fprintf('dqoi before bc %g (forward=%g) \n',d_qoi,use_forward_flux);
end

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
d_qoi = d_qoi + BCqoiLeft + BCqoiRite;

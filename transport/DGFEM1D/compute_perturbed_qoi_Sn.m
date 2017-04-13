function d_qoi = compute_perturbed_qoi_Sn(forward,phia,phi,E,psi,psia,is_sn)

global npar dat snq

sourcePert=dat.sourcePert;
sigaPert=dat.sigaPert;

% source data (volumetric) 
if forward
    qv  = dat.qv_forward;
    if is_sn
        qv = qv / snq.sw;
    end
else
    qv  = dat.qv_adjoint;
end

% type of solution data passed in
[n1,n2]=size(phia);
if n2==1
    phia=reshape(phia,npar.porder+1,npar.nel);
end
[n1,n2]=size(phi);
if n2==1
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
d_qoi=0;

% loop over elements
for iel=1:npar.nel
    my_zone=npar.iel2zon(iel);
    delta_sigs= dat.sigsPert(my_zone);
    delta_sigt= dat.sigtPert(my_zone);
    qext = qv(my_zone)*sourcePert;
    Jac   = npar.dx(iel)/2;
    % assemble
    d_qoi = d_qoi + Jac*dot(m*ones(2,1)*qext,phia(:,iel)) ;
    d_qoi = d_qoi + Jac*delta_sigs/snq.sw*dot(m*phi(:,iel),phia(:,iel));
    % Anglular integration of psi psia product
    for idir=1:snq.n_dir
        d_qoi = d_qoi - delta_sigt*Jac* snq.w(idir)* dot(m*psi(:,iel,idir), psia(:,iel,idir));
    end
end

% % loop over elements
% aux=zeros(9,1);
% aa=0;
% for iel=1:npar.nel
%     my_zone=npar.iel2zon(iel);
%     delta_sigs= dat.sigsPert(my_zone);
%     delta_sigt= dat.sigtPert(my_zone);
%     qext = qv(my_zone)*(1+sourcePert);
%     Jac   = npar.dx(iel)/2;
%     % assemble
% %     qoi = qoi + Jac*dot(m*ones(2,1)*qext,phia(:,iel)) ;
%     aux(9) = aux(9) + Jac*delta_sigs/snq.sw*dot(m*phi(:,iel),phia(:,iel));
%     % Anglular integration of psi psia product
%     for idir=1:snq.n_dir
%         aux(idir) = aux(idir) - delta_sigt*Jac* snq.w(idir)* dot(m*psi(:,iel,idir), psia(:,iel,idir));
%     end
%     for idir=1:snq.n_dir
%         aa = aa + Jac* snq.w(idir)* dot(m*ones(2,1), psia(:,iel,idir));
%     end
% end

% Add boundary terms

% Incident flux on right
iel=npar.nel;
dir_index_rite = 1:snq.n_dir/2; % the first half of the directions are <0
dir_index_left = snq.n_dir/2+1:snq.n_dir; % the second half of the directions are >0
incident_rite = dat.inc_forward(dir_index_rite) + dat.psiIncPert(dir_index_rite);
adjoint_flx_rite = shiftdim(psia(1,iel,dir_index_left),1);
% the reason for using the "other" direction for the adjoint flux is that
% psia(mu) = psi(-mu), or equivalently psi(-mu)=psi(mu) and remember we
% faked a "forward" solve for the adjoint

% Incident flux on left
iel=1;
% total perturbed incident
incident_left = dat.inc_forward(dir_index_left) + dat.psiIncPert(dir_index_left);  
adjoint_flx_left = shiftdim(psia(1,iel,dir_index_rite),1);
% 
%%%% d_qoi=d_qoi - dot(snq.w.*[incident_left incident_rite]',[adjoint_flx_left adjoint_flx_rite]);


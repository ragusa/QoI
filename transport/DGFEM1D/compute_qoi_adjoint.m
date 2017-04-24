function qoi = compute_qoi_adjoint(phia,is_sn,pert,phi,E,psi,psia)
% pert input: 
% 0 -> No Perturbation, just QoI
% 1 -> Perturbed QoI
% 2 -> Just sensitivity
% Need to add in input checks, should not be used if computing QoI from
% forward

global npar dat snq

if nargin < 5
   phi=0;
   psi = 0;
   psia=0;
end

% source data (volumetric)
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


% type of solution data passed in
[n1,n2]=size(phia);
if n2==1 % if DFEM solution store as vector, reshape it into a 2d array (local_dofs, elem)
    phia=reshape(phia,npar.porder+1,npar.nel);
    if phi~=0;
        phi=reshape(phi,npar.porder+1,npar.nel);
    end
end

% shortcuts
porder= npar.porder;
% quadrature
[xq,wq] = GLNodeWt(porder+1);
% store shapeset
[b,dbdx] =feshpln(xq,porder);
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
    if pert==0 || pert==1   
        qoi = qoi + Jac*qext*dot(ones(porder+1,1),m*phia(:,iel));
    end
    if pert==1 || pert==2
        if is_sn
            %Perturbation data
            delta_q=dat.sourcePert(my_zone);
            delta_q = delta_q / snq.sw;
            delta_sigs= dat.sigsPert(my_zone);
            delta_sigt= dat.sigtPert(my_zone);
            %Compute inner products and add to QoI
            qoi = qoi + Jac*delta_q*dot(ones(porder+1,1),m*phia(:,iel));
            qoi = qoi + Jac*delta_sigs/snq.sw*dot(phi(:,iel),m*phia(:,iel));
            for idir=1:snq.n_dir
                qoi = qoi - delta_sigt*Jac* snq.w(idir)* dot(psi(:,iel,idir), m*psia(:,iel,idir));
            end
        else
            %Perturbation data
            delta_siga= dat.sigaPert(my_zone);
            isigtr = 3*dat.cdif(my_zone); % inverse of sigma_tr 
            isigtrp= 3*1./(3*(dat.sigt(my_zone)+dat.sigtPert(my_zone)));
            delta_isigtr=isigtrp-isigtr;
            E0=E(1,iel);
            E1=E(2,iel);
            Eloc = (E1+E0)/2+xq*(E1-E0)/2;
            dEdx = (E1-E0)/2;
            delta_q=dat.sourcePert(my_zone);
            %Compute inner products and add to QoI
            qoi = qoi + Jac*delta_q*dot(ones(porder+1,1),m*phia(:,iel));
            qoi = qoi - delta_siga*Jac*dot(phi(:,iel), m*phia(:,iel));
            qoi = qoi + delta_isigtr*Jac*dot(phi(:,iel).*wq, (Eloc.*dbdx+dEdx.*b)*phia(:,iel));
        end
    end
end
%Boundary terms
%need to check this because our definition of "forward" in this context is confusing.
if is_sn
    %Add boundary terms if using the adjoitn to compute the QoI
    %dir_index_rite = 1:snq.n_dir/2; % the first half of the directions are <0
    %dir_index_left = snq.n_dir/2+1:snq.n_dir; % the second half of the directions are >0
    for idir=1:snq.n_dir
        %right boundary
        if(snq.mu(idir)>0)
            if pert==0 || pert==1 
                %qoi = qoi - snq.w(idir)/snq.sw*psi(porder+1,npar.nel,idir)*psia(porder+1,npar.nel,(snq.n_dir-idir+1));
                qoi = qoi - snq.w(idir)/snq.sw*psi(porder+1,npar.nel,idir)*dat.inc_adjoint(snq.n_dir-idir+1);
            end
        else
            if pert==0 || pert==1 
                %qoi = qoi + snq.w(idir)/snq.sw*psi(porder+1,npar.nel,idir)*psia(porder+1,npar.nel,(snq.n_dir-idir+1));
                qoi = qoi + snq.w(idir)/snq.sw*dat.inc_forward(idir)*psia(porder+1,npar.nel,(snq.n_dir-idir+1));
            end
            if pert==1 || pert==2 
                qoi = qoi + snq.w(idir)/snq.sw*dat.inc_forward_pert(idir)*psia(porder+1,npar.nel,(snq.n_dir-idir+1));
            end
        end
    end
    for idir=1:snq.n_dir
        %left boundary
        if(snq.mu(idir)>0)
            if pert==0 || pert==1 
                %qoi = qoi + snq.w(idir)/snq.sw*psi(1,1,idir)*psia(1,1,(snq.n_dir-idir+1));
                qoi = qoi + snq.w(idir)/snq.sw*dat.inc_forward(idir)*psia(1,1,(snq.n_dir-idir+1));
            end
            if pert==1 || pert==2 
                qoi = qoi + snq.w(idir)/snq.sw*dat.inc_forward_pert(idir)*psia(1,1,(snq.n_dir-idir+1));
            end
        else
            if pert==0 || pert==1 
                %qoi = qoi - snq.w(idir)/snq.sw*psi(1,1,idir)*psia(1,1,(snq.n_dir-idir+1));
                qoi = qoi - snq.w(idir)/snq.sw*psi(1,1,idir)*dat.inc_adjoint(snq.n_dir-idir+1);
            end
        end
    end
else
    x=1;
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

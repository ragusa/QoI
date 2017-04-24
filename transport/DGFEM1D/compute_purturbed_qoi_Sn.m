function qoi = compute_purturbed_qoi_Sn(forward,phia,phi,E,psi,psia)

global npar dat snq

sourcePert=dat.sourcePert;
sigaPert=dat.sigaPert;

% source data (volumetric) 
if forward
    qv  = dat.qv_forward;
else
    qv  = dat.qv_adjoint;
end

% type of solution data passed in
[n1,n2]=size(phia);
if n2==1
    phia=reshape(phia,npar.porder+1,npar.nel);
    phi=reshape(phi,npar.porder+1,npar.nel);
end

% shortcuts
porder= npar.porder;
% quadrature
[xq,wq] = GLNodeWt(porder+1);
% store shapeset
[b,dbdx] =feshpln(xq,porder);

qoi=0;
mass = [2 1;1 2]/6;
% loop over elements
for iel=1:npar.nel
    my_zone=npar.iel2zon(iel);
    delta_sigs= dat.sigsPert(my_zone);
    delta_sigt= dat.sigtPert(my_zone);
    qext = qv(my_zone)*(1+sourcePert);
    dx   = npar.dx(iel);
    % compute local matrices + load vector
    f= mass*(qext.*wq);
    % assemble
    qoi = qoi + dx*dot(f,phia(:,iel)) ;
    qoi = qoi + dx*delta_sigs/snq.sw*dot(mass*phi(:,iel),phia(:,iel));
    % Anglular integration of psi psia product
    for idir=1:snq.n_dir
        qoi = qoi - delta_sigt*dx* snq.w(idir)* dot(mass*psi(:,iel,idir), psia(:,iel,idir));
    end
end
%Add boundary terms
%Incident flux on left
lterm=0;
rterm=0;
for idir=1:snq.n_dir
        dx   = npar.dx(1);
        %m = mass * dx;
        %m=2*dx;
        incident=psi(2,1,idir)+dat.psiIncPert(idir);  %total perturbed incident
        adjointTerm = psia(1,1,idir);
        lterm=lterm+2*(snq.w(idir)*incident*adjointTerm);
end
for idir=1:snq.n_dir
        dx   = npar.dx(npar.nel);
        %m = mass * dx;
        %m=2*dx;
        incident=psi(2,npar.nel,idir)+dat.psiIncPert(idir);  %total perturbed incident
        adjointTerm = psia(2,npar.nel,idir);
        rterm=rterm+2*(snq.w(idir)*incident*adjointTerm);
end

lterm
rterm

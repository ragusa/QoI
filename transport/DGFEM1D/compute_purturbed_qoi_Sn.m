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

angleProduct=zeros(npar.porder+1,npar.nel);

% shortcuts
porder= npar.porder;
% quadrature
[xq,wq] = GLNodeWt(porder+1);
% store shapeset
[b,dbdx] =feshpln(xq,porder);

qoi=0;
ProductSum=0;
mass = [2 1;1 2]/6;
%mass = [6 1;1 6]/30;
%mass = [12 5;5 12]/60;
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
    psigs= dat.sigsPert(my_zone);
    psigt= dat.sigtPert(my_zone);
    qext = qv(my_zone)*(1+sourcePert);
    dx   = npar.dx(iel);
    % compute local matrices + load vector
    for i=1:porder+1
        f(i)= dot(qext.*wq, b(:,i));
        g(i)= (1/(4*3.14159))*psigs*dot(phi(:,iel).*wq, b(:,i));
        %Anglular integration of psi psia product
        %m = psigt*Jac*b;
    end
    m = mass * dx * psigt;
    % assemble
    qoi = qoi + Jac*dot(f,phia(:,iel)) ;
    qoi = qoi + Jac*dot(g,phia(:,iel)) ;
    for idir=1:snq.n_dir
        qoi=qoi-(snq.w(idir)*((psia(:,iel,idir))'*m*psi(:,iel,idir)));
    end
end
%Add boundary terms
%Incident flux on left
for idir=(snq.n_dir/2)+1:snq.n_dir
        dx   = npar.dx(1);
        m = mass * dx;
        incident=dat.inc_forward+dat.psiIncPert;  %total perturbed incident
        adjointTerm = m*psia(:,1,idir);
        qoi=qoi-(snq.w(idir)*incident(idir)*adjointTerm(1));
end
for idir=1:snq.n_dir/2
        dx   = npar.dx(npar.nel);
        m = mass * dx;
        incident=dat.inc_forward+dat.psiIncPert;  %total perturbed incident
        adjointTerm = m*psia(:,npar.nel,idir);
        qoi=qoi-(snq.w(idir)*incident(idir)*adjointTerm(2));
end

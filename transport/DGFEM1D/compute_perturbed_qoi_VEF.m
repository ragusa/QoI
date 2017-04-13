function qoi = compute_purturbed_qoi_VEF(forward,phia,phi,E,Ebd)

global npar dat

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
    psiga= dat.sigaPert(my_zone);
    qext = qv(my_zone)*(1+sourcePert);
    %get some data for sigt pert
    isigtr = 3*dat.cdif(my_zone); % inverse of sigma_tr 
    isigtrp= 3*1./(3*(dat.sigt(my_zone)+dat.sigtPert(my_zone)));
    deltaisigtr=isigtrp-isigtr;
    E0=E(1,iel);
    E1=E(2,iel);
    Eloc = (E1+E0)/2+xq*(E1-E0)/2;
    dEdx = (E1-E0)/2;
    % compute local matrices + load vector
    for i=1:porder+1
        f(i)= dot(qext.*wq, b(:,i));
        g(i)= psiga*dot(phi(:,iel).*wq, b(:,i));
        h(i)= deltaisigtr*dot(phi(:,iel).*wq, Eloc.*dbdx(:,i)+dEdx*b(:,i));
    end
    % assemble
    qoi = qoi + Jac*dot(f,phia(:,iel)) ;
    qoi = qoi - Jac*dot(g,phia(:,iel)) ;
    qoi = qoi + Jac*dot(h,phia(:,iel));
end


%boundary terms
%Left Boundary
    Jac = npar.dx(1)/2;
    E0=E(1,1);
    E1=E(2,1);
    Eloc = (E1+E0)/2+xq*(E1-E0)/2;
    dEdx = (E1-E0)/2;
    %my_zone=npar.iel2zon(1)
    isigtr = 3*dat.cdif(my_zone);
    for i=1:porder+1
        for j=1:porder+1
            k(i,j) = dot(wq.*dbdx(:,i) , Eloc.*dbdx(:,j)+dEdx*b(:,j));
        end
    end
    term1=isigtr*k/Jac;
    term2=term1*phi(:,1);
    testx=1;
%Right Boundry
    E0=E(1,end);
    E1=E(2,end);
    Eloc = (E1+E0)/2+xq*(E1-E0)/2;
    dEdx = (E1-E0)/2;

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


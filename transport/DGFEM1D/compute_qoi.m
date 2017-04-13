function qoi = compute_qoi(forward,phi,is_sn)

global npar dat snq 

% source data (volumetric) 
if forward
    qv  = dat.qv_forward;
    if is_sn
        qv = qv / snq.sw;
    end
else
    % this is q^\dagger(x,mu) = function given as input in load input. no need to tweak
    qv  = dat.qv_adjoint;
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


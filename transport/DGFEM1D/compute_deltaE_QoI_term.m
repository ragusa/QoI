function [d_qoi,BCdqoiLeft,BCdqoiRite,volValues] = compute_deltaE_QoI_term(phia,phi_unpert,phi_pert,E,E_pert,B,B_pert)

global npar dat snq IO_opts
d_qoi=0;
deltaE=E_pert-E;
deltaB=B_pert-B;
volValues=[];

% type of solution data passed in
[n1,n2]=size(phia);
if n2==1
    phia=reshape(phia,npar.porder+1,npar.nel);
end
[n1,n2]=size(phi_unpert);
if n2==1
    phi_unpert=reshape(phi_unpert,npar.porder+1,npar.nel);
end
[n1,n2]=size(phi_pert);
if n2==1
    phi_pert=reshape(phi_pert,npar.porder+1,npar.nel);
end

% shortcuts
porder= npar.porder;
% quadrature
[xq,wq] = GLNodeWt(porder+1);
% store shapeset
[b,dbdx] =feshpln(xq,porder);
% compute elementary mass matrix
md=zeros(porder+1,porder+1);
mdd=zeros(porder+1,porder+1);
for i=1:porder+1
    for j=1:porder+1
        md(i,j)= dot(wq.*dbdx(:,i), b(:,j));
        mdd(i,j)= dot(wq.*dbdx(:,i), dbdx(:,j));
    end
end
% loop over elements
for iel=1:npar.nel
    my_zone=npar.iel2zon(iel);
    isigtr = 3*dat.cdif(my_zone); % inverse of sigma_tr 
    DE0=deltaE(1,iel);
    DE1=deltaE(2,iel);
    Eloc = (DE1+DE0)/2;%+xq*(DE1-DE0)/2;
    dDEdx = (DE1-DE0)/2;
    Jac   = npar.dx(iel)/2;
    % assemble
    segment= - Jac*isigtr*dot(phi_unpert(:,iel), (Eloc*mdd+dDEdx*md)*phia(:,iel));
    volValues=[volValues segment];
    d_qoi = d_qoi +segment;
end

BCdqoiLeft=phia(1,1)*phi_pert(1,1)*deltaB(2);
BCdqoiRite=phia(npar.porder+1,npar.nel)*phi_pert(npar.porder+1,npar.nel)*deltaB(1);

return
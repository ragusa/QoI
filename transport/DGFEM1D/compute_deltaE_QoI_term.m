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
m=zeros(porder+1,porder+1);
md=zeros(porder+1,porder+1);
k=zeros(porder+1,porder+1);
stif=zeros(porder+1,porder+1);
for i=1:porder+1
    for j=1:porder+1
        m(i,j)= dot(wq.*b(:,i), b(:,j));
        md(i,j)= dot(wq.*dbdx(:,i), b(:,j));
        stif(i,j)= dot(wq.*dbdx(:,i), dbdx(:,j));
    end
end
% loop over elements
for iel=1:npar.nel
    my_zone=npar.iel2zon(iel);
    isigtr = 3*dat.cdif(my_zone); % inverse of sigma_tr 
    isigtrp= 3*1./(3*(dat.sigt(my_zone)+dat.sigtPert(my_zone)));
    DE0=deltaE(1,iel);
    DE1=deltaE(2,iel);
    Eloc = (DE1+DE0)/2+xq*(DE1-DE0)/2;
    dDEdx = (DE1-DE0)/2;
    Jac   = npar.dx(iel)/2;
    % assemble
    %segment= - isigtr/Jac*dot(phi_pert(:,iel), (stif*Eloc.*phia(:,iel)));
    for i=1:porder+1
        for j=1:porder+1
            k(i,j)= dot(wq.*dbdx(:,i) , Eloc.*dbdx(:,j)+dDEdx*b(:,j));
        end
    end
    %segment= - isigtr/Jac*dot(phi_pert(:,iel), (stif*Eloc.*phia(:,iel)+dDEdx*md*phia(:,iel)));
    %segment= - isigtr/Jac*dot(phi_pert(:,iel), (stif*(Eloc.*phia(:,iel))));
    segment= - isigtrp/Jac*dot(phia(:,iel), k*phi_unpert(:,iel));
    volValues=[volValues segment];
    d_qoi = d_qoi +segment;
end

BCdqoiLeft=-phia(1,1)*phi_unpert(1,1)*deltaB(2);
BCdqoiRite=-phia(npar.porder+1,npar.nel)*phi_unpert(npar.porder+1,npar.nel)*deltaB(1);

return
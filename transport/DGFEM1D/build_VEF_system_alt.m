function [A,rhs] = build_VEF_system_alt(forward,E,Ebd)

global npar dat

% source data (volumetric) 
if forward
    qv  = dat.qv_forward;
else
    qv  = dat.qv_adjoint;
end

dsa_method='SIP';
npar.Ckappa = 10;
npar.Ckappa_bd=2*npar.Ckappa;

% Ckappa = 1;
% Ckappa_bd = 1/400;
% shortcuts
porder= npar.porder;
gn    = npar.gn;
nel   = npar.nel;
% ideally, we would analyze the connectivity to determine nnz
nnz=(porder+3)*nel; %this is an upperbound, not exact
% n: linear system size
n=nel*(porder+1);
% allocate memory
A=spalloc(n,n,nnz);
K=A; % stiffness matrix
rhs=zeros(n,1);

% compute local matrices
% load Gauss Legendre quadrature (GLQ is exact on polynomials of degree up to 2n-1,
% while using only integrand evaluations (n-point quadrature).
% estimate the max poly order (it comes from the mass matrix  when coef. are
% piecewise constant and the material mesh matches the computational mesh
% poly_max=2*porder;
[xq,wq] = GLNodeWt(porder+1);
% initialize local matrices/vectors
m=zeros(porder+1,porder+1);
k=m;
f=zeros(porder+1,1);
% store shapeset
[b,dbdx] =feshpln(xq,porder);
% store shape set at edges
[be,dbedx] =feshpln([-1;1],porder);

% definition of the weak form:
% int_domain (grad u D grad b) - int_bd_domain (b D grad u n) ...
%      = int_domain (b rhs)

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
    isigtr = 3*dat.cdif(my_zone); % inverse of sigma_tr 
    siga = dat.siga(my_zone);
    qext = qv(my_zone);
    % get E values in the interval
    E0=E(1,iel);
    E1=E(2,iel);
    Eloc = (E1+E0)/2+xq*(E1-E0)/2;
    dEdx = (E1-E0)/2;
    % compute local matrices + load vector
    for i=1:porder+1
        for j=1:porder+1
            %Switched Eloc.*dbdx(:,j) to Eloc(j)*dbdx(:,j) - 8/10
            k(i,j) = dot(wq.*dbdx(:,i) , Eloc.*dbdx(:,j)+dEdx*b(:,j));
            m(i,j) = dot(wq.*b(:,i)    , b(:,j));
        end
        f(i)= dot(qext.*wq, b(:,i));
    end
    % assemble
    A(gn(iel,:),gn(iel,:)) = A(gn(iel,:),gn(iel,:)) + siga*m*Jac;
    K(gn(iel,:),gn(iel,:)) = K(gn(iel,:),gn(iel,:)) + isigtr*k/Jac;
    rhs(gn(iel,:)) = rhs(gn(iel,:)) + f*Jac;
end
% add VEF 
% Emat = sparse(diag(reshape(E,npar.ndofs,1)));
% A = A + K*Emat;
A=A+K;

% loop on interior edges
for ie=1:(npar.nel-1)
    % left (interior)/right (exterior) elements
    ieli = ie;
    iele = ie+1;
    % values for 2 cells 1= left of edge, 2= right of edge
    ce1=npar.iel2zon(ieli);
    ce2=npar.iel2zon(iele);
    % element lengths
    d1=npar.dx(ie)  ; Jac1=d1/2;
    d2=npar.dx(ie+1); Jac2=d2/2;
    % conductivity values
    cdif1=3*dat.cdif(ce1);%*E(2,ie)  ;
    cdif2=3*dat.cdif(ce2);%*E(1,ie+1);
    % compute local matrices + load vector
%     mee=(cdif2/d2)*dbedx(1,:)'*be(1,:);
%     mei=(cdif1/d1)*dbedx(end,:)'*be(1,:);
%     mie=(-cdif2/d2)*dbedx(1,:)'*be(end,:);
%     mii=(-cdif1/d1)*dbedx(end,:)'*be(end,:);
%     kap=max( Ckappa*(cdif1/d1+cdif2/d2) , 0.25 );
    % the 0.5 comes from the mean-value operator {{...}}
    mee=0.5*(cdif2/Jac2) *be(1,:)'     *dbedx(1,:)  ;
    mei=0.5*(cdif1/Jac1) *be(1,:)'     *dbedx(end,:);
    mie=0.5*(-cdif2/Jac2)*be(end,:)'   *dbedx(1,:)  ;
    mii=0.5*(-cdif1/Jac1)*be(end,:)'   *dbedx(end,:);
    % penalty coefficent:
    % assemble
    Ee=diag(E(:,ie+1));
    Ei=diag(E(:,ie  ));
    [ kap, sg ] = penalty_coef( dsa_method,cdif1*E(2,ie),cdif2*E(1,ie+1),d1,d2,npar.alpha );
    A(gn(ieli,:),gn(ieli,:)) = A(gn(ieli,:),gn(ieli,:)) + (mii + sg*mii')*Ei;
    A(gn(ieli,:),gn(iele,:)) = A(gn(ieli,:),gn(iele,:)) + (mie + sg*mei')*Ee;
    A(gn(iele,:),gn(ieli,:)) = A(gn(iele,:),gn(ieli,:)) + (mei + sg*mie')*Ei;
    A(gn(iele,:),gn(iele,:)) = A(gn(iele,:),gn(iele,:)) + (mee + sg*mee')*Ee;
%     % assemble
%     A(gn(ieli,:),gn(ieli,:)) = A(gn(ieli,:),gn(ieli,:)) + mii' + mii;
%     A(gn(ieli,:),gn(iele,:)) = A(gn(ieli,:),gn(iele,:)) + mie' + mei;
%     A(gn(iele,:),gn(ieli,:)) = A(gn(iele,:),gn(ieli,:)) + mei' + mie;
%     A(gn(iele,:),gn(iele,:)) = A(gn(iele,:),gn(iele,:)) + mee' + mee;
    % assemble
    A(gn(ieli,end),gn(ieli,end)) = A(gn(ieli,end),gn(ieli,end))+kap;
    A(gn(ieli,end),gn(iele,1))   = A(gn(ieli,end),gn(iele,1))  -kap;
    A(gn(iele,1),gn(ieli,end))   = A(gn(iele,1),gn(ieli,end))  -kap;
    A(gn(iele,1),gn(iele,1))     = A(gn(iele,1),gn(iele,1))    +kap;
end

% apply natural BC
% element lengths
d1=npar.dx(1);   Jac1=d1/2;
dn=npar.dx(end); Jacn=dn/2;
% conductivity values
cdif1=3*dat.cdif(npar.iel2zon(1))  *E(1,1);
cdifn=3*dat.cdif(npar.iel2zon(end))*E(2,end);
% compute load vector
% penalty coefficent:
switch dsa_method
    case 'MIP'
        kap1 = max( npar.Ckappa_bd*cdif1/d1 , npar.alpha);
        kapn = max( npar.Ckappa_bd*cdifn/dn , npar.alpha);
    case {'SIP','NIP'}
        kap1 = npar.Ckappa_bd*cdif1/d1;
        kapn = npar.Ckappa_bd*cdifn/dn;
    case 'M4S'
        kap1=npar.alpha;
        kapn=npar.alpha;
end
% kap1=max( Ckappa_bd*cdif1/d1 , 0.25);
% kapn=max( Ckappa_bd*cdifn/dn , 0.25);

% in this test code, we always have vacuum or incoming flux for transport,
% which means that we never need to apply a Neumann type bc in DSA.
% However, we still have a choice to implement the zero incoming current in
% dsa: Robin or Dirichlet.
% kap1 = max( npar.Ckappa_bd*cdif1/d1 , npar.alpha);
% kapn = max( npar.Ckappa_bd*cdifn/dn , npar.alpha);
% %%%%%%%%%%% [kap1 kapn]
% if npar.Ckappa_bd*cdif1/d1  > npar.alpha
%     kap1=npar.alpha*2;
% else
%     kap1=npar.alpha; % a factor 2 here would lead to the large hump for large MFPs
% end
% kapn=kap1;


% bc type:         0=neumann, 1=robin,               2=dirichlet
% C is defined as: kdu/dn=C   // u+k/hcv*du/dn =C // u=C
% u+k/hcv*du/dn = C
% u/4 + D/2*du/dn = Jinc
% u/2 +   D*du/dn = 2*Jinc. we let C=Jinc=0 in DSA
establish_bc_for_VEF( forward );

switch dat.bcVEF.left.type
    case 0 % Neumann, int_bd_domain (b D grad u n) is on the RHS
        rhs(1)=rhs(1)+dat.bcVEF.left.C;
    case 1 % Robin
        A(1,1)=A(1,1) + Ebd(1); % 1/2; %kap1; %1/2;
        rhs(1)=rhs(1) + 2*dat.bcVEF.left.C;
    case 2 % Dirichlet
        minus_one=1;
        A(gn(1,:),gn(1,:))=A(gn(1,:),gn(1,:))+kap1*be(minus_one,:)'*be(minus_one,:)...
            +0.5*cdif1*be(minus_one,:)'   *dbedx(minus_one,:)/Jac1... % 0.5 is to be equivalent to {{...}} on the bd
            +0.5*cdif1*dbedx(minus_one,:)'*be(minus_one,:)   /Jac1*sg;
%         A(gn(1,:),gn(1,:))=A(gn(1,:),gn(1,:))+kap1*be(minus_one,:)'*be(minus_one,:)...
%             +2*cdif1*be(minus_one,:)'   *dbedx(minus_one,:)/d1...
%             +2*cdif1*dbedx(minus_one,:)'*be(minus_one,:)   /d1;
        rhs(gn(1,:))=rhs(gn(1,:))+dat.bcVEF.left.C*...
            ( kap1*be(minus_one,:)' +2*cdif1*dbedx(minus_one,:)'/d1 );
end
switch dat.bcVEF.rite.type
    case 0 % Neumann, int_bd_domain (b D grad u n) is on the RHS
        rhs(n)=rhs(n)+dat.bcVEF.rite.C;
    case 1 % Robin
        A(n,n)=A(n,n) + Ebd(2); %+ 1/2; %kapn; %1/2;
        rhs(n)=rhs(n) + 2*dat.bcVEF.rite.C;
    case 2 % Dirichlet
        plus_one =2;
        A(gn(end,:),gn(end,:)) = A(gn(end,:),gn(end,:))...
            +kapn*be(plus_one,:)'*be(plus_one,:)...
            -0.5*cdifn*be(plus_one,:)'   *dbedx(plus_one,:)/Jacn...  % 0.5 is to be equivalent to {{...}} on the bd
            -0.5*cdifn*dbedx(plus_one,:)'*be(plus_one,:)   /Jacn*sg;
%         A(gn(end,:),gn(end,:))=A(gn(end,:),gn(end,:))...
%             +kapn*be(plus_one,:)'*be(plus_one,:)...
%             -2*cdifn*be(plus_one,:)'   *dbedx(plus_one,:)/dn...
%             -2*cdifn*dbedx(plus_one,:)'*be(plus_one,:)   /dn;
        rhs(gn(end,:))=rhs(gn(end,:))+dat.bcVEF.rite.C*...
            ( kapn*be(plus_one,:)' -2*cdifn*dbedx(plus_one,:)'/dn );
end

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


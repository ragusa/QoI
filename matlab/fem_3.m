function F=driver

% clear the console screen
clc; close all; 
% load the data structure with info pertaining to the physical problem
dat.material_type='constant_coefficients';
dat.response_type='constant';
dat.source_type='constant';

dat.diff=@diffu;
dat.siga=@siga;
dat.esrc=@esrc;
dat.resp=@resp;

dat.width=10;
% BCs
bc.left.type=2; %0=neumann, 1=robin, 2=dirichlet
bc.left.C=0; % (that data is C in: -Ddu/dn=C // u/4+D/2du/dn=C // u=C)
bc.rite.type=1;
bc.rite.C=0;
dat.bc=bc; clear bc;

% load the numerical parameters, npar, structure pertaining to numerics
% number of elements
npar.nel = 50;
% domain
npar.x = linspace(0,dat.width,npar.nel+1);
% polynomial degree
npar.porder=1;
% quadrature order
npar.q_pts=2;
% nbr of dofs per variable
npar.ndofs = npar.porder*npar.nel+1;
% connectivity
gn=zeros(npar.nel,npar.porder+1);
gn(1,:)=linspace(1,npar.porder+1,npar.porder+1);
for iel=2:npar.nel
    gn(iel,:)=[gn(iel-1,end) , gn(iel-1,2:end)+npar.porder ];
end
npar.gn=gn; clear gn;

% assemble the matrix and the rhs
[A,b,npar]=assemble_system(npar,dat);
% solve system
u=A\b;
% plot
figure(1)
plot(npar.x,u,'.-'); hold all
% verification is always good
% verif_diff_eq(dat)
[QoI_f] = QoIforward(dat.resp,npar,dat,u)
% u'*b/dat.esrc*dat.siga

% % us=u/dat.esrc*dat.siga;
% % bs=b/dat.esrc*dat.siga;
% % bs(1)  =dat.bc.left.C;
% % bs(end)=dat.bc.rite.C;
% % us=A\bs;

% <us|b_mod>=<us|A_mod.u>
%           =<us.(A_mod)'|u>
%           =<u|(A_mod)'.us>
%           =<u|r>
% assemble the matrix and the rhs
dat2=dat; dat2.esrc=respo; 
% we pass homog Neumann BC soo that the adjoint rhs 
% is untouched by the application of the BCs
dat2.bc.left.type=0; dat2.bc.rite.type=0;
dat2.bc.left.C=0;    dat2.bc.rite.C=0;
[As,bs,npar]=assemble_system(npar,dat2);
% we do not care about As
clear As As_nobc;
% us=As\bs;
us=A\bs;
plot(npar.x,us,'r-'); hold all


QoI_u_bstar = u'*bs
QoI_ustar_b = us'*b
% QoI_u_Astar_ustar = u'*(As*us)
% us'*(As*u)


% QoI_a = us'*b + us'*(As-A)*u
% QoI_a = us'*(A*u)+ us'*(As-A)*u
% u'*(A*us)

% us'*(A_nobc*u)
% u'*(As_nobc*us)

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % function u=solve_fem2(dat,npar)
% % 
% % % initial guess
% % % (not needed if a direct method is used to solve the linear system)
% % u=ones(npar.ndofs,1);
% % 
% % % assemble the matrix and the rhs
% % [A,b,npar]=assemble_system(npar,dat);
% % 
% % % solve
% % u=A\b;
% % 
% % return
% % end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [QoIf] = QoIforward(fct_ptr,npar,dat,u)

porder = npar.porder;
gn     = npar.gn;
[xq,wq] = GLNodeWt(npar.q_pts);
[b,dbdx] = feshpln(xq,porder);

QoIf = 0;
% loop over elements
for iel=1:npar.nel
    % element extremities
    x0=npar.x(iel);
    x1=npar.x(iel+1);
    % jacobian of the transformation to the ref. element
    Jac=(x1-x0)/2;
    % get x values in the interval
    x=(x1+x0)/2+xq*(x1-x0)/2;
    % compute function at xq
    value = fct_ptr(x,dat)*Jac;
    local_u  = b(:,:) * u(gn(iel,:));
    % assemble
    QoIf = QoIf + dot(value.*wq, local_u );
end

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,rhs,npar]=assemble_system(npar,dat)
% assemble the matrix, the rhs, apply BC
% shortcuts
porder= npar.porder;
gn    = npar.gn;
nel   = npar.nel;
% ideally, we would analyze the connectivity to determine nnz
nnz=(porder+1)*nel; %this is an upperbound, not exact
% n: linear system size
n=nel*porder+1;
% allocate memory
A=spalloc(n,n,nnz);
rhs=zeros(n,1);
if nargin<3
    field=zeros(n,1);
end

% compute local matrices
% load Gauss Legendre quadrature (GLQ is exact on polynomials of degree up to 2n-1,
% while using only integrand evaluations (n-point quadrature).
% estimate the max poly order (it comes from the mass matrix  when coef. are
% piecewise constant and the material mesh matches the computational mesh
[xq,wq] = GLNodeWt(npar.q_pts);
% initialize local matrices/vectors
m=zeros(porder+1,porder+1);
% k=m;
f=zeros(porder+1,1);
% store shapeset
[b,dbdx] =feshpln(xq,porder);
% % % compute local matrices + load vector
% % for i=1:porder+1
% %     for j=1:porder+1
% %         m(i,j)= dot(wq.*b(:,i)    , b(:,j));
% %         k(i,j)= dot(wq.*dbdx(:,i) , dbdx(:,j));
% %     end
% %     f(i)= dot(wq, b(:,i));
% % end
% % npar.m=m;
% % npar.k=k;
% % npar.f=f;

% definition of the weak form:
% int_domain (grad u D grad b) - int_bd_domain (b D grad u n) ...
%      + int_domain( XSa u b) = int_domain (b rhs)

% loop over elements
for iel=1:npar.nel
    % element extremities
    x0=npar.x(iel);
    x1=npar.x(iel+1);
    % jacobian of the transformation to the ref. element
    Jac=(x1-x0)/2;
    % get x values in the interval
    x=(x1+x0)/2+xq*(x1-x0)/2;
    % b and dbdx(:,:) are of size (nbr of xq values) x (porder+1),
    % field is of length porder+1, 2/dx is the 1d jacobian
%     local_field     = b(:,:) * field(gn(iel,:));
%     local_dfield_dx = dbdx(:,:) * field(gn(iel,:));
    d=dat.diff(x,field,dat)/Jac;
    s=dat.siga(x,field,dat)*Jac;
    q=dat.esrc(x,dat)*Jac;
    % compute local residual
    for i=1:porder+1
        for j=1:porder+1
            m(i,j) =  dot(s.*wq, b(:,i).*b(:,j)) + ...
                      dot(d.*wq, dbdx(:,i).*dbdx(:,j));
        end
        f(i) = dot(q.*wq, b(:,i));
    end
    % assemble
    A(gn(iel,:),gn(iel,:)) = A(gn(iel,:),gn(iel,:)) + m ;
    rhs(gn(iel,:)) = rhs(gn(iel,:)) + f;
end


% apply natural BC
Dirichlet_nodes=[];
Dirichlet_val=[];
switch dat.bc.left.type
    case 0 % Neumann, int_bd_domain (b D grad u n) is on the RHS
        dat.bc.left.C;
        rhs(1)=rhs(1)+dat.bc.left.C;
    case 1 % Robin
        A(1,1)=A(1,1)+1/2;
        rhs(1)=rhs(1)+2*dat.bc.left.C;
    case 2 % Dirichlet
        Dirichlet_nodes=[Dirichlet_nodes 1];
        Dirichlet_val=[Dirichlet_val dat.bc.left.C];
end
switch dat.bc.rite.type
    case 0 % Neumann, int_bd_domain (b D grad u n) is on the RHS
        rhs(n)=rhs(n)+dat.bc.rite.C;
    case 1 % Robin
        A(n,n)=A(n,n)+1/2;
        rhs(n)=rhs(n)+2*dat.bc.rite.C;
    case 2 % Dirichlet
        Dirichlet_nodes=[Dirichlet_nodes n];
        Dirichlet_val=[Dirichlet_val dat.bc.rite.C];
end
% apply Dirichlet BC
for i=1:length(Dirichlet_nodes);% loop on the number of constraints
    id=Dirichlet_nodes(i);      % extract the dof of a constraint
    bcval=Dirichlet_val(i);
    rhs=rhs-bcval*A(:,id);  % modify the rhs using constrained value
    A(id,:)=0; % set all the id-th row to zero
    A(:,id)=0; % set all the id-th column to zero (symmetrize A)
    A(id,id)=1;            % set the id-th diagonal to unity
    rhs(id)=bcval;         % put the constrained value in the rhs
end

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [shapefun,dhdx]=feshpln (xv,p)
% computes the finite element basis functions and their first derivatives
% for any order p
% here, we use Lagrangian FE functions

xd=linspace(-1,1,p+1);

shapefun=zeros(length(xv),p+1);
dhdx    =zeros(length(xv),p+1);

% shape function
for i=1:p+1
    num=1.;
    den=1.;
    for j=1:p+1
        if(j~=i)
            num=num.*(xv-xd(j));
            den=den.*(xd(i)-xd(j));
        end
    end
    shapefun(:,i)=num./den;
end

% derivative of the shape function
for i=1:p+1
    sum=0.;
    den=1.;
    for j=1:p+1
        if(j~=i)
            num=1;
            for k=1:p+1
                if((k~=i)&&(k~=j))
                    num=num.*(xv-xd(k));
                end
            end
            sum=sum+num;
            den=den.*(xd(i)-xd(j));
        end
    end
    dhdx(:,i)=sum./den;
end

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,w] = GL_quad(n)

if(n~=2), error('quadrature entered only for n=2'); end

x=[ -0.577350269189626; 0.577350269189626];
w=[1;1];

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,w] = GLNodeWt(n)
% GLNodeWt  Nodes and weights for Gauss-Legendre quadrature of arbitrary order
%           obtained by solving an eigenvalue problem
%
% Synopsis:  [x,w] = GLNodeWt(n)
%
% Input:     n = order of quadrature rule
%
% Output:    x = vector of nodes
%            w = vector of weights

%  Algorithm based on ideas from Golub and Welsch, and Gautschi.  For a
%  condensed presentation see H.R. Schwarz, "Numerical Analysis: A
%  Comprehensive Introduction," 1989, Wiley.  Original MATLAB
%  implementation by H.W. Wilson and L.H. Turcotte, "Advanced Mathematics
%  and Mechanics Applications Using MATLAB," 2nd ed., 1998, CRC Press

beta   = (1:n-1)./sqrt(4*(1:n-1).^2 - 1);
J      = diag(beta,-1) + diag(beta,1);    % eig(J) needs J in full storage
[V,D]  = eig(J);
[x,ix] = sort(diag(D));  %  nodes are eigenvalues, which are on diagonal of D
w      = 2*V(1,ix)'.^2;  %  V(1,ix)' is column vector of first row of sorted V

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function verif_diff_eq(dat)

sa=dat.siga; cd=dat.diff; src=dat.esrc; L=dat.width;

k=1/sqrt(cd/sa);
% particular solution
part=src/sa;
% general form of the solution:
%  phi =A sh(kx) + B ch(kx) + src/abso
%  dphi/dx= k (A ch(kx) + B sh(kx) )
switch dat.bc.left.type
    case 0 % Neumann
        mat(1,1:2) =k*[ cosh(0), sinh(0)];
        b(1) = -dat.bc.left.C / cd;
    case 1 % Robin
        mat(1,1:2) =[ sinh(0), cosh(0)]/4 -cd/2*k*[ cosh(0), sinh(0)];
        b(1) = dat.bc.left.C -part/4;
    case 2 % Dirichlet
        mat(1,1:2) =[ sinh(0), cosh(0)];
        b(1) = dat.bc.left.C -part;
end
switch dat.bc.rite.type
    case 0 % Neumann
        mat(2,1:2) =k*[ cosh(k*L), sinh(k*L)];
        b(2) = dat.bc.rite.C / cd;
    case 1 % Robin
        mat(2,1:2) =[ sinh(k*L), cosh(k*L)]/4 +cd/2*k*[ cosh(k*L), sinh(k*L)];
        b(2) = dat.bc.rite.C -part/4;
    case 2 % Dirichlet
        mat(2,1:2) =[ sinh(k*L), cosh(k*L)];
        b(2) = dat.bc.rite.C -part;
end
% get coefficient for the analytical solution
a=mat\b';
x=linspace(0,L);
y=a(1)*sinh(k*x)+a(2)*cosh(k*x)+part;
plot(x,y,'r-');hold all

% ff= @(x) dat.siga*(a(1)*sinh(k*x)+a(2)*cosh(k*x)+part);
% num_qoi=quad(ff,0,dat.width)

return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=diffu(x,field,dat);
if(strcmp(dat.material_type,'constant_coefficients'))
    y=3/10;
else
    error('non constant diffu not entered yet');
end
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=siga(x,field,dat);
if(strcmp(dat.material_type,'constant_coefficients'))
    y=3/10;
else
    error('non constant siga not entered yet');
end
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=nothing(x,field);
y=0;
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=esrc(x,dat);
if(strcmp(dat.source_type,'constant'))
    y=1;
else
    error('non constant esrc not entered yet');
end

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=resp(x,dat);
% make the problem-data a global variable
if(strcmp(dat.response_type,'constant'))
    y=1;
else
    error('non constant response not entered yet');
end

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

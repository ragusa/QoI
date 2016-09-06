function F=driver
% clear the console screen
clc; close all; clf
qoi_real = 1;
% Flat Solutions using the 3 BC types
[qoif, qoia] = QOI(1000, 1,1,1,1,1,0,0,0,0);
F_Err = abs(qoif - qoi_real)
FvsA = abs(qoif - qoia)

[qoif, qoia] = QOI(1000, 1,1,1,1,1,1,1/4,1,1/4);
F_Err = abs(qoif - qoi_real)
FvsA = abs(qoif - qoia)

[qoif, qoia] = QOI(1000, 1,1,1,1,1,2,1,2,1);
F_Err = abs(qoif - qoi_real)
FvsA = abs(qoif - qoia)
%/ % load the data structure with info pertaining to the physical problem
%/ dat.diff=1;
%/ dat.siga=3/10;
%/ dat.esrc=3;
%/ dat.width=10;
%/ bc.left.type=1; %0=neumann, 1=robin, 2=dirichlet
%/ bc.left.C=0; % (that data is C in: -Ddu/dn=C // u/4+D/2du/dn=C // u=C)
%/ bc.rite.type=1;
%/ bc.rite.C=0;
%/ dat.bc=bc; clear bc;
%/ 
%/ % load the numerical parameters, npar, structure pertaining to numerics
%/ % number of elements
%/ npar.nel = 5000;
%/ % domain
%/ npar.x = linspace(0,dat.width,npar.nel+1);
%/ % polynomial degree
%/ npar.porder=1;
%/ % nbr of dofs per variable
%/ npar.ndofs = npar.porder*npar.nel+1;
%/ % connectivity
%/ gn=zeros(npar.nel,npar.porder+1);
%/ gn(1,:)=linspace(1,npar.porder+1,npar.porder+1);
%/ for iel=2:npar.nel
%/     gn(iel,:)=[gn(iel-1,end) , gn(iel-1,2:end)+npar.porder ];
%/ end
%/ npar.gn=gn; clear gn;
%/ 
%/ % assemble the matrix and the rhs
%/ [A,b,npar,A_nobc]=assemble_system(npar,dat);
%/ % solve system
%/ u=A\b;
%/ % plot
%/ figure(1)
%/ plot(npar.x,u,'.-'); hold all
%/ % verification is always good
%/ % verif_diff_eq(dat)
%/ respo=1;
%/ [QoIf] = QoIforward(respo,npar,u)
%/ % u'*b/dat.esrc*dat.siga
%/ 
%/ % % us=u/dat.esrc*dat.siga;
%/ % % bs=b/dat.esrc*dat.siga;
%/ % % bs(1)  =dat.bc.left.C;
%/ % % bs(end)=dat.bc.rite.C;
%/ % % us=A\bs;
%/ 
%/ % assemble the matrix and the rhs
%/ dat2=dat; dat2.esrc=respo; 
%/ dat2.bc.left.type=1; dat2.bc.rite.type=0;
%/ dat2.bc.left.C=0;    dat2.bc.rite.C=0;
%/ [As,bs,npar,As_nobc]=assemble_system(npar,dat2);
%/ us=As\bs;
%/ plot(npar.x,us,'r-'); hold all
%/ 
%/ QoI_u_bstar       = u'*bs
%/ QoI_u_Astar_ustar = u'*(As*us)
%/ us'*(As*u)
%/ 
%/ 
%/ QoIa = us'*b + us'*(As-A)*u
%/ QoIa = us'*(A*u)+ us'*(As-A)*u
%/ % u'*(A*us)
%/ 
%/ % us'*(A_nobc*u)
%/ % u'*(As_nobc*us)
%/ 
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % function [qoif, qoia] = QOI(D, SIGMA, SRC, Width, RES, BC_L_T, BC_L_C, BC_R_T, BC_R_C)
% % return qoia and qoif

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [qoif, qoia] = QOI(N, D, SIGMA, SRC, WIDTH, RES, BC_L_T, BC_L_C, BC_R_T, BC_R_C)
	dat.diff=D;
	dat.siga=SIGMA;
	dat.esrc=SRC;
	dat.width=WIDTH;
	bc.left.type=BC_L_T; %0=neumann, 1=robin, 2=dirichlet
	bc.left.C=BC_L_C; % (that data is C in: -Ddu/dn=C // u/4+D/2du/dn=C // u=C)
	bc.rite.type=BC_R_T;
	bc.rite.C=BC_R_C;
	dat.bc=bc; clear bc;

	% load the numerical parameters, npar, structure pertaining to numerics
	% number of elements
	npar.nel = N;
	% domain
	npar.x = linspace(0,dat.width,npar.nel+1);
	% polynomial degree
	npar.porder=1;
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
	[A,b,npar,A_nobc]=assemble_system(npar,dat);
	% solve system
	u=A\b;
	% verification is always good
	% verif_diff_eq(dat)
	respo=RES;
	[qoif] = QoIforward(respo,npar,u);

	% assemble the matrix and the rhs
	dat2=dat; dat2.esrc=respo; 
	dat2.bc.left.type=0; dat2.bc.rite.type=0;
	dat2.bc.left.C=0;    dat2.bc.rite.C=0;
	[As,bs,npar,As_nobc]=assemble_system(npar,dat2);
	us=As\bs;

	qoia = us'*b + us'*(As-A)*u;
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

function [QoIf] = QoIforward(value,npar,u)

QoIf = 0;
f=npar.f;
% loop over elements
for iel=1:npar.nel
    % element extremities
    x0=npar.x(iel);
    x1=npar.x(iel+1);
    % jacobian of the transformation to the ref. element
    Jac=(x1-x0)/2;
    % assemble
    QoIf = QoIf + value*f'*Jac*u(npar.gn(iel,:));
end

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,rhs,npar,A_nobc]=assemble_system(npar,dat)
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

% compute local matrices
% load Gauss Legendre quadrature (GLQ is exact on polynomials of degree up to 2n-1,
% while using only integrand evaluations (n-point quadrature).
% estimate the max poly order (it comes from the mass matrix  when coef. are
% piecewise constant and the material mesh matches the computational mesh
q_pts=2;
[xq,wq] = GL_quad(q_pts);
% initialize local matrices/vectors
m=zeros(porder+1,porder+1);
k=m;
f=zeros(porder+1,1);
% store shapeset
[b,dbdx] =feshpln(xq,porder);
% compute local matrices + load vector
for i=1:porder+1
    for j=1:porder+1
        m(i,j)= dot(wq.*b(:,i)    , b(:,j));
        k(i,j)= dot(wq.*dbdx(:,i) , dbdx(:,j));
    end
    f(i)= dot(wq, b(:,i));
end
npar.m=m;
npar.k=k;
npar.f=f;

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
    % assemble
    A(gn(iel,:),gn(iel,:)) = A(gn(iel,:),gn(iel,:)) + ...
        dat.siga*m*Jac + dat.diff*k/Jac;
    rhs(gn(iel,:)) = rhs(gn(iel,:)) + dat.esrc*f*Jac;
end
A_nobc=A;

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
% here, we use Lagrangian FE functions

if(p~=1), error('only linear FEM'); end
if(length(xv)~=2), error('only 2-pt quadrature'); end

shapefun=zeros(length(xv),p+1);
dhdx    =zeros(length(xv),p+1);

shapefun=[
         0.788675134594813         0.211324865405187;
         0.211324865405187         0.788675134594813];
dhdx=[                -0.5                       0.5;
                      -0.5                       0.5];


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

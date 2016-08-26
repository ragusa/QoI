function hc_ss_x_CG
% Solves the heat conduction equation in 1-D x-geometry using CFEM
% without T gap.
% An arbitrary number of material zones can be used but the analytical
% solution assumes 3 zones are used. The conductivities and the volumetric
% sources can be spatially dependent.

% clear the console screen
clc; close all;

% load the data and numerical parameters
pert_k    = 1e-1 *1;
pert_s    = 1e-1 *0;
pert_bc_L = 1e-1 *0;
pert_bc_R = 1e-1 *0;
[dat,npar] = load_simulation_data(pert_k,pert_s,pert_bc_L,pert_bc_R);

% type of data used to run the simulation
npar.pert_status = 'perturbed';

% assemble the matrix and the rhs
npar.adjoint=false;
npar.pert_status = 'unperturbed';
[Au,qu]=assemble_system(npar,dat);
npar.pert_status = 'perturbed';
[Ap,qp]=assemble_system(npar,dat);
% solve forward system
Tu=Au\qu;
Tp=Ap\qp;

% assemble the matrix and the rhs
npar.adjoint=true; 
npar.pert_status = 'unperturbed';
[Aa,r]=assemble_system(npar,dat);
% solve forward system
phi=Aa\r;

dT = Tp-Tu;
dq = qp-qu;
dA = Ap-Au;

J_for_unpert = Tu'*r
J_for_pert   = Tp'*r
dJ_for = dT'*r

J_adj_unpert = phi'*qu
dJ_adj = phi'*(dq - dA*Tu)
dJ_adj_exact = phi'*(dq-dA*Tu-dA*dT)

% (Tp,Aa.phi) = (Tp,r) = (Ap^{-1}qp,r) = (qp,Ap^{-T}r)
% plot solution
figure(1)
subplot(1,2,1)
plot(npar.xf,Tu,'.-',npar.xf,Tp,'+-');
legend(['Tmp';'Adj']);
subplot(1,2,2)
plot(npar.xf,phi,'r-');
legend(['Tmp';'Adj']);

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,rhs]=assemble_system(npar,dat)

% assemble the matrix, the rhs, apply BC

% shortcuts
porder= npar.porder;
gn    = npar.gn;
nel   = npar.nel;
% ideally, we would analyze the connectivity to determine nnz
nnz=(porder+3)*nel; %this is an upperbound, not exact
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
poly_max=2*porder;
[xq,wq] = GLNodeWt(porder+1);
% initialize local matrices/vectors
m=zeros(porder+1,porder+1);
k=m;
f=zeros(porder+1,1);
% store shapeset
[b,dbdx] =feshpln(xq,porder);

% definition of the weak form:
% int_domain (grad u D grad b) - int_bd_domain (b D grad u n) ...
%      = int_domain (b rhs)

% switch between forward and adjoint calculations
if npar.adjoint
    src = dat.asrc;
    bc  = dat.bc_adj;
    kond = dat.k;
else
    bc  = dat.bc_for;
    % what kind of material properties to use
    switch npar.pert_status
        case 'unperturbed'
            kond = dat.k;
            src = dat.fsrc;
        case 'perturbed'
            for i=1:length(dat.k)
                kond{i} = @(x) dat.k{i}(x)    + dat.dk{i}(x);
                src{i}  = @(x) dat.fsrc{i}(x) + dat.dfsrc{i}(x);
            end
            bc.left.C = bc.left.C + bc.left.dC;
            bc.rite.C = bc.rite.C + bc.rite.dC;
        case 'delta'
            kond = dat.dk;
            src = dat.dfsrc;
            bc.left.C =  bc.left.dC;
            bc.rite.C =  bc.rite.dC;
    end
end

% loop over elements
for iel=1:npar.nel
    % element extremities
    x0=npar.x(iel);
    x1=npar.x(iel+1);
    % jacobian of the transformation to the ref. element
    Jac=(x1-x0)/2;
    % get x values in the interval
    x=(x1+x0)/2+xq*(x1-x0)/2;
    my_zone=npar.iel2zon(iel);
    d=kond{my_zone}(x);
    q=src{my_zone}(x);
    % compute local matrices + load vector
    for i=1:porder+1
        for j=1:porder+1
            k(i,j)= dot(d.*wq.*dbdx(:,i) , dbdx(:,j));
        end
        f(i)= dot(q.*wq, b(:,i));
    end
    % assemble
    A(gn(iel,:),gn(iel,:)) = A(gn(iel,:),gn(iel,:)) + k/Jac;
    rhs(gn(iel,:)) = rhs(gn(iel,:)) + f*Jac;
end

% apply natural BC
Dirichlet_nodes=[];
Dirichlet_val=[];
switch bc.left.type
    case 0 % Neumann, int_bd_domain (b D grad u n) is on the RHS
        rhs(1)=rhs(1)+bc.left.C;
    case 1 % Robin
        A(1,1)=A(1,1)+dat.hcv;
        rhs(1)=rhs(1)+dat.hcv*bc.left.C;
    case 2 % Dirichlet
        Dirichlet_nodes=[Dirichlet_nodes 1];
        Dirichlet_val=[Dirichlet_val bc.left.C];
end
switch bc.rite.type
    case 0 % Neumann, int_bd_domain (b D grad u n) is on the RHS
        rhs(n)=rhs(n)+bc.rite.C;
    case 1 % Robin
        A(n,n)=A(n,n)+dat.hcv;
        rhs(n)=rhs(n)+dat.hcv*bc.rite.C;
    case 2 % Dirichlet
        Dirichlet_nodes=[Dirichlet_nodes n];
        Dirichlet_val=[Dirichlet_val bc.rite.C];
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
%  condensed presentation see H.R. Schwarz, Numerical Analysis: A
%  Comprehensive Introduction, 1989, Wiley.  Original MATLAB
%  implementation by H.W. Wilson and L.H. Turcotte, Advanced Mathematics
%  and Mechanics Applications Using MATLAB, 2nd ed., 1998, CRC Press

beta   = (1:n-1)./sqrt(4*(1:n-1).^2 - 1);
J      = diag(beta,-1) + diag(beta,1);    % eig(J) needs J in full storage
[V,D]  = eig(J);
[x,ix] = sort(diag(D));  %  nodes are eigenvalues, which are on diagonal of D
w      = 2*V(1,ix)'.^2;  %  V(1,ix)' is column vector of first row of sorted V

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dat,npar]=load_simulation_data(pert_k,pert_s,pert_bc_L,pert_bc_R)

triga = false;
% load the data structure with info pertaining to the physical problem
if triga
    dat.k{1}=@k_zr; % W/m-K
    dat.k{2}=@k_fuel; % W/m-K
    dat.k{3}=@k_fuel; % W/m-K
    dat.k{4}=@k_clad; % W/m-K
    % forward calculation source
    dat.fsrc{1}=@zero_function; % W/m3
    dat.fsrc{2}=@esrc; % W/m3
    dat.fsrc{3}=@esrc; % W/m3
    dat.fsrc{4}=@zero_function; % W/m3
    % adjoint calculation source
    dat.asrc{1}=@response_function;
    dat.asrc{2}=@zero_function;
    dat.asrc{3}=@zero_function;
    dat.asrc{4}=@zero_function; % W/m3
else
    dat.k{1}=@(x) 2150/200*(1+0*x);
    dat.fsrc{1}=@(x) 10000*(1+0*x);
    dat.asrc{1}=@(x) 2*(1+0*x);
end
% create perturbations
for i=1:length(dat.k)
    dat.dk{i}    = @(x) pert_k*dat.k{i}(x);
    dat.dfsrc{i} = @(x) pert_s*dat.fsrc{i}(x);
end

% dimensions
if triga
    % dat.hcv = 1612.414; % W/m^2-C
    dat.hcv = 100; % W/m^2-C
    % dat.width = [0.003175 0.01 0.0174115 0.0179195]; % m
    dat.width = [1 3 4 7]; % m
    % mesh resolution per region
    nel_zone = [ 5 5 5 2];
else
    dat.hcv = 16;
    dat.width = 0.5;
    nel_zone = [ 10 ];
end
% forward bc
bc.left.type=0; %0=neumann, 1=robin, 2=dirichlet
bc.left.C=0; % (that data is C in: kdu/dn=C // u+k/hcv*du/dn =C // u=C)
if triga
    bc.rite.type=2;
    bc.rite.C=15;
else
    bc.rite.type=1;
    bc.rite.C=300;
end
dat.bc_for=bc;
% create perturbations
dat.bc_for.left.dC = pert_bc_L*dat.bc_for.left.C;
dat.bc_for.rite.dC = pert_bc_R*dat.bc_for.rite.C;
% adjoint bc
bc.left.C=0;
bc.rite.C=0;
dat.bc_adj=bc;
clear bc;

% load the numerical parameters, npar, structure pertaining to numerics
% number of elements
if length(dat.width)~=length(nel_zone)
    error('not the same number of zones in dat.width and nel_zone');
end
if length(dat.k)~=length(nel_zone)
    error('not the same number of zones in dat.k and nel_zone');
end
if length(dat.fsrc)~=length(nel_zone)
    error('not the same number of zones in dat.fsrc and nel_zone');
end
if length(dat.asrc)~=length(nel_zone)
    error('not the same number of zones in dat.asrc and nel_zone');
end
npar.nel = sum(nel_zone);

% domain
tmp_width = [ 0 dat.width];
x=[];
iel2zon=[];
for i_zone=1:length(nel_zone)
    x_zone = linspace(tmp_width(i_zone),tmp_width(i_zone+1),nel_zone(i_zone)+1);
    if isempty(x)
        x = x_zone;
        iel2zon=i_zone*ones(nel_zone(i_zone),1);
    else
        x = [x x_zone(2:end)];
        iel2zon =[ iel2zon; i_zone*ones(nel_zone(i_zone),1)];
    end
end
npar.x=x;
npar.iel2zon=iel2zon;

% polynomial degree
npar.porder=2;
% nbr of dofs per variable
npar.ndofs = npar.porder*npar.nel+1;

% connectivity
gn=zeros(npar.nel,npar.porder+1);
gn(1,:)=linspace(1,npar.porder+1,npar.porder+1);
for iel=2:npar.nel
    gn(iel,:)=[gn(iel-1,end) , gn(iel-1,2:end)+npar.porder ];
end
npar.gn=gn; clear gn;

% create x values for plotting FEM solution
npar.xf=zeros(npar.nel*npar.porder+1,1);
for iel=1:npar.nel
    ibeg = (iel-1)*npar.porder+1;
    iend = (iel  )*npar.porder+1;
    npar.xf(ibeg:iend)=linspace(npar.x(iel),npar.x(iel+1),npar.porder+1) ;
end

return
end
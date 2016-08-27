function hc_ss_x_CG_nonlinear
% Solves the heat conduction equation in 1-D x-geometry using CFEM
% without T gap.
% problem definition
% -d/dx (k(T) dT/dx) = q  with k(T)=a/(b+T)+c
% dTdx|_0=0 and -k.dTdx|_L=h(T(L)-T_inf)

% clear the console screen
clc; close all;

% load the data and numerical parameters
pert.a    = 1e-1 *0.;
pert.b    = 1e-1 *0.;
pert.c    = 0;
pert.q    = 1e-1 *0;
pert.bc_L = 1e-1 *0;
pert.bc_R = 1e-1 *1;
[dat,npar] = load_simulation_data(pert);

% solve unperturbed forward problem
npar.adjoint=false;
npar.pert_status = 'unperturbed';
npar.do_analytical=false;
[Tu,Au,qu]=solve_system(npar,dat);

% solve perturbed forward problem
npar.adjoint=false;
npar.pert_status = 'perturbed';
[Tp,Ap,qp]=solve_system(npar,dat,Tu);
%
plot(npar.xf,Tu,npar.xf,Tp);

% solve unperturbed adjoint problem
npar.adjoint=true; 
npar.pert_status = 'unperturbed';
[phiu,Aau,ru]=solve_system(npar,dat,Tu,Tu);

% solve perturbed adjoint problem
npar.adjoint=true; 
npar.pert_status = 'perturbed';
[phip,Aap,rp]=solve_system(npar,dat,Tu,Tu);

% QoI summary
J_for_unpert = Tu'*ru;
J_for_pert   = Tp'*ru;
J_adj_unpert = phiu'*qu;
fprintf('-------------QoI------------\n')
fprintf('J forward unperturbed\t%14.8g \n',J_for_unpert);
fprintf('J adjoint unperturbed\t%14.8g \n',J_adj_unpert);
% sensitivity
dT = Tp-Tu;
dq = qp-qu; % this is how the bc mods to the rhs get eliminated. need to check. especially for dirichlet !!!
% dA = Ap-Au; caveat: application of bc disappear whend oing this !!!!
dJ_for = dT'*ru;
fprintf('J forward perturbed \t%14.8g \n',J_for_pert);
fprintf('dJ forward \t%14.8g \t%14.8g\n',(J_for_pert-J_for_unpert)/pert.bc_R,dJ_for/pert.bc_R);
% assembly system to compute sitffness matrix with conductivity=dkdp
npar.adjoint=false; 
npar.pert_status = 'delta_p';
[dAdp,~]=assemble_system(npar,dat,Tu);
dJ_adj = phip'*(dq - dAdp*Tu);
fprintf('dJ adjoint \t%14.8g \n',dJ_adj/pert.bc_R);
% dJ_adj = phiu'*(dq - dAdp*Tu);
% fprintf('dJ adjoint \t%14.8g \n',dJ_adj);

phip'*dAdp*Tu



% plot solution
figure(1)
subplot(1,2,1);
plot(npar.xf,Tu,'.-',npar.xf,Tp,'+-')
legend(['Tmp';'Adj']);
subplot(1,2,2);
plot(npar.xf,phiu,'.-',npar.xf,phip,'+-')
legend(['Tmp';'Adj']);

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [T,A,q]=solve_system(npar,dat,T_init,T_for)

%check
if nargin==4 && npar.adjoint==false
    error('Only adjoint calculations require 4 input arguments in solve_system');
end
% initial guess
if nargin==2
    T=100*ones(npar.ndofs,1);
else
    T=T_init;
end

% nonlinear solve
for iter=1:npar.max_nl_iter
    % assemble system
    if npar.adjoint
        T_for_assembly=T_for;
    else
        T_for_assembly=T;
    end
    [A,q]=assemble_system(npar,dat,T_for_assembly);
    % save old values
    T_old=T;
    % solve linear system
    T=A\q;
    % compute error
    err=norm(T_old-T);
    % check convergence
    fprintf('Picard=%g, error=%d\n',iter,err);
    if err<npar.nl_tol
        fprintf('Converged\n\n');
        break; % done
    else
        if iter==npar.max_nl_iter
            warning('Exiting Picard without converging to tolerance.  err=%d, tol=%d',err,npar.nl_tol);
        end
    end
end

% plot numerical solution and analytical one
% analytical solution
if npar.do_analytical
    a=dat.conductivity_constants.a;
    b=dat.conductivity_constants.b;
    L=sum(dat.width);
    T_inf=dat.bc_for.rite.C;
    h=dat.hcv;
    q=dat.src_strength;
    T_L=q*L/h+T_inf;
    % expression
    fh_T = @(x) -b+(b+T_L)*exp(q/2/a*(L^2-x.^2));
    figure(1);
    x=linspace(0,L,10000);
    plot(x,fh_T(x),'.-',npar.xf,T,'r+-'); axis tight; title('Temperature');legend(['analytical';'numerical ']);
end


return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,rhs]=assemble_system(npar,dat,T)

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
                kond{i} = @(T) dat.k{i}(T)    + dat.dka{i}(T)+ dat.dkb{i}(T)+ dat.dkc{i}(T);
                src{i}  = @(x) dat.fsrc{i}(x) + dat.dfsrc{i}(x);
            end
            bc.left.C = bc.left.C + bc.left.dC;
            bc.rite.C = bc.rite.C + bc.rite.dC;
        case 'delta_p'
            for i=1:length(dat.k)
                kond{i} = @(T) dat.dka{i}(T)+ dat.dkb{i}(T)+ dat.dkc{i}(T);
            end
            src = dat.dfsrc;
            bc.left.type =  0; % Neumann
            bc.rite.type =  0;
            bc.left.C =  0;
            bc.rite.C =  0;
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
    % get the mat zone ID
    my_zone=npar.iel2zon(iel);
    % evaluate external source at qp
    q=src{my_zone}(x);
    % evaluate previous temperature at qp
    T_qp = b(:,:) * T(gn(iel,:));
    % evaluate
    d=kond{my_zone}(T_qp);
    % compute local matrices + load vector
    for i=1:porder+1
        for j=1:porder+1
            k(i,j) = dot(d.*wq.*dbdx(:,i) , dbdx(:,j)) /Jac;
        end
        f(i) = dot(q.*wq, b(:,i));
    end
    if npar.adjoint && strcmp(npar.pert_status,'perturbed')
        dTdx = dbdx(:,:) * T(gn(iel,:));
        dkdT=dat.dkdT{my_zone}(T_qp);
        for i=1:porder+1
            for j=1:porder+1
                k(i,j) = k(i,j) + dot(dkdT.*dTdx.*wq.*b(:,i) , dbdx(:,j))/Jac;
            end
        end
    end
    % assemble
    A(gn(iel,:),gn(iel,:)) = A(gn(iel,:),gn(iel,:)) + k;
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

function [dat,npar]=load_simulation_data(pert)

% data 
T_inf=300; h=16; L=0.5;

a=2150; b=200;c=1.05*0; d=1;

cond = @(T) a./(d*T+b)+c;

dat.conductivity_constants.a=a;
dat.conductivity_constants.b=b;
% dat.conductivity_constants.c=c;

q=10000;  sq = @(x) q*(1+0*x);
dat.src_strength=q;

res_funct = @(x) 1/L*(1+0*x);

% load the data structure with info pertaining to the physical problem
dat.k{1}=cond; % W/m-K
% forward calculation source
dat.fsrc{1}=sq; 
% adjoint calculation source
dat.asrc{1}=res_funct; 
% create perturbations. at the input level, pert must contain fractional
% changes. pert.a=-0.1 means that the value of a will be modified by -10%
% k = a/(b+T) + c;
% dk/da = 1/(b+T);
% dk/db = -a/(b+T)^2;
% dk/dc = 1;
for i=1:length(dat.k)
    dat.dka{i}   = @(T) 1./(d*T+b)*pert.a*a;
    dat.dkb{i}   = @(T) -a./((d*T+b).^2)*pert.b*b;
    dat.dkc{i}   = @(T) pert.c*c*(1+0*T);
    dat.dkdT{i}  = @(T) -a./((d*T+b).^2)*d;
    dat.dfsrc{i} = @(x) pert.q*dat.fsrc{i}(x);
end

% dimensions
dat.hcv = h; 
dat.width = L; 
% mesh resolution per region
nel_zone = [ 10 ];

% forward bc
bc.left.type=0; %0=neumann, 1=robin, 2=dirichlet
bc.left.C=0; % (that data is C in: kdu/dn=C // u+k/hcv*du/dn =C // u=C)
bc.rite.type=1;
bc.rite.C=T_inf;
dat.bc_for=bc; 
% create perturbations
dat.bc_for.left.dC = pert.bc_L*dat.bc_for.left.C;
dat.bc_for.rite.dC = pert.bc_R*dat.bc_for.rite.C;
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

npar.max_nl_iter=30;
npar.nl_tol=1e-6;

return
end

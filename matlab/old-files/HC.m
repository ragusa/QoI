function F=HC

% clear the console screen
clc; close all; 
% make the problem-data a global variable
global dat
% load the data structure with info pertaining to the physical problem
dat.diff=@cond_for;
dat.siga=@nothing;
dat.esrc=@esrc;
dat.pb_type='hc';
if(strcmp(dat.pb_type,'rad'))
    dat.width=40;
    dat.alpha=4;
else
    dat.width=0.4;
end
bc.left.type=0; %0=neumann, 1=robin, 2=dirichlet
bc.left.C=0; % (that data is C in: -Ddu/dn=C // u/4+D/2du/dn=C // u=C)
bc.rite.type=2;
bc.rite.C=8020;
dat.bc=bc; clear bc;

% load the numerical parameters, npar, structure pertaining to numerics
% number of elements
npar.nel = 50;
% domain
npar.x = linspace(0,dat.width,npar.nel+1);
% polynomial degree
npar.porder=1;
% nbr of dofs per variable
npar.ndofs = npar.porder*npar.nel+1;
% for newton's solve
npar.max_newton = 250;
npar.atol_newton = 1e-11;
npar.rtol_newton = 1e-11;
npar.tol_newton_lin = 1e-9;
% 1=numjac  + LU
% 2=numjac  + Gmres
% 3=numjac  + Precond Gmres
% 4=matfree + Gmres
% 5=matfree + Precond Gmres
npar.newton_solve_option = 3;
myoptP=1; 
optP=0;  if(npar.newton_solve_option==3 || npar.newton_solve_option==5), optP=myoptP; end
npar.prec_opt=optP;
% scale variables
npar.scale=1;
dat.bc.left.C=dat.bc.left.C/npar.scale;
dat.bc.rite.C=dat.bc.rite.C/npar.scale;

% connectivity
gn=zeros(npar.nel,npar.porder+1);
gn(1,:)=linspace(1,npar.porder+1,npar.porder+1);
for iel=2:npar.nel
    gn(iel,:)=[gn(iel-1,end) , gn(iel-1,2:end)+npar.porder ];
end
npar.gn=gn; clear gn;

%%%%
%%%% solve nonlinear forward problem
%%%%
dat.adjoint=false;
[u,uf,xf] = solve_forward(npar);
% plot
figure(1)
x=linspace(0,dat.width,npar.ndofs);
plot(x,u,'.'); hold on
plot(xf,uf,'-'); hold all
%%% compute QoI:
QoI_forward = QoI(@response,u,npar)

resi = compute_residual(u,npar); 

%%%%
%%%% solve nonlinear forward problem analytically
%%%%
Tc=T_analytical(0)
QoI_analytical=(quad(@T_analytical,0,0.4))/dat.width

dat.diff=@nothing;
resi2 = compute_residual(u,npar); resi2=-resi2;

%%%%
%%%% solve adjoint problem
%%%%
dat.adjoint=true;
dat.diff_adjoint=@cond_adj;
dat.esrc_adjoint=@response;
dat.forward_sol=u;
dat.bc.rite.C=13;
[us,usf,xsf] = solve_adjoint(npar);
% plot
figure(2)
plot(x,us,'.r'); hold on
plot(xf,usf,'r-'); hold all
%%% compute QoI
QoI_adjoint = QoI(dat.esrc,us,npar)

resia = compute_residual(us,npar); 
epsi=1e-2;
pert=zeros(npar.ndofs,1);pert(end)=epsi;
resia_pert = compute_residual(us+pert,npar);
alpha=((resia_pert-resia)/epsi);alpha=alpha(end-1);
Ja = compute_jacobian(us,npar,resia);

dat.diff_adjoint=@nothing3;
resia2 = compute_residual(us,npar); resia2=-resia2;

% forward
h=(dat.width/npar.nel);
resia2'*u+ h/dat.width/2*u(end)
% resia2(end)=resia2(end) + h/dat.width/2;
% resia2'*u

resi3=resi2;
resi3(end-1)=resi3(end-1)-alpha*u(end);
% resi3(end-1)=resi3(end-1)-alpha*us(end);
resi3'*us +h/dat.width/2*u(end) +alpha*us(end)*u(end-1)
 
% QoI_forward-QoI_adjoint

% % % boundary terms
% % [b,dbdx] =feshpln([1],npar.porder);
% % gn    = npar.gn;
% % nel   = npar.nel;
% % % jacobian of the transformation to the ref. element
% % x0=npar.x(nel);
% % x1=npar.x(nel+1);
% % Jac=(x1-x0)/2;
% % % local values
% % local_T         = b(:,:)    * u(gn(nel,:));
% % local_Ts        = b(:,:)    * us(gn(nel,:));
% % local_dTsdx     = dbdx(:,:) * us(gn(nel,:));
% % d=dat.diff([1],local_T)/Jac;
% % local_T*d*local_dTsdx


% % % % verification is always good
% % % u(1)
% % %(5/2*3e11*40^2+400^5)^0.2 % rad
% % %3e4*0.4^2/2/5+400         % cond
% % 
% % % compute the L2 norm of the error
% % % create a high-order quadrature
% % [xq,wq] = GLNodeWt(25);
% % % store shapeset
% % [b,dbdx] = feshpln(xq,npar.porder);
% % local_err = 0;
% % a=[];
% % aa=[];
% % % loop over elements
% % for iel=1:npar.nel
% %     % element extremities
% %     x0=npar.x(iel);
% %     x1=npar.x(iel+1);
% %     % jacobian of the transformation to the ref. element
% %     Jac=(x1-x0)/2;
% %     % get x values in the interval
% %     x=(x1+x0)/2+xq*(x1-x0)/2;
% %     % b and dbdx(:,:) are of size (nbr of xq values) x (porder+1),
% %     % T is of length porder+1, 2/dx is the 1d jacobian
% %     % local_T is of length the # of points in the quadrature
% %     local_T    = b(:,:) * u(npar.gn(iel,:));
% %     if(strcmp(dat.pb_type,'rad'))
% %         local_Exa = ((dat.alpha+1).*(dat.esrc(x,local_T)).*(dat.width^2-x.^2)/2+dat.bc.rite.C^(dat.alpha+1)).^(1/(dat.alpha+1));
% %         a=[a local_Exa'];
% %         aa=[aa x'];
% %     else
% %         kdata=[1.05 2150 200 1];
% %         options = optimset('Display','off');
% %         for i=1:length(x)
% %             ci= @(z)  kdata(1)*(z-dat.bc.rite.C) + kdata(2)* log( (kdata(3)+kdata(4)*z    ) ...
% %                 ./(kdata(3)+kdata(4)*dat.bc.rite.C) )  -3e4*(dat.width^2-x(i)^2)/2;
% %             local_Exa(i,1)=fsolve(ci,dat.bc.rite.C,options);
% %         end
% %         a=[a local_Exa'];
% %         aa=[aa x'];
% %     end
% %     % compute local error
% %     local_err =  local_err + Jac*dot(wq , (local_T-local_Exa).^2 ) ;
% % end    
% % plot(aa,a,'r-')
% % % L2 norm
% % local_err=sqrt(local_err)
% % % L2 norm divided by mesh size ^ (p+1)
% % x0=npar.x(iel);x1=npar.x(iel+1);
% % local_err/(x1-x0)^(npar.porder+1)

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [T,Tf,xf]=solve_forward(npar)

% perform the Newton solve
T = newton(npar);

% finer representation for plotting
[Tf,xf] = finer(T,npar);

% resi = compute_residual(T,npar);
% J = compute_jacobian(T,npar,resi);
% 
% global dat 
% 
% siga_save=dat.siga; dat.siga=@nothing;
% diff_save=dat.diff; dat.diff=@nothing;
% rhs = compute_residual(T,npar);
% 
% dat.siga=siga_save;
% dat.diff=diff_save;

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ta,Taf,xf]=solve_adjoint(npar)

% perform the Newton solve
Ta = newton(npar);

% finer representation for plotting
[Taf,xf] = finer(Ta,npar);

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  T=newton(npar)

global dat

% initial guess
if ~dat.adjoint
    T=ones(npar.ndofs,1)*400;
    T=linspace(1000,dat.bc.rite.C,npar.ndofs)';
else
    T=zeros(npar.ndofs,1);
    T=linspace(1000,dat.bc.rite.C,npar.ndofs)';
end
T=T/npar.scale;

% compute the residual
resi = compute_residual(T,npar);
% if dat.adjoint
%      J = compute_jacobian(T,npar,resi);
%      dat.diff_adjoint=@nothing3;
%      rhs = compute_residual(T,npar);
%      dat.diff_adjoint=@cond_adj;
% end
% compute its norm
norm_resi=norm(resi);
% compute stopping criterion
tol_newton = npar.atol_newton + npar.rtol_newton*norm_resi;

eigstudy=false;
if(npar.newton_solve_option>3),eigstudy=false;end

inewton=0;
while (inewton<npar.max_newton && norm_resi>tol_newton)
    inewton=inewton+1;
    fprintf('Newton iteration %i, ',inewton);
    % compute the jacobian matrix
    J = compute_jacobian(T,npar,resi);
    P = compute_precond(T,npar,npar.prec_opt);
    if(eigstudy)
%         P = compute_precond(T,npar,dat,1);
        figure(2)
        subplot(2,1,1);
        eigplot(full(J),2);
        subplot(2,1,2);
        JiP=J*inv(P);
        eigplot(full(JiP),2);
        fprintf('\n\tEst. of cond number for J: %g, \t for JiP: %g\n',condest(J),condest(JiP));
        figure(3);
        subplot(2,1,1);   spy(J);
        subplot(2,1,2);   spy(JiP);
    end
    % newton solve
    switch npar.newton_solve_option
        case {1}
            dT = -J\resi;
        case{2,3}
            % gmres(A,b,restart,tol,maxit,M1,M2,x0)
            restart=[];restart=ceil(length(J)/3);
            tol_newton_lin=npar.tol_newton_lin;
            [dT,flag,relres,iter,resvec] = gmres(J,-resi,restart,tol_newton_lin,[],P,[],[]);
            switch flag
                case{0}
                    fprintf('\tgmres converged');
                case{1}
                    fprintf('\tgmres did not converge');
                case{2}
                    fprintf('\tgmres: ill-conditioned');
                case{3}
                    fprintf('\tgmres stagnation');
            end
            fprintf('\titer statistics in gmres: %i, %i\n',iter(1),iter(2));
        case {4,5}
            restart=10;
            tol_newton_lin=npar.tol_newton_lin;
            [w,flag,relres,iter,resvec] = gmres(@Jv,-resi,restart,tol_newton_lin,150,[],[],[],... % the next line is for @afun arguments
                T,resi,npar,P);
            switch flag
                case{0}
                    fprintf('\tgmres converged');
                case{1}
                    fprintf('\tgmres did not converge');
                case{2}
                    fprintf('\tgmres: ill-conditioned');
                case{3}
                    fprintf('\tgmres stagnation');
            end
            fprintf('\titer statistics in gmres: %i, %i\n',iter(1),iter(2));
            dT=P\w;    
    end
    % new Newton iterate
    T=T+dT;        
    % compute the residual
    resi = compute_residual(T,npar);
    % compute its norm
    norm_resi=norm(resi);
    if(norm_resi<tol_newton)
        fprintf(' CONVERGED with ');
    end
    fprintf('Nonlinear residual error=%g, delta solution %g \n\n',norm_resi,norm(dT));
end

% re-scale
T=T*npar.scale;

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F]=compute_residual(T,npar)
% assemble the residual and apply BC

% make the problem-data a global variable
global dat

% shortcuts
porder= npar.porder;
gn    = npar.gn;
nel   = npar.nel;
% n: linear system size
n=npar.ndofs;
% allocate memory
F=zeros(n,1);

% compute local matrices
% load Gauss Legendre quadrature (GLQ is exact on polynomials of degree up to 2n-1,
% while using only integrand evaluations (n-point quadrature).
% estimate the max poly order (it comes from the mass matrix  when coef. are
% piecewise constant and the material mesh matches the computational mesh
[xq,wq] = GLNodeWt(2*(porder+1));
% initialize local residual vector
local_res=zeros(porder+1,1);
% store shapeset
[b,dbdx] =feshpln(xq,porder);

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
    % T is of length porder+1, 2/dx is the 1d jacobian
    local_T         = b(:,:)    * T(gn(iel,:));
    local_dTdx      = dbdx(:,:) * T(gn(iel,:));
    s=dat.siga(npar.scale*local_T,x)*Jac;
    if ~dat.adjoint
        d=dat.diff(npar.scale*local_T,x)/Jac;
        q=dat.esrc(npar.scale*local_T,x)*Jac/npar.scale;
    else
        local_T_forward = b(:,:) * dat.forward_sol(gn(iel,:));
        d=dat.diff_adjoint(npar.scale*local_T,x,local_T_forward)/Jac;
        q=dat.esrc_adjoint(npar.scale*local_T,x)*Jac/npar.scale;
    end
    % compute local residual
    for i=1:porder+1
        local_res(i) =  dot(s.*wq.*b(:,i)    , local_T(:)) ...
            + dot(d.*wq.*dbdx(:,i) , local_dTdx(:)) ...
            - dot(q.*wq, b(:,i));
    end
    % assemble
    F(gn(iel,:)) = F(gn(iel,:)) + local_res;
end

% apply natural BC
Dirichlet_nodes=[];
Dirichlet_val=[];
switch dat.bc.left.type
    case 0 % Neumann, int_bd_domain (b D grad u n) is on the RHS
        F(1)=F(1)-dat.bc.left.C;
    case 1 % Robin
        F(1)=F(1)+1/2*T(1);
        F(1)=F(1)-2*dat.bc.left.C;
    case 2 % Dirichlet
        Dirichlet_nodes=[Dirichlet_nodes 1];
        Dirichlet_val=[Dirichlet_val dat.bc.left.C];
end
switch dat.bc.rite.type
    case 0 % Neumann, int_bd_domain (b D grad u n) is on the RHS
        F(n)=F(n)-dat.bc.rite.C;
    case 1 % Robin
        F(n)=F(n)+1/2*T(n);
        F(n)=F(n)-2*dat.bc.rite.C;
    case 2 % Dirichlet
        Dirichlet_nodes=[Dirichlet_nodes n];
        Dirichlet_val=[Dirichlet_val dat.bc.rite.C];
end
% apply Dirichlet BC
for i=1:length(Dirichlet_nodes);% loop on the number of constraints
    id=Dirichlet_nodes(i);      % extract the dof of a constraint
    bcval=Dirichlet_val(i);
    F(id)=T(id)-bcval;         % put the constrained value in the rhs
end

return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J=compute_jacobian(x,npar,F)

% make the problem-data a global variable
global dat

% compute Jacobian
n=length(x);
nnz_=(2*npar.porder+1)*npar.nel; % upperbound
J=spalloc(n,n,nnz_);
% choose epsi
epsi=choose_epsi(2,x);
for i = 1:n
    if((i==1 && dat.bc.left.type==2)||(i==n && dat.bc.rite.type==2))
        J(:,i)=0; J(i,i)=1;
        continue
    end
    % choose a component in vector epsi (if it is a vector)
    if(length(epsi)==1)
        i_=1;
    else
        i_=i;
    end
    % perturbation to build numerical Jacobian column by column
    v=zeros(n,1);
    v(i)   = 1;
    x_pert = x + epsi(i_)*v;
    % evaluate perturbed residual
    Feps = compute_residual(x_pert,npar);
    J(:,i) = (Feps - F) /epsi(i_);
end
% add transient + SS contributions
% J=speye(nn)-0.5*dt*J;

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [out]=Jv(xkrylov,T,F,npar,P)

% make the problem-data a global variable
global dat

if(dat.bc.left.type==2)
    xkrylov(1)=0;
end
if(dat.bc.rite.type==2)
    xkrylov(end)=0;
end
% % choose epsi
% epsi=choose_epsi(3,T,norm(xkrylov));
% if(length(epsi)>1), error('epsi in Jv should be of length 1'); end
% preconditioning step
xkrylov=P\xkrylov;
% choose epsi
epsi=choose_epsi(3,T,norm(xkrylov));
if(length(epsi)>1), error('epsi in Jv should be of length 1'); end
% perturb newton state
T_pert = T + epsi*xkrylov;
% perturbed residual
Feps = compute_residual(T_pert,npar);
% FD formula
out =(Feps - F) /epsi;

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function epsi=choose_epsi(opt,vec,norm_kry)

if(nargin==3)
    if(norm_kry<eps),norm_kry=1; end
end
b=sqrt(eps);
switch opt
    case 0
        epsi = b; % returns a scalar
    case 1
        epsi = (1+vec)*b; % returns a vector
    case 2 % column by column
        epsi=b*max(abs(vec),1).*sign_(vec); % returns a vector
    case 3
        if(nargin<3), error('Not enough argument in choose_epsi'); end
        epsi=b/(length(vec)*norm_kry)*norm(1+vec,1); % returns a scalar
end

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sg=sign_(x)

sg=sign(x);
if(x==0), sg=1; end

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A]=compute_precond(T,npar,optP)
% assemble the P matrix, with various options:
% optP=0: Identity
% optP=1: compute k(T) at quad points
% optP=1a: compute k(T) at quad points using T^(0)
% optP=2: compute kave per element
% optP=2a: compute kave per element using T^(0)
% the options with a are actually managed outside of this function

% make the problem-data a global variable
global dat

if(optP==0)
    A=speye(length(T));
    return
end

% shortcuts
porder= npar.porder;
gn    = npar.gn;
nel   = npar.nel;
% ideally, we would analyze the connectivity to determine nnz
nnz=(porder+1)*nel; %this is an upperbound, not exact
% n: linear system size
n=npar.ndofs;
% allocate memory
A=spalloc(n,n,nnz);
% initialize local matrices/vectors
m=zeros(porder+1,porder+1);
k=m;
local_mat=m;

% compute local matrices
% load Gauss Legendre quadrature (GLQ is exact on polynomials of degree up to 2n-1,
% while using only integrand evaluations (n-point quadrature).
% estimate the max poly order (it comes from the mass matrix  when coef. are
% piecewise constant and the material mesh matches the computational mesh
if(optP==1)
    poly_max=2*(porder+1);
    [xq,wq] = GLNodeWt(poly_max);
    % store shapeset
    [b,dbdx] =feshpln(xq,porder);
else
    poly_max=porder+1;
    [xq,wq] = GLNodeWt(poly_max);
    % store shapeset
    [b,dbdx] =feshpln(xq,porder);
    % compute local matrices + load vector
    for i=1:porder+1
        for j=1:porder+1
            m(i,j)= dot(wq.*b(:,i)    , b(:,j));
            k(i,j)= dot(wq.*dbdx(:,i) , dbdx(:,j));
        end
    end
    % increase quadrature order to compute kave over each element
    poly_max=2*(porder+1);
    [xq,wq] = GLNodeWt(poly_max);
    % store shapeset
    [b,dbdx] =feshpln(xq,porder);
end

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
    if(optP==1)
        % b and dbdx(:,:) are of size (nbr of xq values) x (porder+1),
        local_T    = b(:,:) * T(gn(iel,:));
        s=dat.siga(npar.scale*local_T,x);
        if ~dat.adjoint
            d=dat.diff(npar.scale*local_T,x)/Jac;
        else
            local_T_forward = b(:,:) * dat.forward_sol(gn(iel,:));
            d=dat.diff_adjoint(npar.scale*local_T,x,local_T_forward)/Jac;
        end
        % compute local matrices
        for i=1:porder+1
            for j=1:porder+1
                m(i,j)= dot(s.*wq.*b(:,i)    , b(:,j));
                k(i,j)= dot(d.*wq.*dbdx(:,i) , dbdx(:,j));
            end
        end
        local_mat = m*Jac + k/Jac;
    else
        % first compute ave_cond
        % b and dbdx(:,:) are of size (nbr of xq values) x (porder+1),
        if ~dat.adjoint
            local_T    = b(:,:) * T(gn(iel,:));
            d=dat.diff(npar.scale*local_T,x);
        else
            local_T_forward = b(:,:) * dat.forward_sol(gn(iel,:));
            d=dat.diff_adjoint(npar.scale*local_T,x,local_T_forward)/Jac;
        end
        ave_cond = dot(wq, d)/sum(wq);
        % then compute local mat
        local_mat =  0*m*Jac + ave_cond*k/Jac;
    end
    % assemble
    A(gn(iel,:),gn(iel,:)) = A(gn(iel,:),gn(iel,:)) + local_mat;
end

% apply natural BC
Dirichlet_nodes=[];
Dirichlet_val=[];
switch dat.bc.left.type
    case 0 % Neumann, int_bd_domain (b D grad u n) is on the RHS
%         rhs(1)=rhs(1)+dat.bc.left.C;
    case 1 % Robin
        A(1,1)=A(1,1)+1/2;
%         rhs(1)=rhs(1)+2*dat.bc.left.C;
    case 2 % Dirichlet
        Dirichlet_nodes=[Dirichlet_nodes 1];
        Dirichlet_val=[Dirichlet_val dat.bc.left.C];
end
switch dat.bc.rite.type
    case 0 % Neumann, int_bd_domain (b D grad u n) is on the RHS
%         rhs(n)=rhs(n)+dat.bc.rite.C;
    case 1 % Robin
        A(n,n)=A(n,n)+1/2;
%         rhs(n)=rhs(n)+2*dat.bc.rite.C;
    case 2 % Dirichlet
        Dirichlet_nodes=[Dirichlet_nodes n];
        Dirichlet_val=[Dirichlet_val dat.bc.rite.C];
end
% apply Dirichlet BC
for i=1:length(Dirichlet_nodes);% loop on the number of constraints
    id=Dirichlet_nodes(i);      % extract the dof of a constraint
    bcval=Dirichlet_val(i);
%     rhs=rhs-bcval*A(:,id);  % modify the rhs using constrained value
    A(id,:)=0; % set all the id-th row to zero
    A(:,id)=0; % set all the id-th column to zero (symmetrize A)
    A(id,id)=1;            % set the id-th diagonal to unity
%     rhs(id)=bcval;         % put the constrained value in the rhs
end

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [QoI]=QoI(r,T,npar)
% compute the QoI

% shortcuts
porder= npar.porder;
gn    = npar.gn;
nel   = npar.nel;
% n: linear system size
n=npar.ndofs;

% compute local matrices
% load Gauss Legendre quadrature (GLQ is exact on polynomials of degree up to 2n-1,
% while using only integrand evaluations (n-point quadrature).
% estimate the max poly order (it comes from the mass matrix  when coef. are
% piecewise constant and the material mesh matches the computational mesh
[xq,wq] = GLNodeWt(2*(porder+1));
% initialize local residual vector
local_res=zeros(porder+1,1);
% store shapeset
[b,dbdx] =feshpln(xq,porder);

QoI=0;
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
    % T is of length porder+1, 2/dx is the 1d jacobian
    local_T    = b(:,:) * T(gn(iel,:));
    q=r(x,local_T)*Jac;
    % compute local response
    QoI = QoI +  dot(q.*wq , local_T(:)) ;
end

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Tf,xf]=finer(T,npar)

cut=10;
% shortcuts
porder= npar.porder;
gn    = npar.gn;
nel   = npar.nel;
% n: linear system size
n=nel*(porder+cut)+1;
% allocate memory
F=zeros(n,1);

xq=linspace(-1,1,porder+1+cut);
[b,dbdx] =feshpln(xq,porder);
% loop over elements
xf=[]; Tf=[];
for iel=1:npar.nel
    % element extremities
    x0=npar.x(iel);
    x1=npar.x(iel+1);
    % get x values in the interval
    x=(x1+x0)/2+xq*(x1-x0)/2;
    % b and dbdx(:,:) are of size (nbr of xq values) x (porder+1),
    % T is of length porder+1, 2/dx is the 1d jacobian
    local_T    = b(:,:) * T(gn(iel,:));
    if(iel==1)
        xf=[xf x];
        Tf=[Tf; local_T];
    else
        xf=[xf x(2:end)];
        Tf=[Tf; local_T(2:end)];
    end
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
function y=cond_for(T,x)
% conductivity function. may depend on space and temperature

% make the problem-data a global variable
global dat

if(strcmp(dat.pb_type,'rad'))
    y=T.^dat.alpha;
else
    kdata=[1.05 2150 200 1];
    % kdata=[5 0 200 1];
    y = kdata(1) + kdata(2)./( kdata(3)+kdata(4)*T);
end

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=cond_adj(T,x,Tforward)
% conductivity function used in adjoint eqs.
y=cond_for(Tforward,x);
% y(:)=1;
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=diffu(T,x)
y=1+0*0.5*x.^2;
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=siga(T,x)
y=x;
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=nothing(~,~)
y=0;
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=nothing3(~,~,~)
y=0;
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=esrc(T,x)

% make the problem-data a global variable
global dat

if(strcmp(dat.pb_type,'rad'))
    y=3e11;
else
    y=3e4;
end

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=response(T,x)
global dat
y=1/dat.width;
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T=T_analytical(xx)
global dat
for i=1:length(xx)
    x=xx(i);
    Tc=1200;
    dT=100;
    while (abs(dT)>1e-6)
        int=quad(dat.diff,dat.bc.rite.C,Tc)-3e4*(dat.width^2-x^2)/2;
        k=dat.diff(Tc,1);
        dT=-int./k;
        Tc=Tc+dT;
    end
    T(i)=Tc;
end
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

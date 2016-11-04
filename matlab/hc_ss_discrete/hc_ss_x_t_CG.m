function hc_ss_x_t_CG
%solver for time dependent PDE using backwards Euler. Includes a
%time-dependent heat equation and a simple du/dt=t case for testing.

% clear the console screen
clc; close all;
 
% load the data and numerical parameters
pert_k    = 1e-1 *1;
pert_s    = 1e-1 *0;
pert_bc_L = 1e-1 *0;
pert_bc_R = 1e-1 *0;
[dat,npar] = load_simulation_data(pert_k,pert_s,pert_bc_L,pert_bc_R);


% assemble the matrix and the rhs
npar.adjoint=false;
npar.pert_status = 'unperturbed';
[Tu,qu]=solve_time_system(npar,dat);
% 
% assemble the matrix and the rhs
npar.adjoint=true;
npar.pert_status = 'unperturbed';
[Phiu,ru]=solve_time_system(npar,dat);
% 
npar.adjoint=false;
npar.pert_status = 'perturbed';
[Tp,qp]=solve_time_system(npar,dat);

QoIArray=[];
for ii=1:length(Tu(1,:))
    QoIArray=[QoIArray dot(Tu(:,ii ), ru(:,ii))];
end
%QoIArray
TuQoI=sum(QoIArray)
QoIArray=[];
for ii=1:length(Tp(1,:))
    QoIArray=[QoIArray dot(Tp(:,ii ), ru(:,ii))];
end
TpQoI=sum(QoIArray)

switch dat.systemSelector
    case {1} %Heat Equation
        figure(1) 
        surf(Tu)
        figure(2) 
        surf(Phiu)
        figure(3) 
        surf(Tp)
        figure(4)
        surf(ru)

        figure(97);hold all;
        for k=1:length(Tu(1,:))
            plot(Tu(:,k));
        end
        figure(98);hold all;
        for k=1:length(Phiu(1,:))
            plot(Phiu(:,k));
        end
        figure(99);hold all;
        for k=1:length(Tp(1,:))
            plot(Tp(:,k));
        end
    case {2} 
        figure(1);
        plot(Tu);
end
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [final,finalSource]=solve_time_system(npar,dat)
%Store time data in dat
dat.finalTime=100;
dat.timeSteps=10;
dat.dt=dat.finalTime/dat.timeSteps;
%Set 
max_steps=dat.timeSteps;
dt=dat.dt;
tfinal=dat.finalTime;
if npar.adjoint==false
    [Tprev,solutionForward,qArray]=getInitial(dat,npar);
    for ii=1:max_steps
        t=ii*dt;
        [T,q]=solveOneStep(npar,dat,t,Tprev);
        solutionForward=[solutionForward T];
        qArray=[qArray q];
        Tprev=T;
    end
    final=solutionForward;
    finalSource=qArray;
end
if npar.adjoint==true
    [Phiprev,solutionAdjoint,rArray]=getInitial(dat,npar);
    t = tfinal;
    for ii=1:max_steps
        t = t -dt;
        [Phi,r]=solveOneStep(npar,dat,t,Phiprev);
        solutionAdjoint=[Phi solutionAdjoint];
        rArray=[r rArray];
        Phiprev=Phi;
    end
    final=solutionAdjoint;
    finalSource=rArray;
end
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [initialGuess,solutionArray,sourceArray]=getInitial(dat,npar)
max_steps=dat.timeSteps;
dt=dat.dt;
tfinal=dat.finalTime;
switch dat.systemSelector
    case {1} %Heat Equation
        if npar.adjoint==false
            initialGuess=ones(npar.ndofs,1);
            [~,initialSource]=solveOneStep(npar,dat,0,initialGuess);
        end
        if npar.adjoint==true
            initialGuess=zeros(npar.ndofs,1);
            [~,initialSource]=solveOneStep(npar,dat,tfinal,initialGuess);
        end
    solutionArray = zeros(npar.ndofs,1);
    sourceArray = zeros(npar.ndofs,1);
    solutionArray(:,1) = initialGuess;
    sourceArray(:,1)=initialSource;
    case {2} %Test system du/dt=a*t
        if npar.adjoint==false
            initialGuess=0;
            initialSource=0;
        end
        if npar.adjoint==true
            initialGuess=0;
            initialSource=0;
        end
    solutionArray =[initialGuess];
    sourceArray=[initialSource];
    otherwise 
        error('wong test selector in %s',mfilename);
end

end

function [returnValue,realSource]=solveOneStep(npar,dat,t_end,Tprev)
% systemSelector is used to determine which system is to be solved:
% 1 => Transient heat equation
% 2 => du/dt = t
dt=dat.dt;
switch dat.systemSelector
    case {1}
        [A,sourceTerm,Tdir,AN,realSource]=assemble_system(npar,dat,t_end,Tprev);
        T=A\sourceTerm;
        returnValue=T;
    case {2}
        [u,realSource]=test_system(npar,dat,t_end,Tprev);
        returnValue=u;
    otherwise 
        error('wong test selector in %s',mfilename);
end

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u,source]=test_system(npar,dat,t_end,uprev)
%Simple test case of solving du/dt = f(t)

% define the rhs f function
a = 1;
my_f = @(t) a*t;

% solve for new u:
u = dat.dt*my_f(t_end) + uprev;

source=[dat.dt]; %Just a dummy for now
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,rhs, out, A_extremities_before_bc,realSource]=assemble_system(npar,dat,t_end,Tprev)

% assemble the matrix, the rhs, apply BC
dt=dat.dt;

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
M=spalloc(n,n,nnz);
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
                kond{i} = @(x,t) dat.k{i}(x,t)    + dat.dk{i}(x,t);
                src{i}  = @(x,t) dat.fsrc{i}(x,t) + dat.dfsrc{i}(x,t);
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
    d=kond{my_zone}(x,t_end);
    q=src{my_zone}(x,t_end);
    % compute local matrices + load vector
    for i=1:porder+1
        for j=1:porder+1
            k(i,j)= dot(d.*wq.*dbdx(:,i) , dbdx(:,j));
            m(i,j)= dot(wq.*b(:,i),b(:,j));
        end
        f(i)= dot(q.*wq, b(:,i));
    end
    % assemble
    A(gn(iel,:),gn(iel,:)) = A(gn(iel,:),gn(iel,:)) + k/Jac;
    M(gn(iel,:),gn(iel,:)) = M(gn(iel,:),gn(iel,:)) + m*Jac;
    rhs(gn(iel,:)) = rhs(gn(iel,:)) + f*Jac;
end
realSource=rhs;
if npar.adjoint 
    A = M/dt +A;
    rhs = rhs + M*Tprev/dt;
else 
    A = M/dt +A;
    rhs = rhs + M*Tprev/dt;
end

% apply natural BC
Dirichlet_nodes=[];
Dirichlet_val=[];
r_functional=[];
switch bc.left.type
    case 0 % Neumann, int_bd_domain (b D grad u n) is on the RHS
        rhs(1)=rhs(1)+bc.left.C;
    case 1 % Robin
        A(1,1)=A(1,1)+dat.hcv;
        rhs(1)=rhs(1)+dat.hcv*bc.left.C;
    case 2 % Dirichlet
        Dirichlet_nodes=[Dirichlet_nodes 1];
        Dirichlet_val=[Dirichlet_val bc.left.C];
        r_functional=[r_functional rhs(1)];
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
        r_functional=[r_functional rhs(n)];
end
% save pieces of A
A_extremities_before_bc=A(Dirichlet_nodes,:);
% apply Dirichlet BC
for i=1:length(Dirichlet_nodes) % loop on the number of constraints
    id=Dirichlet_nodes(i);      % extract the dof of a constraint
    bcval=Dirichlet_val(i);
    rhs=rhs-bcval*A(:,id);  % modify the rhs using constrained value
    A(id,:)=0; % set all the id-th row to zero
    A(:,id)=0; % set all the id-th column to zero (symmetrize A)
    A(id,id)=1;            % set the id-th diagonal to unity
    rhs(id)=bcval;         % put the constrained value in the rhs
end

if npar.adjoint
    out=r_functional;
else
    out=Dirichlet_val;
end


return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [QoI,volumetric_integral] = qoi_num_eval(fem_sol,npar,dat,AN)

% calculation of the QoI using definition and numerical quadrature
% dx=diff(npar.xf);
% JJJJUUUU = dot( (Tu(1:end-1)+Tu(2:end))/2 , dx )/sum(dx)
% JJJJPPPP = dot( (Tp(1:end-1)+Tp(2:end))/2 , dx )/sum(dx)

% shortcuts
porder= npar.porder;
gn    = npar.gn;
nel   = npar.nel;

% select quadrature
poly_max=2*porder;
[xq,wq] = GLNodeWt(porder+1);
% store shapeset
[b,~] =feshpln(xq,porder);

% use adjoint source (response function)
use_response_function = ~npar.adjoint;

% switch between forward and adjoint calculations
if use_response_function % when the forward solution is passed as argument
    src = dat.asrc;
    % no need for anything else: we use only T and r to compute the QoI !!!
    %     bc  = dat.bc_adj;
    %     kond = dat.k;
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
            %         case 'delta_p'
            %             for i=1:length(dat.k)
            %                 kond{i} = @(T) dat.dka{i}(T)+ dat.dkb{i}(T)+ dat.dkc{i}(T);
            %             end
            %             src = dat.dfsrc;
            %             bc.left.type =  0; % Neumann
            %             bc.rite.type =  0;
            %             bc.left.C =  0;
            %             bc.rite.C =  0;
    end
end

QoI = 0;
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
    q_or_r_qp=src{my_zone}(x);
    % evaluate previous temperature at qp
    fem_sol_qp = b(:,:) * fem_sol(gn(iel,:));
    %     % evaluate
    %     d=kond{my_zone}(fem_sol_qp);
    % add to QoI
    QoI = QoI + dot( wq.*fem_sol_qp ,q_or_r_qp)*Jac;
end

volumetric_integral=QoI; % we are done at this point if the forward solution was passed.

% add missing part to QoI when evaluating it with the adjoint function
% recall 0=neumann, 1=robin, 2=dirichlet
% (that data is C in: kdu/dn=C // u+k/hcv*du/dn =C // u=C)
if ~use_response_function
    fprintf('Adding bc stuff to volumetric integral\n\t npar.adjoint=%g, QoI so far=%g\n',npar.adjoint,QoI);
    % store shapeset at element extremities
    [b,dbdx] = feshpln([-1 1],porder);
    switch npar.pert_status
        case 'unperturbed'
            switch bc.left.type
                case 0 % Neumann
                    fem_sol_extremity = b(:,:) * fem_sol(gn(1,:));
                    fprintf('Neumann left, bc value %g\n\tadd to phi.q: ',bc.left.C);
                    fprintf('%g \n',bc.left.C*fem_sol_extremity(1));
                    QoI = QoI +bc.left.C*fem_sol_extremity(1);
                case 1 % Robin
                    fem_sol_extremity = b(:,:) * fem_sol(gn(1,:));
                    fprintf('Robin left, bc value %g\n\tadd to phi.q: ',bc.left.C);
                    fprintf('%g \n',dat.hcv*bc.left.C*fem_sol_extremity(1));
                    QoI = QoI +dat.hcv*bc.left.C*fem_sol_extremity(1);
                case 2 % Dirichlet
                    if npar.use_matrix_terms
                        QoI = QoI -dot(AN(1,:),fem_sol)*bc.left.C;
                    else
                        iel=1;
                        Jac = (npar.x(iel+1)-npar.x(iel))/2;
                        dfem_sol_extremity = dbdx(:,:) * fem_sol(gn(iel,:)) / Jac;
                        my_zone=npar.iel2zon(iel);
                        d=kond{my_zone}(npar.x(iel));
                        fprintf('dirichlet left, bc value %g\n\tadd to phi.q: ',bc.left.C)
                        fprintf('%g \n',d*bc.left.C*dfem_sol_extremity(1))
                        %                     dfem_sol_extremity(1)
                        QoI = QoI +d*bc.left.C*dfem_sol_extremity(1);
                    end
            end
            switch bc.rite.type
                case 0 % Neumann
                    fem_sol_extremity = b(:,:) * fem_sol(gn(end,:));
                    fprintf('Neumann right, bc value %g\n\tadd to phi.q: ',bc.rite.C);
                    fprintf('%g \n',bc.rite.C*fem_sol_extremity(end));
                    QoI = QoI +bc.rite.C*fem_sol_extremity(end);
                case 1 % Robin
                    fem_sol_extremity = b(:,:) * fem_sol(gn(end,:));
                    fprintf('Robin right, bc value %g\n\tadd to phi.q: ',bc.rite.C);
                    fprintf('%g \n',dat.hcv*bc.rite.C*fem_sol_extremity(end));
                    QoI = QoI +dat.hcv*bc.rite.C*fem_sol_extremity(end);
                case 2 % Dirichlet
                    if npar.use_matrix_terms
                        QoI = QoI -dot(AN(end,:),fem_sol)*bc.rite.C;
                    else
                        iel=npar.nel;
                        Jac = (npar.x(iel+1)-npar.x(iel))/2;
                        dfem_sol_extremity = dbdx(:,:) * fem_sol(gn(iel,:)) / Jac;
                        my_zone=npar.iel2zon(iel);
                        d=kond{my_zone}(npar.x(iel+1));
                        fprintf('dirichlet right, bc value %g\n\tadd to phi.q: ',bc.rite.C);
                        fprintf('%g \n',-d*bc.rite.C*dfem_sol_extremity(end));
                        QoI = QoI -d*bc.rite.C*dfem_sol_extremity(end);
                    end
            end
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
dat.systemSelector=1;


% dimensions
if triga
    % dat.hcv = 1612.414; % W/m^2-C
    dat.hcv = 100; % W/m^2-C
    % dat.width = [0.003175 0.01 0.017415 0.0179195]; % m
    dat.width = [1 3 4 7]; % m
    % mesh resolution per region
    nel_zone = [ 5 5 5 2];
else
    dat.hcv = 100;
    dat.width = 0.5*5;
    nel_zone = [ 100 ];
end

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
    dat.k{1}=@(x,t) 2150/200*(1+0*x+0*t/100);
    dat.fsrc{1}=@(x,t) (t<=40)*(10000*(1+0*x+0*t))+(t>40)*(t<=60)*(10000*(1+(t-40)/20))+(t>60)*(10000*(2));
    %dat.fsrc{1}=@(x,t) (10000*(1+0*x+0*t)); %For debuging
    dat.asrc{1}=@(x,t) (1+0*x+0*t)/dat.width(end);
end
% create perturbations
for i=1:length(dat.k)
    dat.dk{i}    = @(x,t) pert_k*dat.k{i}(x,t)*(1+0*t);
    if isequal(dat.fsrc{i},@zero_function)
        dat.dfsrc{i} = @zero_function;
    else
        dat.dfsrc{i} = @(x,t) pert_s*(1+0*x)*(1+0*t);
    end
end

% forward bc
bc.left.type=1; %0=neumann, 1=robin, 2=dirichlet
bc.left.C=0; % (that data is C in: kdu/dn=C // u+k/hcv*du/dn =C // u=C)
if triga
    bc.rite.type=2;
    bc.rite.C=15;
else
    bc.rite.type=2;
    bc.rite.C=100;
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

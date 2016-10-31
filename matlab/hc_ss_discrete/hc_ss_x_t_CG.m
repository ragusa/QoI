function hc_ss_x_t_CG
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


% assemble the matrix and the rhs
npar.adjoint=false;
npar.pert_status = 'unperturbed';
[Tu]=solve_time_system(npar,dat);
%return
%[Au,qu,Tdiru,AN]=assemble_system(npar,dat);
% solve forward system
%Tu=Au\qu;

% assemble the matrix and the rhs
npar.adjoint=true;
npar.pert_status = 'unperturbed';
[Phiu]=solve_time_system(npar,dat);

npar.adjoint=false;
npar.pert_status = 'perturbed';
[Tp]=solve_time_system(npar,dat);
return
%[Aa,r,r_functional_u,~]=assemble_system(npar,dat);
% solve forward system
phi=Aa\r;

% calculation of the QoI using inner products
J_for_unpert = Tu'*r   + dot(r_functional_u,Tdiru);
J_adj_unpert = phi'*qu + dot(r_functional_u,Tdiru);

% calculation of the QoI using definition and numerical quadrature
npar.adjoint=false;
[J_num_for_unpert,vol_integral_Tr] = qoi_num_eval(Tu,npar,dat,AN);
% [J_num_for_pert,~]   = qoi_num_eval(Tp,npar,dat);
npar.adjoint=true;
npar.use_matrix_terms=false;
[J_num_adj_unpert,vol_integral_Phiq] = qoi_num_eval(phi,npar,dat,AN);
J_num_adj_unpert = J_num_adj_unpert + (npar.porder==1)*dot(r_functional_u,Tdiru);
npar.use_matrix_terms=true;
[J_num_adj_unpert_AN,~] = qoi_num_eval(phi,npar,dat,AN);
J_num_adj_unpert_AN = J_num_adj_unpert_AN + dot(r_functional_u,Tdiru);
% trapezoidal rule
dx=diff(npar.xf);
JJJJUUUU  = dot( (Tu(1:end-1)+Tu(2:end))/2   , dx )/sum(dx);
JJJJaUUUU = dot( (phi(1:end-1)+phi(2:end))/2 , dx )*10000;

fprintf('QoI unperturbed: Forward \n------------------------\n');
fprintf('\tinner: %g\n',J_for_unpert);
fprintf('\t%s\t %g \n','Tu''*r',Tu'*r);
fprintf('\t%s %g \n','dot(r_functional_u,Tdiru)',dot(r_functional_u,Tdiru));
fprintf('\tForward: num. : %g\n',J_num_for_unpert);
fprintf('\tvol.int. Tr: %g\n',vol_integral_Tr);
fprintf('\tForward: Trapezoidal rule: %g\n',JJJJUUUU);

fprintf('QoI unperturbed: Adjoint \n-----------------------\n');
extra_term = -dot(AN,phi)*dat.bc_for.rite.C;
fprintf('\tinner: %g\n',J_adj_unpert);
fprintf('\t%s\t %g \n','phi''*qu',phi'*qu);
fprintf('\t%s %g \n','dot(r_functional_u,Tdiru)',dot(r_functional_u,Tdiru));
fprintf('\tAdjoint: num. : %g\n',J_num_adj_unpert);
fprintf('\tAdjoint: num.AN: %g\n',J_num_adj_unpert_AN);
fprintf('\tAdjoint: vol.int. Phiq: %g\n',vol_integral_Phiq);
fprintf('\tAdjoint: Trapezoidal vol.int.: %g\n',JJJJaUUUU);
fprintf('\tAdjoint: vol.int.+Extra: %g\n',vol_integral_Phiq+extra_term+dot(r_functional_u,Tdiru));

% plot solution
figure(1);
subplot(1,2,1); hold all; plot(npar.xf,Tu,'.-'); legend('Tu');
subplot(1,2,2); hold all; plot(npar.xf,phi,'.-');legend('Adj');

%%%%%%%%%%%% perturbations

% assemble the matrix and the rhs
npar.adjoint=false;
npar.pert_status = 'perturbed';
[Ap,qp,Tdirp,AN]=assemble_system(npar,dat);
% solve forward system
Tp=Ap\qp;


% perturbed QoI computed suing forward solution
J_for_pert   = Tp'*r  + dot(r_functional_u,Tdirp);
JJJJPPPP = dot( (Tp(1:end-1)+Tp(2:end))/2   , dx )/sum(dx);

% calculation of the QoI using definition and numerical quadrature
npar.adjoint=false;
[J_num_for_pert,vol_integral_Tr] = qoi_num_eval(Tp,npar,dat,AN);


fprintf('QoI perturbed: Forward \n----------------------\n');
fprintf('\tinner: %g\n',J_for_pert);
fprintf('\t%s\t %g \n','Tp''*r',Tp'*r);
fprintf('\t%s %g \n','dot(r_functional_u,Tdirp)',dot(r_functional_u,Tdirp));
fprintf('\tForward: num. : %g\n',J_num_for_pert);
fprintf('\tvol.int. Tr: %g\n',vol_integral_Tr);
fprintf('\tForward: Trapezoidal rule: %g\n',JJJJPPPP);

% plot solution
figure (1)
subplot(1,2,1); plot(npar.xf,Tp,'.-'); legend(['Tu';'Tp'],'Location','Best');


dT = Tp-Tu;
dq = qp-qu;
dA = Ap-Au;

dJ_for = dT'*r        +dot(r_functional_u,Tdirp-Tdiru);

dJ_adj       = phi'*(dq - dA*Tu)        +dot(r_functional_u,Tdirp-Tdiru);
dJ_adj_exact = phi'*(dq  -dA*Tu -dA*dT) +dot(r_functional_u,Tdirp-Tdiru);

fprintf('SENSITIVITY: \n----------------\n');
fprintf('\tForward: inner: %g\n',dJ_for);
fprintf('\t%s %g \n','dot(r_functional_u,Tdirp-Tdiru)',dot(r_functional_u,Tdirp-Tdiru));
fprintf('\tAdjoint: inner: %g\n',dJ_adj);
fprintf('\tAdjoint: phi*dq: %g\n',phi'*dq);
fprintf('\tAdjoint: phi*dA*Tu: %g\n',phi'*dA*Tu);
fprintf('\tAdjoint: inner exact: %g\n',dJ_adj_exact);

npar.adjoint=true;
npar.use_matrix_terms=false;
[J_num_adj_pert,vol_integral_Phiq] = qoi_num_eval(phi,npar,dat,AN);
J_num_adj_pert = J_num_adj_pert + (npar.porder==1)*dot(r_functional_u,Tdirp);
npar.use_matrix_terms=true;
[J_num_adj_pert_AN,~] = qoi_num_eval(phi,npar,dat,AN);
J_num_adj_pert_AN = J_num_adj_pert_AN + dot(r_functional_u,Tdirp);
fprintf('\tAdjoint: num. : %g\n',J_num_adj_pert-J_num_adj_unpert*0);
fprintf('\tAdjoint: num.AN: %g\n',J_num_adj_pert_AN-J_num_adj_unpert_AN*0);

% (Tp,Aa.phi) = (Tp,r) = (Ap^{-1}qp,r) = (qp,Ap^{-T}r)

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [final]=solve_time_system(npar,dat)
tfinal=dat.finalTime;
dt=dat.dt;
max_steps=dat.timeSteps;
if npar.adjoint==false
    Tprev=ones(npar.nel*npar.porder+1,1);
    for ii=1:max_steps
        t=ii*dt;
        [A,q,Tdir,AN]=assemble_system(npar,dat,t,Tprev);
        T=A\q;
        pause(.001)
        figure(2)
        thing=strcat('T(x,',num2str(t),')');
        plot(npar.xf,T,'.-'); legend(thing);
        Tprev=T;
    end
    final=T;
end
if npar.adjoint==true
    Phiprev=zeros(npar.nel*npar.porder+1,1);
    for ii=1:max_steps
        t=tfinal-ii*dt;
        [A,r,Tdir,AN]=assemble_system(npar,dat,t,Phiprev);
        Phi=A\r;
        pause(.05)
        figure(3)
        thing=strcat('Phi(x,',num2str(t),')');
        plot(npar.xf,Phi,'.-'); legend(thing);
        Phiprev=Phi;
    end
    final=Phi;
end
return
end

function [A,rhs, out, A_extremities_before_bc]=assemble_system(npar,dat,t,Tprev)

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
B=spalloc(n,n,nnz);
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
    d=kond{my_zone}(x,t);
    q=src{my_zone}(x,t);
    % compute local matrices + load vector
    for i=1:porder+1
        for j=1:porder+1
            k(i,j)= dot(d.*wq.*dbdx(:,i) , dbdx(:,j));
            m(i,j)=dot(wq.*b(:,i),b(:,j));
        end
        f(i)= dot(q.*wq, b(:,i));
    end
    % assemble
    A(gn(iel,:),gn(iel,:)) = A(gn(iel,:),gn(iel,:)) + k/Jac;
    B(gn(iel,:),gn(iel,:)) = B(gn(iel,:),gn(iel,:)) + m;
    rhs(gn(iel,:)) = rhs(gn(iel,:)) + f*Jac;
end
if npar.adjoint 
    A=A+(B/dt);
else 
    A=A+(B/dt);
end
rhs=rhs+(Tprev/dt);

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

dat.finalTime=100;
dat.timeSteps=100;
dat.dt=dat.finalTime/dat.timeSteps;

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

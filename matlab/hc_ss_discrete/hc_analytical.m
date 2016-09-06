function hc_analytical
clc; close all;

% problem definition
a=0; b=0.5;
k=2150/200; q=10000; h=100; T_dir=200; T_inf=50;
% region of interest
x1=a; x2=b;
% qoi response. r is used in compute adjoint. fh_r is used for the QoI
% using fh_T. Always make sure that r and fh_r(x) are consistent!
r=1/(x2-x1);
fh_r = @(x) (x<x1).*0 + (x>=x1&x<x2).*r + (x>=x2).*0;
% 
do_plot = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve forward problem for unperturbed temperature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fh_T,fh_dTdx] = compute_temperature(a,b,k,q,T_dir,h,T_inf);
% plot temperature
if do_plot
    x=linspace(a,b,1000);
    figure(1); subplot(1,2,1); plot(x,fh_T(x)); axis tight; title('Temperature');
end
% compute functional using forward solution
J_for = integral(@(x) fh_T(x).*fh_r(x),x1,x2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve adjoint problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi_dir=0; phi_inf=0; 
[fh_phi,fh_dphidx] = compute_adjoint(a,b,x1,x2,k,r,phi_dir,h,phi_inf);
% plot adjoint
if do_plot
    subplot(1,2,2); plot(x,fh_phi(x)); axis tight; title('Adjoint');
end
% compute functional using adjoint
J_adj_vol = q * integral(fh_phi,a,b);
J_adj_dir = k*T_dir*fh_dphidx(a);
J_adj_rob = h*T_inf*fh_phi(b);
J_adj     = J_adj_vol + J_adj_rob + J_adj_dir;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QoI values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Unperturbed QoI values: \nForward %g\nAdjoint %g\n',J_for,J_adj);
fprintf('Splitting the adjoint values:\n\tvolume  \t%g\n\tRobin bc\t%g\n\tDirichlet bc\t%g\n',...
    J_adj_vol,J_adj_rob,J_adj_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create perturbation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save values
fh_T0 = fh_T; fh_dT0dx = fh_dTdx; J_for0 = J_for;
% save unperturbed parameters
k0=k; q0=q; h0=h; T_dir0=T_dir; T_inf0=T_inf;
% perturbed values
pert = 1e-1*0;
q     = q0     * (1 + pert*1); 
k     = k0     * (1 + pert*1);
h     = h0     * (1 + pert*1);
T_dir = T_dir0 * (1 + pert*1);
T_inf = T_inf0 * (1 + pert*1);
%
q1 = q-q0; h1 = h-h0; k1 = k-k0;
T_dir1=T_dir-T_dir0; T_inf1=T_inf-T_inf0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve forward problem for perturbed temperature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fh_T,fh_dTdx] = compute_temperature(a,b,k,q,T_dir,h,T_inf);
if do_plot
    subplot(1,2,1); hold all; plot(x,fh_T(x)); axis tight; hold off; 
end
J_for = integral(@(x) fh_T(x).*fh_r(x),x1,x2);
fprintf('\nPerturbed QoI values: \nForward %g\n',J_for);
% forward sensitivity
dJ_for = (J_for - J_for0)  ;
% verif
% [(J_for - J_for0) integral(@(x) T(x)-T0(x),x1,x2)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adjoint-based sensitivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sensivitiy due to heat source
dJ_adj_vol_src = q1 * integral(fh_phi,a,b);
% sensivitiy due to conductivity, first-order and exact
dJ_adj_vol_cond       = -k1 * integral(@(x) fh_dT0dx(x).*fh_dphidx(x),a,b);
dJ_adj_vol_cond_exact = -k1 * integral(@(x) fh_dTdx(x) .*fh_dphidx(x),a,b);
dJ_adj_dirichlet      = k0*fh_dphidx(a)*T_dir1;
dJ_adj_robin          = -h1*fh_T0(b)*fh_phi(b) + (h*T_inf-h0*T_inf0)*fh_phi(b);
dJ_adj_robin_exact    = -h1*fh_T(b) *fh_phi(b) + (h*T_inf-h0*T_inf0)*fh_phi(b);

dJ_adj       = dJ_adj_vol_src + dJ_adj_vol_cond       + dJ_adj_dirichlet +dJ_adj_robin      ;
dJ_adj_exact = dJ_adj_vol_src + dJ_adj_vol_cond_exact + dJ_adj_dirichlet +dJ_adj_robin_exact;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sensitivity values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Sensitivity values:\nForward: %g\n',dJ_for);
fprintf('Adjoint:\n');
fprintf('\tAll    \t\t%g\n',dJ_adj);
fprintf('\tAll E  \t\t%g\n',dJ_adj_exact);
fprintf('\tvol.src  \t%g\n',dJ_adj_vol_src);
fprintf('\tvol.cond \t%g\n',dJ_adj_vol_cond);
fprintf('\tvol.condE\t%g\n',dJ_adj_vol_cond_exact);
fprintf('\tDirichlet\t%g\n',dJ_adj_dirichlet);
fprintf('\tRobin bc \t%g\n',dJ_adj_robin);
fprintf('\tRobin bcE\t%g\n',dJ_adj_robin_exact);

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%

function [fh_T,fh_dTdx] = compute_temperature(a,b,k,q,T_dir,h,T_inf)
% Problem statement: -d/dx(k.dT/dx) = q
%                    T(a) = T_dir
%                    -k.dT/dx|_b = h(T(b) - T_inf)
%
% Analytical solution:
%      T(x) = alpha (x-a)(x-b) + beta(x-b) + gamma
%      with alpha = -q/k/2

mat = [ (a-b) 1 ; k  h];
rhs = [T_dir;  q/2*(b-a) + h*T_inf];
coef = mat\rhs;

alpha = -q/k/2; beta = coef(1); gamma = coef(2);
% function handle for temperature
fh_T = @(x) alpha*(x-a).*(x-b) + beta*(x-b) + gamma;
% function handle for derivative of temperature
fh_dTdx = @(x) alpha*(2*x-a-b) + beta;

% varargout = cell(1,nargout);
% % function handle for T and dTdx
% varargout{1} = fh_T;
% varargout{2} = fh_dTdx;
% switch nargout
%     case {1,2}
%         % nothing to do
%     case 3
%         % value at x=a
%         varargout{3} = fh_T(b);
%     case 4
%         % value at x=b
%         varargout{3} = fh_T(b);
%         % value of derivative at at x=a
%         varargout{4} = alpha*(2*a-a-b)+beta;
%     otherwise
%         error('wrong nargout in compute T');
% end        

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%

function [fh_phi, fh_dphidx] = compute_adjoint(a,b,x1,x2,k,r,phi_dir,h,phi_inf)
% Problem statement: -d/dx(k.dphi/dx) = r
%                    phi(a) = phi_dir
%                    -k.dphi/dx|_b = h(phi(b) - phi_inf)
%
% Analytical solution:
% if 3 regions (r(x).ne.0 only inside the domain [a,b]:
%    phi{1}(x) = alpha{1}(x-a) + beta{1}
%    phi{2}(x) = alpha{2}(x-x1)(x-x2) + beta{2}(x-x2) + gamma{2}
%    phi{3}(x) = alpha{3}(x-b) + beta{3}
% if 2 regions (case x1=a)
%    phi{1}(x) = alpha{1}(x-a)(x-x2) + beta{1}(x-a) + gamma{1}
%    phi{2}(x) = alpha{2}(x-b) + beta{2} 
% if 2 regions (case x2=b)
%    phi{1}(x) = alpha{1}(x-a) + beta{1}
%    phi{2}(x) = alpha{2}(x-x1)(x-b) + beta{2}(x-b) + gamma{2}
% if 1 region. same as solving for T routine

my_case='not_defined';

if ( abs(a-x1)<eps && abs(b-x2)<eps )
    my_case = 'both_bd';
end
if ( abs(a-x1)<eps && abs(b-x2)>eps )
    my_case = 'left_bd';
end
if ( abs(a-x1)>eps && abs(b-x2)<eps )
    my_case = 'right_bd';
end
if ( abs(a-x1)>eps && abs(b-x2)>eps )
    my_case = 'inside';
end
% fprintf('my case is %s\n',my_case);

% analytical solution
switch my_case
    case 'inside'
        % set up linear system
        mat = zeros(6,6); rhs = zeros(6,1);
        % coef = [ a1,b1,b2,g2,a3,b3]
        alpha2 = -r/k/2;
        % fill system
        % phi{1}(a) = phi_dir = beta{1}
        mat(1,2) = 1; rhs(1) = phi_dir;
        % -k dphi{3}dx(b) = h (phi{3}(b)-phi_inf)
        % - k a3 = h(b3-phi_inf)
        mat(2,5:6) = [k h ]; rhs(2) = h*phi_inf;
        % interfaces: values
        % phi{1}(x1) = phi{2}(x1)
        % a1(x1-a)+b1 = b2(x1-x2)+g2
        mat(3,1:4) = [(x1-a) 1 -(x1-x2) -1];
        % phi{2}(x2) = phi{3}(x2)
        % g2 = a3(x2-b) +b3
        mat(4,4:6) = [1 -(x2-b) -1];
        % interfaces: derivatives
        % a1 = a2(2x-x1-x2)+b2
        mat(5,1) = 1; mat(5,3) = -1; rhs(5)= alpha2*(x1-x2);
        % a3 = a2(2x-x1-x2)+b2
        mat(6,5) = 1; mat(6,3) = -1; rhs(6)=-alpha2*(x1-x2);
        % solve system
        coef = mat\rhs;
        % assign coefficients
        alpha1 = coef(1);
        beta1  = coef(2);
        beta2  = coef(3);
        gamma2 = coef(4);
        alpha3 = coef(5);
        beta3  = coef(6);
        % define solution per region
        phi{1} = @(x) alpha1*(x-a) + beta1;
        phi{2} = @(x) alpha2*(x-x1).*(x-x2) + beta2*(x-x2) + gamma2;
        phi{3} = @(x) alpha3*(x-b) + beta3;
        % create full-domain solution
        fh_phi = @(x) (x<x1).*phi{1}(x) + (x>=x1&x<x2).* phi{2}(x) + (x>=x2).*phi{3}(x);
        % define solution derivative per region
        dphi{1} = @(x) alpha1;
        dphi{2} = @(x) alpha2*(2*x-x1-x2) + beta2;
        dphi{3} = @(x) alpha3;
        % create full-domain solution
        fh_dphidx = @(x) (x<x1).*dphi{1}(x) + (x>=x1&x<x2).* dphi{2}(x) + (x>=x2).*dphi{3}(x);
        
    case 'left_bd'
        % set up linear system
        mat = zeros(4,4); rhs = zeros(4,1);
        % coef = [ b1,g1,a2,b2]
        alpha1 = -r/k/2;
        % fill system
        % phi{1}(a) = phi_dir = gamma{1}
        mat(1,2) = 1; rhs(1) = phi_dir;
        % -k dphi{2}dx(b) = h (phi{2}(b)-phi_inf)
        % - k a2 = h(b2-phi_inf)
        mat(2,3:4) = [k h ]; rhs(2) = h*phi_inf;
        % interfaces: values
        % phi{1}(x2) = phi{2}(x2)
        % b1(x2-a)+g1 = a2(x2-b)+b2
        mat(3,1:4) = [(x2-a) 1 -(x2-b) -1];
        % interfaces: derivative at x2
        % a1(2x-a-x2)+b1 = a1(x2-a)+b1 = a2
        mat(4,1) = 1; mat(4,3) = -1; rhs(4)= -alpha1*(x2-a);
        % solve system
        coef = mat\rhs;
        % assign coefficients
        beta1  = coef(1);
        gamma1 = coef(2);
        alpha2 = coef(3);
        beta2  = coef(4);
        % define solution per region
        phi{1} = @(x) alpha1*(x-a).*(x-x2) + beta1*(x-a) + gamma1;
        phi{2} = @(x) alpha2*(x-b) + beta2;
        % create full-domain solution
        fh_phi = @(x) (x<x2).*phi{1}(x)  + (x>=x2).*phi{2}(x);
        % define solution derivative per region
        dphi{1} = @(x) alpha1*(2*x-a-x2) + beta1;
        dphi{2} = @(x) alpha2;
        % create full-domain solution
        fh_dphidx = @(x) (x<x2).*dphi{1}(x) + (x>=x2).*dphi{2}(x);
        
    case 'right_bd'
%    phi{1}(x) = alpha{1}(x-a) + beta{1}
%    phi{2}(x) = alpha{2}(x-x1)(x-b) + beta{2}(x-b) + gamma{2}
       % set up linear system
        mat = zeros(4,4); rhs = zeros(4,1);
        % coef = [ a1,b1,b2,g2]
        alpha2 = -r/k/2;
        % fill system
        % phi{1}(a) = phi_dir = beta{1}
        mat(1,2) = 1; rhs(1) = phi_dir;
        % -k dphi{2}dx(b) = h (phi{2}(b)-phi_inf)
        % - k(a2(2x-x1-b)+b2) = h(g2-phi_inf)
        mat(2,3:4) = [k h ]; rhs(2) = -k*alpha2*(b-x1) + h*phi_inf;
        % interfaces: values
        % phi{1}(x1) = phi{2}(x1)
        % a1(x1-a)+b1 = b2(x1-b)+g2
        mat(3,1:4) = [(x1-a) 1 -(x1-b) -1];
        % interfaces: derivatives at x1
        % a1 = a2(2x-x1-b)+b2 = a2(x1-b)+b2
        mat(4,1) = 1; mat(4,3) = -1; rhs(4)= alpha2*(x1-b);
        % solve system
        coef = mat\rhs;
        % assign coefficients
        alpha1 = coef(1);
        beta1  = coef(2);
        beta2  = coef(3);
        gamma2 = coef(4);
        % define solution per region
        phi{1} = @(x) alpha1*(x-a) + beta1;
        phi{2} = @(x) alpha2*(x-x1).*(x-b) + beta2*(x-b) + gamma2;
        % create full-domain solution
        fh_phi = @(x) (x<x1).*phi{1}(x)  + (x>=x1).*phi{2}(x);
        % define solution derivative per region
        dphi{1} = @(x) alpha1;
        dphi{2} = @(x) alpha2*(2*x-x1-b);
        % create full-domain solution
        fh_dphidx = @(x) (x<x1).*dphi{1}(x) + (x>=x1).*dphi{2}(x);
                
    case 'both_bd'
        [fh_phi, fh_dphidx] = compute_temperature(a,b,k,r,phi_dir,h,phi_inf);
    otherwise
        error('unknown case %s',my_case);
end

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%

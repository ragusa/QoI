function hc_nonlinear_analytical
clc; close all;

% problem definition: k(T)=T^q
% -d/dx (T^a dT/dx) = q = -1/(a+1).d^2(T^(a+1))/dx^2
% dTdx|_0=0 T(L)=T_dir
L=0.5; q=10000; T_dir=100; a=1;
% region of interest
% "whole domain" or "dirac at 0"
% 
do_plot = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve forward problem for unperturbed temperature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% int_0^x (eq) = -1/(a+1).d(T^(a+1))/dx = q.x
% int_x^L (new.eq) => T^(a+1)(x) = T_dir^(a+1) + (a+1)/2.q.(L^2-x^2)
%
ap1=a+1;
fh_T = @(x) (T_dir^ap1 + ap1/2*q*(L^2-x.^2)).^(1/ap1);
fh_dTdx = @(x) (-q*x)./(fh_T(x).^a);
% plot temperature
xx=linspace(0,L,10000);
if do_plot
    figure(1); subplot(1,2,1); plot(xx,fh_T(xx)); axis tight; title('Temperature');
end
% compute functional using forward solution
response = 'dirac';
switch response
    case 'dirac'
        % response = dirac at x=0
        J_for = fh_T(0);
    case 'integrated'
        % response = 1 for all x
        J_for = integral(fh_T,0,L);
    case 'averaged'
        % response = 1/L for all x
        J_for = integral(fh_T,0,L)/L;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve adjoint problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bc : phi(L)=phi_dir=0, 
%
% if response = constant, it is easy to see that: 
% phi(x) is proportional to (T(x)-T_dir). phi(x)=C.(T(x)-T_dir)
%  -d/dx (T^a dphi/dx) = r = q.r/r 
% <==> -q/r d/dx (T^a dphi/dx) = q 
% <==> -C.q/r d/dx (T^a dT/dx) = q 
% <==> -C.q/r/(a+1) d^2[T^(a+1)]/dx^2 = q 
% => C = r/q
%
% if response = dirac @ 0 :
% -d/dx (T^k dphi/dx) = dirac(x) = 0
% source condition at x=0: -2T^k(0) dphidx|_0 = 1
% integrate(eq) from x to L: (rhs=0 in that range)
% ==> T^k.dphidx = C. Apply source condition. C = -1/2
% so dphi/dx=-1/2/T^a(x)
% integrate between x and L:
% phi(L) - phi(x) = (-1/2) int_x^L 1/T^a(y) dy = -phi(x)
%
switch response
    case 'dirac'
        for i=1:length(xx)
            phi(i) = integral(@(x) 1./(fh_T(x).^a), xx(i),L)/2;
        end
    case {'integrated','averaged'}
        if strcmpi(response,'integrated')
            rr=1;
        else
            rr=1/L;
        end
        fh_phi = @(x) rr/q*(fh_T(x)-T_dir);
        phi = fh_phi(xx);
end
% plot adjoint
if do_plot
    subplot(1,2,2); plot(xx,phi); axis tight; title('Adjoint');
end
% compute functional using adjoint
% first, the volumetric part:
J_adj_vol = q*dot( (phi(1:end-1)+phi(2:end))/2 , diff(xx) ); 
fprintf('QoI: Adjoint\n  volumetric term: %g\n',J_adj_vol);
if strcmpi(response,'integrated') || strcmpi(response,'averaged')
    J_adj_vol = q*integral(fh_phi,0,L);
    fprintf('QoI: Adjoint\n  volumetric term (exact): %g\n',J_adj_vol);
end
% now, the bc parts:
switch response
    case 'dirac'
        % it is easier to think of this one as integral from -L to +L
        % (1) kdT/dx.phi:  phi(+/-L)=0, so nothing
        % (2) -kdphi/dx.T: kdphi/dx=-1/2 everythwere
        % because of the dirac, the vol+bc contributions calculated as 
        % integral from -L to +L are the QoI, not twice the QoI
       J_adj_dir = T_dir;
       J_adj_vol = 2*J_adj_vol;
    case {'integrated','averaged'}
        % we need to add 
        % (1) kdT/dx.phi. dTdx=0 at x=0+, and phi(L)=0, so nothing
        % (2) -kdphi/dx.T.
        % only term left is: -kdphi/dx.T at x=L, so T_dir^(a+1)*dphi/dx
        % but dphidx = C dT/dx = r/q (-qL T^(-a))
       J_adj_dir = L*rr*T_dir;
end
J_adj = J_adj_vol + J_adj_dir;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QoI values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Unperturbed QoI values: \nForward %g\nAdjoint %g\n',J_for,J_adj);
fprintf('Splitting the adjoint values:\n\tvolume  \t%g\n\tDirichlet bc\t%g\n',...
    J_adj_vol,J_adj_dir);

disp('What is below is a copy-paste form the linear analytical solution');
warning('Amend first if you want to use');
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create perturbation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save values
fh_T0 = fh_T; fh_dT0dx = fh_dTdx; phi0=phi; J_for0 = J_for;
% save unperturbed parameters
a0=a; q0=q; T_dir0=T_dir;
% perturbed values
pert = 1e-1;
q     = q0     * (1 + pert*1); 
a     = a0     * (1 + pert*1);
T_dir = T_dir0 * (1 + pert*1);
%
a1 = a-a0; q1 = q-q0; T_dir1=T_dir-T_dir0; 

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

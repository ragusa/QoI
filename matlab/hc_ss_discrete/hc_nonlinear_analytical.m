function hc_nonlinear_analytical
clc; close all;

% problem definition: k(T)=T^q
% -d/dx (T^a dT/dx) = q = -1/(a+1).d^2(T^(a+1))/dx^2
% dTdx|_0=0 T(L)=T_dir
L=0.5; q=10000; T_dir=0; a=1;
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
xx=linspace(0,L,100);
if do_plot
    figure(1); subplot(1,2,1); hold all;
    plot(xx,fh_T(xx)); axis tight; title('Temperature');
end
% compute functional using forward solution
response = 'averaged';
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
        fh_dphidx = @(x) -1/2*(fh_dTdx(x)).^(-a);
    case {'integrated','averaged'}
        if strcmpi(response,'integrated')
            rr=1;
        else
            rr=1/L;
        end
        fh_phi = @(x) rr/q*(fh_T(x)-T_dir);
        phi = fh_phi(xx);
        fh_dphidx = @(x) rr/q*fh_dTdx(x);
end
% plot adjoint
if do_plot
    subplot(1,2,2); hold all; plot(xx,phi); axis tight; title('Adjoint');
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
a     = a0     * (1 + pert*0);
T_dir = T_dir0 * (1 + pert*0);
%
a1 = a-a0; q1 = q-q0; T_dir1=T_dir-T_dir0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve forward problem for perturbed temperature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fh_T = @(x) (T_dir^ap1 + ap1/2*q*(L^2-x.^2)).^(1/ap1);
fh_dTdx = @(x) (-q*x)./(fh_T(x).^a);
% plot temperature
if do_plot
    figure(1); subplot(1,2,1); plot(xx,fh_T(xx),'r-'); legend(['TempU';'TempP']);
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

fprintf('\nPerturbed QoI values: \nForward %g\n',J_for);
dJ_for = (J_for - J_for0)  ;
fprintf('\nSensitivity, dJ: \nForward %g\n',dJ_for);
% forward sensitivity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adjoint-based sensitivity
% first-order formula:
% dJ = int[dq,phi] + int[dk gradphi gradT] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sensivitiy due to heat source
dJ_adj_vol_src = q1*dot( (phi(1:end-1)+phi(2:end))/2 , diff(xx) ); 
fprintf('QoI: Adjoint\n  volumetric term: %g\n',dJ_adj_vol_src);
if strcmpi(response,'integrated') || strcmpi(response,'averaged')
    dJ_adj_vol_src = q1*integral(fh_phi,0,L);
    fprintf('QoI: Adjoint\n  volumetric term (exact): %g\n',dJ_adj_vol_src);
end
% sensivitiy due to conductivity, first-order and exact
dJ_adj_vol_cond       = -k1 * integral(@(x) fh_dT0dx(x).*fh_dphidx(x),0,L);
dJ_adj_vol_cond_exact = -k1 * integral(@(x) fh_dTdx(x) .*fh_dphidx(x),0,L);
dJ_adj_dirichlet      = -k0*fh_dphidx(L)*T_dir1;
% dJ_adj_robin          = -h1*fh_T0(b)*fh_phi(b) + (h*T_inf-h0*T_inf0)*fh_phi(b);
% dJ_adj_robin_exact    = -h1*fh_T(b) *fh_phi(b) + (h*T_inf-h0*T_inf0)*fh_phi(b);

% dJ_adj       = dJ_adj_vol_src + dJ_adj_vol_cond       + dJ_adj_dirichlet +dJ_adj_robin      ;
% dJ_adj_exact = dJ_adj_vol_src + dJ_adj_vol_cond_exact + dJ_adj_dirichlet +dJ_adj_robin_exact;
dJ_adj       = dJ_adj_vol_src + dJ_adj_vol_cond       + dJ_adj_dirichlet ;
dJ_adj_exact = dJ_adj_vol_src + dJ_adj_vol_cond_exact + dJ_adj_dirichlet ;

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
% fprintf('\tRobin bc \t%g\n',dJ_adj_robin);
% fprintf('\tRobin bcE\t%g\n',dJ_adj_robin_exact);

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


return
end
%%%%%%%%%%%%%%%%%%%%%%%%%

function [fh_phi, fh_dphidx] = compute_adjoint(a,b,x1,x2,k,r,phi_dir,h,phi_inf)
% Problem statement: -d/dx(k.dphi/dx) = r
%                    phi(a) = phi_dir
%                    -k.dphi/dx|_b = h(phi(b) - phi_inf)
%
% Analytical solution:


return
end
%%%%%%%%%%%%%%%%%%%%%%%%%

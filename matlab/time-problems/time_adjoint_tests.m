function time_adjoint_tests
clc
close all
global pbID

% intiial value
pbID=1;
y0=1;


fprintf('Entering reference solution computation \n');
switch pbID
    case {1}
        tend=10;
        % tolerances for odesolvers
        rtol = 1e-13; atol = 1e-13;
        options = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',1e-10,'OutputFcn',@odeplot,'OutputSel',[1]);
        tspan=[0 tend];
        [t_ref,y_ref]=ode15s(@myfunc,tspan,y0,options);

        QoI_trap_ref=dot(diff(t_ref),(y_ref(1:end-1)+y_ref(2:end))/2);
        fprintf('Using the ODE solver data, the QoI using trapezoidal rule is    %g \n',QoI_trap_ref);
        
        tol=1e-10;
        % QoI_quad_ref=quad(@(time)myintegrant(time,y0),0,tend,tol);
        % QoI_quad_ref=3.27793900571645; % tol 1e-8
        QoI_quad_ref=3.27793899125144; % tol 1e-10
        fprintf('Using the quad function and calls to the ODE solver, the QoI is %g  (tol=%g)\n',QoI_quad_ref,tol);

    case {2}
        tend=1;
        t_ref=linspace(0, tend, 100);
        y_ref =(5*exp(2*t_ref)-2*t_ref-1)/4;
        QoI_quad_ref=(5*exp(2)-9)/8;
        QoI_trap_ref=QoI_quad_ref;
        plot(t_ref,y_ref)

%         rtol = 1e-13; atol = 1e-13;
%         options = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',1e-10,'OutputFcn',@odeplot,'OutputSel',[1]);
%         tspan=[0 tend];
%         [t_ref2,y_ref2]=ode15s(@myfunc,tspan,y0,options);
% 
%         QoI_trap_ref2=dot(diff(t_ref2),(y_ref2(1:end-1)+y_ref2(2:end))/2);
%         tol=1e-12;
%         QoI_quad_ref2=quad(@(time)myintegrant(time,y0),0,tend,tol);
%         [QoI_quad_ref QoI_quad_ref2 QoI_trap_ref2]'
%         [QoI_quad_ref QoI_quad_ref2 QoI_trap_ref2]'-QoI_quad_ref
        
end
        
% initial number of time steps 
nsteps0=20;
nbr_time_refinement=10;
for i=1:nbr_time_refinement
    % fprintf('CN evaluation %g out of %g\n',i,nbr_time_refinement);
    nsteps(i)=nsteps0*2^(i-1);
    [dt(i), yforward{i}, ~, ~] = forward_crank_nicholson(nsteps(i),y0,tend);
    QoI_trap(i)=dt(i)*sum((yforward{i}(1:end-1)+yforward{i}(2:end))/2);
    y_err(i)= abs(yforward{i}(end)-y_ref(end));
end
% hold all
% time=(0:nsteps(end))*dt(end);
% plot(time,yforward{end},'+-');legend('exact','best CN')
% hold off

% verify that cn is 2nd order
figure(2); order=2; C=y_err(1)/dt(1)^order/2;
loglog(dt,y_err,'r+-',dt,C*dt.^order,'m-'); axis tight; legend('error',sprintf('slope %g',order))

% plot conv rate of QoI
if pbID==1
    order=1.5; % it is unclear to me why when a(t) and b(t) are function of time, the cv rate goes to 1.5
end
figure(3); C=abs(QoI_trap(1)-QoI_trap_ref)/dt(1)^order;
loglog(dt,abs(QoI_trap'-QoI_trap_ref),'r+-',dt,C*dt.^order,'m-'); axis tight; legend('QoI error',sprintf('slope %g',order))

figure(4); C=abs(QoI_trap(1)-QoI_quad_ref)/dt(1)^order;
loglog(dt,abs(QoI_trap'-QoI_quad_ref),'r+-',dt,C*dt.^order,'m-'); axis tight; legend('QoI error',sprintf('slope %g',order))

QoI_trap'
% QoI_trap'-QoI_trap_ref
% QoI_trap'-QoI_quad_ref

nsteps=nsteps(end);
[dt, yforward, Aforward, bforward] = forward_crank_nicholson(nsteps,y0,tend);
K=speye(nsteps+1);K(1,1)=0.5; K(end,end)=0.5;K=K*dt;
r=ones(nsteps+1,1);

% trapezoidal rule for QoI. should be the same answer as before, 
% just written differently with operators
qoi_ur=dot(yforward,K*r)

% create adjoint operator without caring about the adjoint final conditions
A_adj=K\Aforward'*K;
% solve for adjoint solution
u_adj=A_adj\r;

% verify that the qoi computed the adjoint yields the same answer!
qoi_uadjb=dot(u_adj,K*bforward)

end

%%%%--------------------------------------%%%%
%%%%   forward_crank_nicholson
%%%%  (I-Anew.dt/2)unew = (I+Aold.dt/2)uold + dt/2(bnew+bold)
%%%%--------------------------------------%%%%
function [dt,y, A, bb] = forward_crank_nicholson(nsteps,y0,tend)
% compute time step size
dt=tend/nsteps;
% initialize output structure
y=zeros(nsteps+1,1);y(1)=y0;
if(nargout==2)
    A=[];bb=[]
else
    A=spalloc(nsteps+1,nsteps+1,3*nsteps+2);A(1,1)=1;
    bb=zeros(nsteps+1,1); bb(1)=y0;
end
% loop over time steps
for i=1:nsteps
    % compute time values
    ti  =dt*(i-1);
    tip1=dt*i;
    % compute governing law values
    ai  =a(ti);
    aip1=a(tip1);
    bi  =b(ti);
    bip1=b(tip1);
    % crank nicholson rhs
    rhs=(1+ai*dt/2)*y(i) + dt/2*(bi+bip1);
    y(i+1)=rhs/(1-dt/2*aip1);
    if(nargout~=2)
        A(i+1,i:i+1)=[ -(1+dt/2*ai) (1-dt/2*aip1)];
        bb(i+1)=dt/2*(bi+bip1);
    end
end

end

%%%%--------------------------------------%%%%
%%%%   time dependent function a(t)
%%%%--------------------------------------%%%%
function out=a(time)
global pbID
switch pbID
    case {1}
        out=-2*time;
    case{2}
        out=2;
end
end

%%%%--------------------------------------%%%%
%%%%   time dependent function b(t)
%%%%--------------------------------------%%%%
function out=b(time)
global pbID
switch pbID
    case {1}
        out=sqrt(time);
    case{2}
        out=time;
end
end

%%%%--------------------------------------%%%%
%%%%   rhs of governing law: a*u+b
%%%%--------------------------------------%%%%
function out=myfunc(time,y)
out=a(time)*y+b(time);
end

%%%%--------------------------------------%%%%
%%%%   output u(t) needed for quad integration
%%%%--------------------------------------%%%%
function out=myintegrant(time,y0)
rtol = 1e-13; atol = 1e-13;
options = odeset('RelTol',rtol,'AbsTol',atol,'InitialStep',1e-10);
for i=1:length(time)
    if(time(i)<eps)
        out(i)=y0;
    else
        tspan=[0 time(i)];
        [t,y]=ode15s(@myfunc,tspan,y0,options);
        out(i)=y(end);
    end
end
end


% close all; clear all; clc
function test_dirichlet2
close all;
clc

n=11;

% size
len=10; h=len/(n-1);
% forward
Lx= 2; Rx= 3; vx=1;
% adjoint
Ly= 1; Ry= 11; vy= 1/len;

% u=ax^2+bx+c, u(0)=c=Lx; 2a=-vx; u(L)=(-vx/2)*len^2+b*len+Lx=Rx
alpha=-vx/2; gamma=Lx; beta=(Rx-Lx+vx/2*len^2)/len;
xx=linspace(0,len,n)'; 
z=[-1 1]/sqrt(3);z=(z+1)/2*h;ii=1:n-1;ii=(ii-1)*h;
zz=[ (z(1)+ii ) (z(2)+ii )]; zz=sort(zz);
exa2=alpha*zz.^2+beta*zz+gamma;
xxx=linspace(0,len,100)'; 
exa3=alpha*xxx.^2+beta*xxx+gamma;
QoI_analytical = (6*gamma + len*(3*beta + 2*alpha*len))/6.

[A,b]=make_matrix(n,h,vx);
b_before=b;
[A,b]=apply_bc(A,b,Lx,Rx);
x=A\b;
plot(xx,x); hold all
[As,bs]=symmetrize_bc(A,b,Lx,Rx);
xs=As\bs;
plot(xx,xs,'o')
plot(zz,exa2,'+-')
plot(xxx,exa3,'v-')
% As*x-b-(As-A)*x
% bs = b+(As-A)*x

QoI=(sum(x)-0.5*(x(1)+x(end)))*h/len
% % QoI=0;
% % for i=1:n-1
% %     QoI=QoI+h*(x(i)+x(i+1))/2;
% % end
% % QoI=QoI/len


[M,q]=make_matrix(n,h,vy);
q_before=q;
Ly=q(1);Ry=q(end);
[M,q]=apply_bc(M,q,Ly,Ry);
y=M\q;
plot(xx,y); hold all
[Ms,qs]=symmetrize_bc(M,q,Ly,Ry);
ys=Ms\qs;
plot(xx,ys,'o')

disp('q_before*x')
q_before'*x
dq   =q-q_before;
dq_bc=qs-q;
disp('(dq+dq_bc)*x')
(dq+dq_bc)'*x
disp('y*bs')
y'*bs
y'*bs-(dq+dq_bc)'*x

% [y' *b  x' *q  y'*(As-A)*x x'*(Ms-M)*y]
% [ys'*bs xs'*qs]
% ys(1:end-1)'*bs(1:end-1) + ys(end)*bs(end)
% xs(1:end-1)'*qs(1:end-1) + xs(end)*qs(end)

% aa=ys(1)*(xs(2)-xs(1)) - ys(end)*(xs(end)-xs(end-1))
% bb=xs(1)*(ys(2)-ys(1)) - xs(end)*(ys(end)-ys(end-1))
% aa-bb
% aa+bb

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,b]=make_matrix(n,h,val);

A=zeros(n,n);

for i=1:n
    A(i,i)=2;
    if(i>1), A(i,i-1)=-1; end    
    if(i<n), A(i,i+1)=-1; end
end

A=A/h;

b=h*val*[ 0.5; ones(n-2,1); 0.5];

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,b]=apply_bc(A,b,L,R);

n=length(b);
A(1,:)=0; b(1)=L; A(1,1)=1;
A(n,:)=0; b(n)=R; A(n,n)=1;

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,b]=symmetrize_bc(A,b,L,R);

n=length(b);
b=b-A(:,1).*L; A(:,1)=0; A(1,1)=1; b(1)=L;
b=b-A(:,n).*R; A(:,n)=0; A(n,n)=1; b(n)=R;

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all; clear all; clc
function test_dirichlet
close all;

n=15;

Lx=0; Rx=15; vx=10;
Ly=10; Ry=5; vy=1;

[A,b]=make_matrix(n,vx);
[A,b]=apply_bc(A,b,Lx,Rx)
x=A\b;
plot(x); hold all
[As,bs]=symmetrize_bc(A,b,Lx,Rx)
xs=As\bs;
plot(xs,'o')
% As*x-b-(As-A)*x
% bs = b+(As-A)*x

[M,q]=make_matrix(n,vy);
[M,q]=apply_bc(M,q,Ly,Ry)
y=M\q;
plot(y); hold all
[Ms,qs]=symmetrize_bc(M,q,Ly,Ry)
ys=Ms\qs;
plot(ys,'o')

[y' *b  x' *q  y'*(As-A)*x x'*(Ms-M)*y]
[ys'*bs xs'*qs]

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,b]=make_matrix(n,val);

A=zeros(n,n);

for i=1:n
    A(i,i)=2;
    if(i>1), A(i,i-1)=-1; end    
    if(i<n), A(i,i+1)=-1; end
end

b=val*ones(n,1);

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
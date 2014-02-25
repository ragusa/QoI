% close all; clear all; clc
function test_neumann
close all;

n=20;

% forward
Lx= 10; Rx=.10; vx=10;
% adjoint
Ly= -20; Ry= -7; vy= 1;

[A,b]=make_matrix(n,vx);
[A,b]=apply_bc(A,b,Lx,Rx)
x=A\b;
plot(x); hold all
[As,bs]=symmetrize_bc(A,b,Lx,Rx)
xs=As\bs;
plot(xs,'o')

[M,q]=make_matrix(n,vy);
[M,q]=apply_bc(M,q,Ly,Ry)
y=M\q;
plot(y); hold all
[Ms,qs]=symmetrize_bc(M,q,Ly,Ry)
ys=Ms\qs;
plot(ys,'o')

[y' *b  x' *q  y'*(As-A)*x x'*(Ms-M)*y]
[ys'*bs xs'*qs]

aa=ys(1)*(xs(2)-xs(1)) - ys(end)*(xs(end)-xs(end-1))
bb=xs(1)*(ys(2)-ys(1)) - xs(end)*(ys(end)-ys(end-1))
aa-bb
aa+bb
[y-ys x-xs x xs]

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
A(1,1)=1;
A(n,n)=1;

b=val*ones(n,1);

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,b]=apply_bc(A,b,L,R);

n=length(b);
A(1,:)=0; b(1)=L; A(1,1)=1;
b(n)=b(n)-R; 

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,b]=symmetrize_bc(A,b,L,R);

n=length(b);
b=b-A(:,1).*L; A(:,1)=0; A(1,1)=1; b(1)=L;
% b=b-A(:,n).*R; A(:,n)=0; A(n,n)=1; b(n)=R;

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
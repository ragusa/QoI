function [M,D,omega,weights]=gauss_legendre_1d(ndir,n_moments,logi_galerkin)
% create the 1D Gauss-Legendre quadrature
% note: the standard quadrature IS also the Galerkin quadrature...
% no need to distinguish between the two ...

% sanity check
if mod(ndir,2)~=0
    error('The number of directions (%g) in 1D must be even',ndir);
end

% compute the Gauss-Legendre quadrature in 1 D
[nodes,weights] = GLNodeWt(ndir);

% compute associated Legendre polynomials
for L=0:n_moments
    P{L+1}= legendre(L,nodes,'sch');
end

M    = 0.5*P{1}(1,1)*ones(ndir,n_moments);
sphr =     P{1}(1,1)*ones(ndir,n_moments);
for i=2:n_moments
    for j=1:ndir
        M(j,i)    = P{i}(1,j)*((2*(i-1)+1)/2);
        sphr(j,i) = P{i}(1,j);
    end
end

if(logi_galerkin)
    D = inv(M);
else
    w = diag(weights);
    D      = sphr'*w;
end
omega = nodes';
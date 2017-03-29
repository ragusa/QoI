function sum=do_inner_prod(phi,lower,upper)

global npar dat snq
xvec=npar.xf;
sum=0;
for ii=lower+1:2:upper
    phi_avg = (1/2)*(phi(ii)+phi(ii+1));
    dx=xvec(ii+1)-xvec(ii);
    sum=sum+dx*phi_avg;
end
return 
end
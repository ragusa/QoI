function [ kap,sg ] = penalty_coef( method,cdif1,cdif2,d1,d2,alpha )

global npar 
switch method
    case 'MIP'
        kap=max( npar.Ckappa*(cdif1/d1+cdif2/d2) , alpha );
    case {'SIP','NIP'}
        kap=npar.Ckappa*(cdif1/d1+cdif2/d2);
    case 'M4S'
        kap=alpha;
end
switch method
    case {'MIP','SIP'}
        sg=+1;
    case 'NIP'
        sg=-1;
    case 'M4S'
        sg=0;
end

end


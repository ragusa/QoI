function establish_bc_for_VEF( forward )

global dat snq

% source data (incoming) 
if forward
    inc = dat.inc_forward;
else
    inc = dat.inc_adjoint;
end

% bc type:         0=neumann, 1=robin,               2=dirichlet
% C is defined as: kdu/dn=C   // u+k/hcv*du/dn =C // u=C
% u+k/hcv*du/dn = C
% u/4 + D/2*du/dn = Jinc
% u/2 +   D*du/dn = 2*Jinc. we let C=Jinc=0 in DSA

% left incoming in transport
ndir = snq.n_dir;
ii = find (inc(1:ndir/2)>0);
if isempty(ii) % vacuum
    dat.bcVEF.left.type = 1; % Pick your choice: 1=Robin, 2=Dirichlet
    dat.bcVEF.left.C = 0;
else
    dat.bcVEF.left.type = 1; % 1=Robin
    xa=snq.w(1:ndir/2);
    xb=snq.mu(1:ndir/2);
    xc= inc(1:ndir/2);
    xd=snq.w(1:ndir/2)'.*snq.mu(1:ndir/2);
    Jinc = dot (snq.w(1:ndir/2)'.*snq.mu(1:ndir/2), inc(1:ndir/2) );
    dat.bcVEF.left.C = 2*Jinc;
end
ii = find (inc(ndir/2+1:end)>0);
if isempty(ii) % vacuum
    dat.bcVEF.rite.type = 1; % Pick your choice: 1=Robin, 2=Dirichlet
    dat.bcVEF.rite.C = 0;
else
    dat.bcVEF.rite.type = 1; % 1=Robin
    Jinc = dot (snq.w(ndir/2+1:end)'.*snq.mu(ndir/2+1:end), inc(ndir/2+1:end) );
    dat.bcVEF.rite.C = 2*Jinc;
end

% method = 'SIP';

end


function load_input(pb_ID)
% load XS (total and isotropic scattering) for each material
% and load each volumetric external (isotropic) source intensity
% non-isotropic data is not considered here

global dat snq npar

sn=snq.n_dir;
mu=snq.mu;

%%%%%%%%%%%%%%%%%%%%%% select input data

switch pb_ID
    case 1 % problem A from paper
        % number of elements per zone
        nel_zone = [ 40 ];
        % width of each zone
        width_zone = [ 10 ];
        % sigt/sigs per zone
        sigt=[100];
        sigs=[100];
        % volumetric source value, per zone
        qvf=[0.01];
        % incoming flux values
        incf(1:sn)   = 0;
        
    case 2 % problem B from paper
        % number of elements per zone
        nel_zone = [ 10 ];
        % width of each zone
        width_zone = [ 10 ];
        % sigt/sigs per zone
        epsilon=1e-4;
        sigt=[100];
        sigs=[100];
        % volumetric source value, per zone
        qvf=[0];
        % incoming flux values
        incf(1:sn)   = 0;
        incf(sn/2+1:end) = 1;
        
    case 3 % problem C from paper
        % number of elements per zone
        nel_zone = [ 50 ];
        % width of each zone
        width_zone = [ 10 ];
        % sigt/sigs per zone
        epsilon=1e-4;
        sigt=[100];
        sigs=[100];
        % volumetric source value, per zone
        qvf=[0];
        % incoming flux values
        incf(1:sn)   = 0;
        incf(sn)=1/snq.w(sn);
        
    case 4 % problem D from paper
        % number of elements per zone
        nel_zone = [ 10 ];
        % width of each zone
        width_zone = [ 1 ];
        % sigt/sigs per zone
        epsilon=1e-4;
        sigt=[1/epsilon];
        sigs=[1/epsilon-epsilon];
        % volumetric source value, per zone
        qvf=[epsilon];
        % incoming flux values
        incf(1:sn)   = 0;
        
    case 5 % 2-zone problem
        % number of elements per zone
        nel_zone = [ 50 50 ];
        % width of each zone
        width_zone = [ 0.5 0.5 ];
        % sigt/sigs per zone
        sigt=[100 1e4];
        sigs=sigt;
        % volumetric source value, per zone
        qvf=[1e-2 1e-2];
        % incoming flux values
        incf(1:sn)   = 0;
        
    case 6 % 3-zone problem
        % number of elements per zone
        nel_zone = [ 50 50 10]*3;
        % width of each zone
        width_zone = [ 1/3 1/3 1/3 ];
        % sigt/sigs per zone
        sigt=[100 100 100];
        sigs=[50 50 50]*2;
        % volumetric source value, per zone
        qvf=[1 2 1];
        % incoming flux values
        incf(1:sn)   = 0;
        
    case 7 % Reed (5-zone problem)
        % number of elements per zone
        nel_zone = [ 20 10 20 10 20]*4;
        % width of each zone
        width_zone = [ 2 1 2 1 2 ];
        % sigt/sigs per zone
        sigt=[50 5 1e-8 1 1 ];
        sigs=[ 0 0    0 0.9 0.9];
        % volumetric source value, per zone
        qvf=[50 0 0 1 0];
        % incoming flux values
        incf(1:sn)   = 0;
        % volumetric source value, per zone
        qva=[ 0 0 0 1 0];
        % incoming flux values
        inca(1:sn)   = 0;
        
    case 8 % Adams (PIM)
        % number of elements per zone
        nel_zone = [ 10 ];
        % width of each zone
        width_zone = [ 500 ];
        % sigt/sigs per zone
        sigt=[1];
        sigs=[1];
        % volumetric source value, per zone
        qvf=[0];
        % incoming flux values
        incf(1:sn)   = 0;
        incf(sn/2+1)=1/snq.w(sn/2+1);
        
    case 9 % Adams (PIM) ref
        % number of elements per zone
        nel_zone = [ 1000 1000];
        % width of each zone
        width_zone = [ 5 495 ];
        % sigt/sigs per zone
        sigt=[1 1];
        sigs=[1 1];
        % volumetric source value, per zone
        qvf=[0 0];
        % incoming flux values
        incf(1:sn)   = 0;
        incf(sn/2+1)=1/snq.w(sn/2+1);
        
    case 10 % Test thick diffusive limit
        % number of elements per zone
        nel_zone = [ 100 ];
        % width of each zone
        width_zone = [ 500 ];
        % sigt/sigs per zone
        sigt=[1];
        sigs=[0.9999];
        % volumetric source value, per zone
        qvf=[0];
        % incoming flux values
        incf(1:sn)   = 0;
        
    case 11 % problem A from paper - modified for IP-testing
        % number of elements per zone
        nel_zone = [ 10 10 10 ]*10;
        % width of each zone
        width_zone = [ 10 10 20];
        % sigt/sigs per zone
        sigt=[1 2 3];
        sigs=[0.5 1. 2. ];
        % volumetric source value, per zone
        qvf=[1 0 3];
        % incoming flux values
        incf(1:sn)   = 0;
        
    case 12 %
        % number of elements per zone
        nel_zone = [ 10 10 10 ]*4;
        % width of each zone
        width_zone = [ 2 2 2 ];
        % sigt/sigs per zone
        sigt=[1 1e-8 1];
        sigs=[0.3 0 0.3];
        % volumetric source value, per zone
        qvf=[1 0 0];
        % incoming flux values
        incf(1:sn) = 0;
        % volumetric source value, per zone
        qva=[0 0 1];
        % incoming adj flux values
        inca(1:sn) = 0;
        
    case 13 %
        % number of elements per zone
        nel_zone = [ 100  20 10 20 100 ];
        % width of each zone
        width_zone = [ 20 4 2 4 20 ]/5;
        % sigt/sigs per zone
        sigt=[1 1 1 1 1];
        sigs=[0 0 0 0 0];
        % volumetric source value, per zone
        qvf=[1 1 1 1 1];
        % incoming flux values
        incf(1:sn) = 0;  incf(1:sn/2)=0;
        % volumetric source value, per zone
        qva=[1 1 1 1 1];
        % incoming adj flux values
        inca(1:sn) = 0;
        
    case 14 % flat solution !!!
        % number of elements per zone
        nel_zone = [ 100  20 10 20 100 ]*1;
        % width of each zone
        width_zone = [ 20 4 2 4 20 ]/5;
        % sigt/sigs per zone
        sigt=[1 1 1 1 1];
        sigs=[1 1 1 1 1]*0;
        % volumetric source value, per zone
        qvf=[1 1 1 1 1];
        % incoming flux values
        incf(1:sn) = 1/2;  %incf(1:sn/2)=0.6; incf(sn/2+1:sn)=0.4;
        % volumetric source value, per zone
        qva=[1 1 1 1 1];
        % incoming adj flux values
        inca(1:sn) = 1; inca(1:sn/2)=.1;

    case 99 % Ian's 5 regions
        % number of elements per zone
        nel_zone = [ 100  20 10 20 100 ]*1;
        % width of each zone
        width_zone = [ 20 4 2 4 20 ]/5;
        % sigt/sigs per zone
        sigt=[1 1 1 1 1];
        sigs=[1 1 1 1 1]*0;
        % volumetric source value, per zone
        qvf=[1 1 1 1 1];
        % incoming flux values
        incf(1:sn) = 1/2;  %incf(1:sn/2)=0.6; incf(sn/2+1:sn)=0.4;
        % volumetric source value, per zone
        qva=[1 1 1 1 1];
        % incoming adj flux values
        inca(1:sn) = 1; inca(1:sn/2)=.1;

    otherwise
        error('problem ID %g is unknown',pb_ID);
        
end

% sanity checks
if length(nel_zone) ~= length(width_zone)
    error('length(nel_zone) ~= length(width_zone)');
end
if length(nel_zone) ~= length(qvf)
    error('length(nel_zone) ~= length(qv)');
end
if length(nel_zone) ~= length(sigt)
    error('length(nel_zone) ~= length(sigt)');
end
if length(nel_zone) ~= length(sigs)
    error('length(nel_zone) ~= length(sigs)');
end

% another simple check
aux = sigs./sigt;
i=find(aux>1);
if(~isempty(i))
    error('sigs cannot be > than sigt \nA problem occured with zones %i in file %s',i,mfilename);
end


% save data
dat.pb_ID=pb_ID;
dat.sigt = sigt;
dat.sigs = sigs;
dat.qv_forward =qvf;
dat.inc_forward = incf;
if pb_ID==12 || pb_ID==7 || pb_ID==13 || pb_ID==14
    dat.qv_adjoint =qva;
    dat.inc_adjoint = inca;
else
    error('fix load input');
    dat.qv_adjoint =qvf;
    dat.inc_adjoint = incf;
end

%%%%%%%%%%%%%%%%%%%%%%%%% prepare data for computation

% finish computing some XS data
dat.siga = dat.sigt-dat.sigs;
dat.cdif = 1./(3*dat.sigt);

% total number of elements
nel = sum(nel_zone);
% number of spatial dofs (LD: 2 spatial unknowns/cell)
npar.porder=1;
ndofs=nel*(npar.porder+1);

% delta x and iel2zon
tmp_x = [ 0 ];
for i=1:length(nel_zone)
    tmp_x = [ tmp_x sum(width_zone(1:i)) ];
end
x=[];
iel2zon=[];
for z=1:length(nel_zone)
    x_zone = linspace(tmp_x(z),tmp_x(z+1),nel_zone(z)+1);
    if isempty(x)
        x = x_zone;
        iel2zon=z*ones(nel_zone(z),1);
    else
        x = [x x_zone(2:end)];
        iel2zon =[ iel2zon; z*ones(nel_zone(z),1)];
    end
end
dx=diff(x);

% save data
npar.ndofs=ndofs;
npar.nel=nel;
npar.x=x';
npar.dx=dx';
npar.iel2zon=iel2zon;

% build connectivity array (used in diffusion)
gn=zeros(npar.nel,npar.porder+1);
gn(1,:)=linspace(1,npar.porder+1,npar.porder+1);
for iel=2:npar.nel
    gn(iel,:)=gn(iel-1,:) + npar.porder+1 ;
end
npar.gn=gn; clear gn;

% create x values for plotting FEM solution
npar.xf=zeros(npar.nel*(npar.porder+1),1);
for iel=1:npar.nel
    ibeg = (iel-1)*(npar.porder+1)+1;
    iend = (iel  )*(npar.porder+1)  ;
    npar.xf(ibeg:iend)=linspace(npar.x(iel),npar.x(iel+1),npar.porder+1) ;
end


function load_reed(input_vars)
% load XS (total and isotropic scattering) for each material
% and load each volumetric external (isotropic) source intensity
% non-isotropic data is not considered here

global dat snq npar

sn=snq.n_dir;
mu=snq.mu;

r1q=input_vars(1); %50
r4q=input_vars(2); %1
r1a=input_vars(3); %50
r2a=input_vars(4); %5
r4a=input_vars(5); %0.1
r5a=input_vars(6); %0.1
r4s=input_vars(7); %0.9
r5s=input_vars(8); %0.9


%%%%%%%%%%%%%%%%%%%%%% select input data
% number of elements per zone
nel_zone = [ 10 5 10 5 10]*10;
% width of each zone
width_zone = [ 2 1 2 1 2 ];
% sigt/sigs per zone
siga=[r1a r2a 1e-8 r4a r5a ];
sigs=[ 0 0 0 r4s r5s];
sigt=siga+sigs;
% volumetric source value, per zone
qvf=[r1q 0 0 r4q 0];
% incoming flux values
incf(1:sn)   = 0;
% volumetric response value, per zone
qva=[ 0 0 0 0 1];
% incoming flux values
inca(1:sn)   = 0;
%Regions to be perturbed. Use value of 1 to specify
dat.sigaPertRegion=[1 1 0 0 0];
dat.sigsPertRegion=[0 0 0 1 1];
dat.sourcePertRegion=[1 0 0 1 0];
dat.incPertRegion(1:sn) = 0;
dat.incPertRegion((sn/2)+1:sn) = 0;
dat.name='Reed System';

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
    %error('sigs cannot be > than sigt \nA problem occured with zones %i in file %s',i,mfilename);
end


% save data
dat.sigt = sigt;
dat.sigs = sigs;
dat.siga = siga;
dat.qv_forward =qvf;
dat.inc_forward = incf;
dat.qv_adjoint =qva;
dat.inc_adjoint = inca;

if ~isfield(dat,'name')
    dat.name=['NO NAME'];
end

%%%%%%%%%%%%%%%%%%%%%%%%% prepare data for computation

% finish computing some XS data
%dat.siga = dat.sigt-dat.sigs;
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
end


function [x,pb_data]=load_input(n_ref)

% XS per material (4 medium total)
XS.tot=[50.; 5.; 0.; 1. ];
XS.sca=[ 0.; 0.; 0.; 0.9];
% source strength
s_ext=[50. 1.];
% geometry input: medium ID (all cell width=1cm)
matID=[1 1 2 3 3 4 4 4];
% geometry input: src ID
srcID=[1 1 0 0 0 2 0 0];

if (n_ref<=0), error('n_ref given is %i but it cannot be <=0',n_ref); end
n_mesh = length(matID) * n_ref;
x=linspace(0,length(matID),n_mesh+1);

ID.mat=kron( matID, ones(1,n_ref) );
ID.src=kron( srcID, ones(1,n_ref) );

pb_data.XS =XS;
pb_data.ID =ID;
pb_data.src=s_ext;

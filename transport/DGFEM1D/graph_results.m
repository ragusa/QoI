function graph_results
%Specify Input File Path
file='DGFEM1D_prob34_0605171121.csv';
filename=fullfile('C:\Users\Ian\checkout','output',file);
%Load into matrix M and define data values
M=csvread(filename);
probID=M(:,1);
qoi_sn_f=M(:,2);
qoi_VEF_f=M(:,3);
qoi_sn_a=M(:,4);
qoi_VEF_a=M(:,5);
sigt_pert=M(:,6);
sigs_pert=M(:,7);
source_pert=M(:,8);
inc_pert=M(:,9);
dqoi_sn_f=M(:,10);
dqoi_VEF_f=M(:,11);
dqoi_sn_a=M(:,12);
dqoi_VEF_a=M(:,13);
E_diff=M(:,14);

figure(1)
clf
hold on
xlabel('\sigma_t % change')
ylabel('QoI % response')
plot(sigt_pert,dqoi_sn_f./qoi_sn_f,'-+r')
plot(sigt_pert,dqoi_VEF_f./qoi_sn_f,'--+g')
plot(sigt_pert,dqoi_sn_a./qoi_sn_f,'-+m')
plot(sigt_pert,dqoi_VEF_a./qoi_sn_f,'--+b')
legend('sn forward','VEF forward','sn adjoint','VEF adjoint')
end
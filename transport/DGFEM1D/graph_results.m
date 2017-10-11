function graph_results
%Specify Input File Path
file='DGFEM1D_prob57_1011171553.csv';
filename=fullfile('C:\Users\Ian\checkout','output',file);
%filename=fullfile('E:\Storage\git_checkout','output',file);
%Load into matrix M and define data values
M=csvread(filename);
probID=M(:,1);
qoi_sn_f=M(:,2);
qoi_VEF_f=M(:,3);
qoi_sn_a=M(:,4);
qoi_VEF_a=M(:,5);
siga_pert=M(:,6);
sigs_pert=M(:,7);
source_pert=M(:,8);
inc_pert=M(:,9);
dqoi_sn_f=M(:,10);
dqoi_VEF_f=M(:,11);
dqoi_sn_a=M(:,12);
dqoi_VEF_a=M(:,13);
E_diff=M(:,14);
deltaEterm=M(:,15);

xvalues=sigs_pert;

close all
figure(1)
clf
hold on
title(file,'interpreter','none')
xlabel('\sigma_s % change')
ylabel('QoI % response')
plot(xvalues,dqoi_sn_f./qoi_sn_f,'-+r')
plot(xvalues,dqoi_VEF_f./qoi_sn_f,'--+g')
plot(xvalues,dqoi_sn_a./qoi_sn_f,'-+m')
plot(xvalues,dqoi_VEF_a./qoi_sn_f,'--+b')
legend({'sn forward','VEF forward','sn adjoint','VEF adjoint'},'Position',[0.5 0.80 0.01 0.01])

figure(2)
plot(xvalues,deltaEterm,'--+r')
legend('delta E term')
end
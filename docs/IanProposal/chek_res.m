clear all; close all; clc;

% uniform system with uniform perturbation. qoi in the middle so pretty
% much the infinite solution q/(sigt-sigs).

%%% perturbation in sigt
sigt=2;
sigs=1;

pert=linspace(-0.1,0.1);
sigt_p = sigt + sigt*pert;

qoi = 1/(sigt-sigs);
qoi_p = 1./(sigt_p-sigs);

figure(1);
plot(pert,(qoi_p-qoi)/qoi);

%%% perturbation in sigs
sigt=2;
sigs=1;

pert=linspace(-0.1,0.1);
sigs_p = sigs + sigs*pert;

qoi = 1/(sigt-sigs);
qoi_p = 1./(sigt-sigs_p);

figure(2);
plot(pert,(qoi_p-qoi)/qoi);

%%% perturbation in qvol
% no need to plot: must be a straight line always
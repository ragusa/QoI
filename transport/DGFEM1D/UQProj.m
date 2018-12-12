function UQProj()
% Linear Discontinous FEM code for Sn transport in 1D
% Jean Ragusa, Ian Halvic
close all;
clc; clear variables; clear global;
%UQ Lab Setup
clearvars
uqlab
rng(100, 'twister') % also set the random seed to a fixed value for repeatable results
global npar dat snq IO_opts results
%Start
input_vars=[50 1 50 5 0.1 0.1 0.9 0.9];
%UQLab Inputs
Model.mFile = 'reed_wrapper';
myModel = uq_createModel(Model);
%Region 1 Source
Input.Marginals(1).Name = 'R1_q';
Input.Marginals(1).Type = 'Gaussian';
Input.Marginals(1).Moments = [50 5]; 
%Region 4 Source
Input.Marginals(2).Name = 'R4_q';
Input.Marginals(2).Type = 'Gaussian';
Input.Marginals(2).Moments = [1 0.1]; 
%Region 1 Sig_a
Input.Marginals(3).Name = 'R1_sig_a';
Input.Marginals(3).Type = 'Gaussian';
Input.Marginals(3).Moments = [50 5]; % (cm^-1)
%Region 2 Sig_a
Input.Marginals(4).Name = 'R2_sig_a';
Input.Marginals(4).Type = 'Gaussian';
Input.Marginals(4).Moments = [5 0.5]; % (cm^-1)
%Region 4 Sig_a
Input.Marginals(5).Name = 'R4_sig_a';
Input.Marginals(5).Type = 'Gaussian';
Input.Marginals(5).Moments = [0.1 0.01]; % (cm^-1)
%Region 5 Sig_a
Input.Marginals(6).Name = 'R5_sig_a';
Input.Marginals(6).Type = 'Gaussian';
Input.Marginals(6).Moments = [0.1 0.01]; % (cm^-1)
%Region 4 Sig_s
Input.Marginals(7).Name = 'R4_sig_s';
Input.Marginals(7).Type = 'Gaussian';
Input.Marginals(7).Moments = [0.9 0.09]; % (cm^-1)
%Region 5 Sig_s
Input.Marginals(8).Name = 'R5_sig_s';
Input.Marginals(8).Type = 'Gaussian';
Input.Marginals(8).Moments = [0.9 0.09]; % (cm^-1)
myInput = uq_createInput(Input);
PCESobol.Type = 'Sensitivity';
PCESobol.Method = 'Sobol';
PCESobol.Sobol.Order = 2;
% PCESobolAnalysis = uq_createAnalysis(PCESobol);
% SobolResults = SobolAnalysis.Results;
%uq_print(SobolAnalysis)

% 4.2 Create a PCE of the full model
% select the metamodel tool in UQLab
PCEOpts.Type = 'Metamodel';
% choose the polynomial chaos expansions module
PCEOpts.MetaType = 'PCE';
% Specify the model that will be sampled to create the experimental design
% for PCE
PCEOpts.FullModel = myModel;
% Specify the maximum polynomial degree (default: sparse PCE expansion)
PCEOpts.Degree = 5;
%  Specify the size of the experimental design (total cost of the metamodel)
PCEOpts.ExpDesign.NSamples = 50;









setup_solver()

%  Calculate and store the PCE in UQLab
myPCE = uq_createModel(PCEOpts);
PCESobolAnalysis = uq_createAnalysis(PCESobol);
PCESobolResults = PCESobolAnalysis.Results;
uq_print(PCESobolAnalysis,1)
uq_print(PCESobolAnalysis,2)
out=PCESobolAnalysis.Results.Total;
names=PCESobolAnalysis.Results.VariableNames;  
fprintf('V & SN5 & VEF5 & SN3 & VEF3 & dE \\\\ \n');
for ii=1:length(out)
    fprintf('%s & %g & %g & %g & %g & %g \\\\ \n',names{ii},out(ii,1),out(ii,2),out(ii,3),out(ii,4),out(ii,5));
end
end

function setup_solver()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global npar dat snq IO_opts results
% set variable once for all
dat.forward_flux = true;
dat.adjoint_flux = ~dat.forward_flux;
do_transport_adjoint=false;
IO_opts.show_dqoi_pre_bc=false;
IO_opts.show_dqoi_from_bc=false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nominal_input=[50 1 50 5 0.1 0.1 0.9 0.9];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select angular approx (must be an even number)
sn=8;
n_moments=1; logi_galerkin=false;
[M,D,omega,weights] = gauss_legendre_1d(sn,n_moments,logi_galerkin);
snq.D = D; snq.M = M; snq.n_dir = sn;
snq.mu = omega; snq.w = weights;
%snq.mu = [-1,1];
%snq.mu;
%snq.w;
% sum of the weights
snq.sw = sum(snq.w);
IO_opts.console_io = false;
dat.do_dsa = true;
dat.NTD_sweep = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
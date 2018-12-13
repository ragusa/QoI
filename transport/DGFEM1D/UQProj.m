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
dat.nominal_input=[50 1 50 5 0.1 0.1 0.9 0.9];
dat.use_reduced_input=false;
do_nsamp_plot=true;
%UQLab Inputs
Model.mFile = 'reed_wrapper';
%Model.mHandle = @qoi_model; 
qoiModel = uq_createModel(Model);
if dat.use_reduced_input
    %Region 2 Sig_a
    Input.Marginals(1).Name = 'R2_sig_a';
    Input.Marginals(1).Type = 'Gaussian';
    Input.Marginals(1).Moments = [5 0.5]; % (cm^-1)
    %Region 4 Sig_s
    Input.Marginals(2).Name = 'R4_sig_s';
    Input.Marginals(2).Type = 'Gaussian';
    Input.Marginals(2).Moments = [0.9 0.09]; % (cm^-1)
else
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
end
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
PCEOpts.FullModel = qoiModel;
% Specify the maximum polynomial degree (default: sparse PCE expansion)
PCEOpts.Degree = 2;
%  Specify the size of the experimental design (total cost of the metamodel)
PCEOpts.ExpDesign.NSamples = 10;

setup_solver()

%  Calculate and store the PCE in UQLab
myPCE = uq_createModel(PCEOpts);
PCESobolAnalysis = uq_createAnalysis(PCESobol);
PCESobolResults = PCESobolAnalysis.Results;
%uq_print(PCESobolAnalysis,1)
%uq_print(PCESobolAnalysis,2)
out=PCESobolAnalysis.Results.Total;
%out=PCESobolAnalysis.Results.FirstOrder;
% fprintf('V & SN5 & VEF5 & SN3 & VEF3 & dE \\\\ \n');
% for ii=1:size(out,1)
%     fprintf('%s & %g & %g & %g & %g & %g \\\\ \n',names{ii},out(ii,1),out(ii,2),out(ii,3),out(ii,4),out(ii,5));
% end
fprintf('  & \\xi_1 & \\xi_2 & \\xi_3 & \\xi_4 & \\xi_5 & \\xi_6 & \\xi_7 & \\xi_8 \\\\ \n');
fprintf('S')
for ii=1:size(out,1)
    fprintf('& %g ',out(ii,5));
end
fprintf('\\\\ \n')

ref_vec=[];
model_vec=[];
pert_vec=[];
%%%%%%%%%%%%%Plot some output
for ii=1:0.01:1.1
    test_input = [50 1 50 5 0.1 0.1 ii*0.9 ii*0.9];
    load_reed(test_input)
    [phi_ref,E_ref,Ebd_ref,psi_ref]=solveForwardSn;
    qoi3_sn = compute_qoi(dat.forward_flux,phi_ref,snq.n_dir,psi_ref,psi_ref);
    model_out= uq_evalModel(myPCE,test_input);
    qoi3_model=model_out(4);
    ref_vec=[ref_vec qoi3_sn];
    model_vec=[model_vec qoi3_model];
    pert_vec=[pert_vec ii];
end

% figure(3)
% hold on
% plot(pert_vec,ref_vec,'-+b')
% plot(pert_vec,model_vec,'-+r')

sn=snq.n_dir;
%Transport solve of unperturbed problem: Forward and Sn/VET adjoint
load_reed(dat.nominal_input)
[phi_ref,E_unpert,Ebd_ref,psi_ref]=solveForwardSn;
[phia,Ea,Ebda,psia]=solveAdjointSn;
phiVEFa=solve_VEF_math_adjoint(dat.adjoint_flux,E_unpert,Ebd_ref);
qoi_ref = compute_qoi(dat.forward_flux,phi_ref,snq.n_dir,psi_ref,psi_ref);
E_unpert=reshape(E_unpert,800,1);

%Transport solve of perturbed problem for reference
test_input = [50 1 50 5 0.1 0.1 0.99 0.99];
load_reed(test_input)
[phi_ref,E_pert,Ebd_ref,psi_ref]=solveForwardSn;
E_pert=reshape(E_pert,800,1);

%Model the perturbed eddington
if dat.use_reduced_input
    test_input = [5 0.99];
end
model_out= uq_evalModel(myPCE,test_input);
E_model=model_out(6:805)';


%Start plotting output
figure(1)
do_plot(phi_ref,'$$\phi_{nom}$$',1,dat.forward_flux)

fprintf('Nominal: %d \n',qoi_ref);

% figure(4)
% hold on
% plot(E_pert-E_unpert,'-b')
% plot(E_model-E_unpert,'--r')

figure(5)
hold on
plot(E_pert-E_unpert,'-b')
legstr={'Ref'};
legstr2={'Ref'};

if dat.use_reduced_input
    model_input = [5 0.99];
else
    model_input = [50 1 50 5 0.1 0.1 0.99 0.99];
end


%Plot reference delta QOI solution
ref_qoi_vec=[];
ref_qoi_adj_vec=[];
ref_qoi_adjVET_vec=[];
pert_vec=[];
for ii=-0.1:0.01:0.1
    test_input = [50 1 50 5 0.1 0.1 (1+ii)*0.9 (1+ii)*0.9];
    load_reed(test_input)
    [phi_pert,E_pert,Ebd_pert,psi_pert]=solveForwardSn;
    qoi3_sn = compute_qoi(dat.forward_flux,phi_pert,snq.n_dir,psi_pert,psi_pert);
    ref_qoi_vec=[ref_qoi_vec qoi3_sn-qoi_ref];
    pert_vec=[pert_vec ii];
end

figure(7)
hold on
plot(100*pert_vec,100*ref_qoi_vec/qoi_ref,'-+b')
%plot(pert_vec,model_vec,'-+r')

load_reed(dat.nominal_input)

for ii=7:-1:3
    %Plot the dE approximation
    PCEOpts.ExpDesign.NSamples = 2^ii;
    myPCE = uq_createModel(PCEOpts);
    model_out= uq_evalModel(myPCE,model_input);
    E_model=model_out(6:805)';
    figure(5)
    hold on
    plot(E_model-E_unpert,'-')
    legstr{end+1}=num2str(2^ii);
    legstr2{end+1}=num2str(2^ii);
    %Compute sensitivity using Adj and dE from PC
    delta_E_correction=[];
    for jj=-0.1:0.01:0.1
        if dat.use_reduced_input
           model_in=[5 (1+jj)*0.9];
        else
           model_in=[50 1 50 5 0.1 0.1 (1+jj)*0.9 (1+jj)*0.9];
        end
        model_out= uq_evalModel(myPCE,model_in);
        E_model=model_out(6:805)';
        %Load perturbation for adjoint
        dat.sigsPert = dat.sigsPertRegion.*dat.sigs*jj;
        dat.sigtPert = dat.sigsPert;
        dat.sigaPert = dat.sigaPertRegion.*0.0;
        dat.sourcePert =dat.sourcePertRegion.*0.0;
        dat.psiIncPert = dat.incPertRegion.*0.0;
        %Compute dE Correction
        [deltaE_term,deltaB_L,deltaB_R,~]=compute_deltaE_QoI_term(phiVEFa,phi_ref,phi_ref,reshape(E_unpert,2,400),reshape(E_model,2,400),Ebd_ref,Ebd_ref);
        delta_qoi_VEF_a_Eint=deltaE_term+deltaB_L+deltaB_R;
        delta_E_correction=[delta_E_correction delta_qoi_VEF_a_Eint]
        if ii==7
            dqoi3_sn = compute_perturbed_qoi_Sn(dat.adjoint_flux,phia,phi_ref,psi_ref,psia,sn);
            [dqoi3_VET, ~] = compute_perturbed_qoi_VEF(dat.adjoint_flux,phiVEFa,phi_ref,reshape(E_model,2,400));
            ref_qoi_adj_vec=[ref_qoi_adj_vec dqoi3_sn]
            ref_qoi_adjVET_vec=[ref_qoi_adjVET_vec dqoi3_VET]
        end
    end
    figure(7)
    hold on
    plot(100*pert_vec,100*(ref_qoi_adjVET_vec+delta_E_correction)/qoi_ref,'-+')
end
figure(5)
legend(legstr,'Location','best');
xlabel('x')
ylabel('\delta E')

legstr2{end+1}='Trans_{adj}';
legstr2{end+1}='VET_{adj}';
figure(7)
hold on
plot(100*pert_vec,100*ref_qoi_adj_vec/qoi_ref,'-+')
plot(100*pert_vec,100*ref_qoi_adjVET_vec/qoi_ref,'-+')
legend(legstr2,'Location','best');
xlabel('% \delta \sigma_s')
ylabel('% \delta QOI')


fprintf('Nominal: %d \n',qoi_ref);

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

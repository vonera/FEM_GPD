%% This script is used to test FEM-GPD on a test-case
clearvars
close all

%% Matlab configuration:
% add the framework to Matlab-path 
addpath(genpath('../FEM_Framework'));

%% If gurobi is used uncomment lines 11-14 and adjust line 12 
% save path and setup gurobi solver, then return to the current path
% path = pwd;
% cd /Library/gurobi751/mac64/matlab # PATH_TO_GUROBI/mac64/matlab
% gurobi_setup
% cd(path)

%% set rng and save its configuration
rng('shuffle');
my_rng = rng;

%% Generate Data
[Pcell, NLL_gp, Gamma, input] = generate_data();


%% FEM-GPD data configuration
% N is the number of locations
N = 1;
% if not used initialize random tmp_dist = spatial regularization matrix; 
% if spatial regularization required, set tmp_dist to a meaningful distance
% matrix between locations
tmp_dist = rand(N);
tmp_dist = 0.5*(tmp_dist + tmp_dist');

% xt is the input, i.e., threshold excesses for a single location
% ut is the set of covariates for each location
xt = cell(1,N);
ut = cell(1,N);

for j = 1:N
    xt{j} = input(1).xt;
    ut{j} = input(1).ut;
end

%% Start FEM-GPD fro diff lasso configurations ans save results
% choose the best result wrt to AICc
for lasso_flag = 0%[0 0.0001, 0.001, 0.01, 0.1, 1]
      
    cfg = adjustCFG([], tmp_dist, lasso_flag);
    cfg = Core.Construction.FEMConfig(cfg.modelFit, cfg.gamma, cfg.model);
    
    iHandler = FEM.createInputHandler(cfg, xt, ut);
    
    tic
    rHandler = FEM.fit(iHandler);
    time_needed = toc;
    
    rHandler.computeIC(iHandler);
    rHandler.printAllInfo();
    
    resCell = rHandler.getIterable();
    optRes = rHandler.getBest('AICc');
    
    save(['result_testcase_lasso',num2str(lasso_flag),'_AICc',num2str(optRes.ic.AICc),'.mat'])
end
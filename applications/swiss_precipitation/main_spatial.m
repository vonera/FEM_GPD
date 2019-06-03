%function main_spatial()

% this script is used to test the spatial FEM-BV application to real data
clearvars
close all

addpath(genpath('../../FEM_Framework'));
% path = pwd;
% cd /Library/gurobi604/mac64/matlab
% gurobi_setup
% cd(path)

warning 'off'
rng('shuffle');
my_rng = rng;

N = 17;
S = load('data/daily_accumulated_precipitation_0.98_delayedENSO_diffnorm.mat');
tmp_dist = rand(N);
tmp_dist = 0.5*(tmp_dist + tmp_dist');

xt = cell(1,N);
ut = cell(1,N);

for j = 2:N+1
    % threshold excesses follow gpd, thus Xt - threshold!
    xt{j-1} = S.daily_accumulated_precipitation.data{j,7} - S.daily_accumulated_precipitation.data{j,5};
    ut{j-1} = S.daily_accumulated_precipitation.data{j,8};
end
clear S

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
    
    save(['ksi_const_sinoNaoEnso_potSwiss_QL0.98_lasso',num2str(lasso_flag),'_AICc',num2str(optRes.ic.AICc),'.mat'])
    
end
function cfg = adjustCFG(posnSwitch, distMatrix, lasso_flag)
%ADJUSTCFG is used to adjust the config for the spatial FEM,
% including the model and the corresponding parameters

%%%%%%%%%%%%%%%%%%%%%%%
%choose the model
%cfg = FEM.getdefaultCfg('MCMC','parforSPATIAL_BV_Gurobi');

% attention using matlab will be very slow!
cfg = FEM.getdefaultCfg('MCMC','SPATIAL_BV_Mat');
cfg.model.subtype = 'GPDcpp';

%%%%%%%%%%%%%%%%%%%%%%%%%
% expected same number of combinations in covariates for ksi, si
cfg.model.ut_combis{1} = {[3 5 7 8 11]}; % cell of arrays: {[1 2 3] [2 3 4]} ...
cfg.model.ut_combis{2} = {[3 5 7 8 11]}; 

%%%%%%%%%%%%%%%
% FEM options
cfg.modelFit.subspaceEps = 1e-03;
cfg.modelFit.annealing = 100;
cfg.modelFit.subspaceItr = 1000;
cfg.modelFit.Klist = [1 2];
cfg.modelFit.Ctlist = 10:10:20;% 10:10:100;
cfg.modelFit.Cnlist = inf; % please note inf for BV indicates: no spatial regularisation

%%%%%%%%%%%%%%%%%
% MCMC adjustment
cfg.model.lassoFLAG = lasso_flag; 
cfg.model.ridgeFLAG = 0;
cfg.model.paramTHRESHOLD = 1e-04;

%%%%%%%%%%%%%%%%%
% Default, but can be changed
cfg.model.proposMaxTimes = 200;
cfg.model.samples = 10000;
cfg.model.gaussFLAG = true;
cfg.model.init_noise = 0.0006;
cfg.model.beta_factor = 2.5;
cfg.model.init_beta = 1;
if cfg.model.gaussFLAG
    cfg.model.noise_less = ( 1 / cfg.model.init_noise ) ^ ( 50 / cfg.model.samples);
    cfg.model.noise_more = ( 0.000001 / cfg.model.init_noise )^( 50 / cfg.model.samples );
else
    cfg.model.noise_less = 0.8;
    cfg.model.noise_more = 1.9; 
end

%%%%%%%%%
% gamma options for SPATIAL
cfg.gamma.num_fe   = 600;
cfg.gamma.discrete = [];
cfg.gamma.distMatrix = distMatrix;
cfg.gamma.positionSwitch = [];
 

%%%%%%%%%%%
%options for 'SPATIAL_BV_Gurobi'
cfg.gamma.discrete = 'I'; %'B'-binary; 'I'-integer; 'C'-continous
cfg.gamma.diagnostics = 0; 

%%%%%%%%
% options for SPATIAL_BV_Mat
% cfg.gamma.options = optimset('Display','off','LargeScale', 'off', 'Simplex', 'on', 'TolCon', 2.2204e-5);


end

## FEM-BV-GPD
We present here the FEM-BV-GPD framework for a data-driven spatio-temporal regression based clustering of threshold excesses.



### Functionality: 
Given a series of threshold excesses and covariates for multiple locations the framework provides an optimal sparse solution for the spatio-temporal clustering problem by incorporating information theory based model selection criteria and L1-regularized regression. The framework goes beyond a priori probabilistic and deterministic assumptions typical for standard approaches like Hidden Markov Models (HMM) and Gaussian Mixture Models (GMM). 


The nonstationarity and the non-homogeneity of the data are resolved by describing the underlying behavior by a set of K>=1 locally stationary models and a spatial nonstationary switching process Γ(s, t), where K > 1 indicates the existence of systematically missing covariates.


The realization of framework allows to consider simultaneously different configurations of the problem referring to hyper-parameter tuning of:
 * the number of clusters,
 * number of switches,
 * different configuration of covariates
 * different L1-regularizations 


The resulting spatial FEM-BV-GPD can be employed as a robust exploratory sparse regression analysis tool for spatio-temporal extremes, that enables to identify the significant resolved covariates and to account for the influence from the systematically missing ones.

#### Spatial dependence structure 
The resulting switching process Γ(s, t) provides a posteriori a pragmatic description of the underlying spatial dependence structure by grouping together all locations that exhibit similar behavior in Γ(s, t), e.g., by estimating the Event Synchronization (ES) measure matrix [1] which allows to identify contiguous regions that exhibit similar spatio-temporal behavior. The estimation of the ES-matrix is included in the framework.

[1] N. Malik, B. Bookhagen, N. Marwan, and J. Kurths. Analysis of spatial and temporal extreme monsoonal rainfall over south asia using complex networks.

### Installation

#### Prerequisites
 * GSL library 
 * Matlab and Matlab mex compiler
 * Matlab mex compiler, deployed MCMC solver is implemented in C++, you need to mex compile it in order to use it. The HowTo can be found in FEM_Framework/+Models/+MCMC/+GPD/+GPDcpp/
 * Gurobi (optional, highly recommended for fast runs, free for academics)


#### Installation/configuration:
1. Download FEM_GPD
2. Mex-Compile the C++ extensions: 
	* open Matlab
	* go to FEM_GPD/+Models/+MCMC/+GPD/+GPDcpp/
	* change PATH_TO_GSL in 
		* mex MEX_GPDresid.cpp GSL_Helper.cpp -I/PATH_TO_GSL/include -L/Users/loka/PATH_TO_GSL/lib -lgsl
		* mex MEX_GPDestimate.cpp GSL_Helper.cpp MCMCMetGPD.cpp MCMCMet.cpp -I/PATH_TO_GSL/include -L/PATH_TO_GSL/lib -lgsl
		* MEX_GPDcomputeAcf.cpp  GSL_Helper.cpp MCMCMetGPD.cpp MCMCMet.cpp -I/PATH_TO_GSL/include -L/PATH_TO_GSL/lib -lgsl
	 
	and execute these three lines after each other in the MATLAB command line. If y can not compile the mex files (no gsl, no mex-complier) you can use the MATLAB GPD version (very slow!) by using 
	 	cfg.model.subtype = 'GPDmatlab';
	* Go to FEM_GPD/demo and run main_demo.m
	* In case you use Gurobi, uncomment the following lines in main_demo and adjust the PATH_TO_GUROBI (check gurobi instructions)
	```matlab
		 path = pwd;
		 cd PATH_TO_GUROBI/mac64/matlab
		 gurobi_setup
		 cd(path)
	```
	* Run main_demo.m


### Usage
In Brief: Follow main_spatial.m in FEM_GPD/applications/swiss_precipitation.
In more detail: 
* Create a folder my_application under FEM_GPD/applications
* copy adjustCFG.m and main_spatial.m from FEM_GPD/applications/swiss_precipitation to FEM_GPD/applications/my_application
* adjust adjustCFG.m (refer to wiki for more details on the configuration file)
* adjust main_spatial.m wrt to data loading
* run main_spatial.m
	

#### Data loading
The framework expects  a cell of size N for the excesses and for the covariates.
Each xt cell contains a (1,Tj)-dimensional series; and each ut cell contains a (dimu, Tj) series 
```matlab
xt = cell(1,N);
ut = cell(1,N);

for j = 2:N+1
    % threshold excesses follow gpd, thus Xt - threshold!
    xt{j-1} = my_excesses{j}
    ut{j-1} = my_covariates{j}
end
```

#### Models Configuration  
The framework is configured in adjustCFG.m file: 

```matlab
%%%%%%%%%%%%%%%%%%%%%%%
%choose the model
cfg = FEM.getdefaultCfg('MCMC','SPATIAL_BV_Mat'); % no Gurobi, attention be very slow!
cfg.model.subtype = 'GPDcpp';

% options
%cfg = FEM.getdefaultCfg('MCMC','parforSPATIAL_BV_Gurobi'); % with Gurobi and no spatial regularizaiton
%cfg = FEM.getdefaultCfg('MCMC','SPATIAL_BV_Gurobi'); % with Gurobi and spatial regularizaiton
%cfg.model.subtype = 'GPDmatlab'; % if gsl or mexcomplier can not be used, using matlab will be very slow!

%%%%%%%%%%%%%%%%%%%%%%%%%
% GPD options
% expected same number of combinations in covariates for ksi, si
cfg.model.ut_combis{1} = {[3 5 7 8 11]}; % cell of arrays: e.g. {[1 2 3], [2 3 4], []} ...
cfg.model.ut_combis{2} = {[3 5 7 8 11]}; 

%%%%%%%%%%%%%%%
% FEM options
cfg.modelFit.subspaceEps = 1e-03; % float
cfg.modelFit.annealing = 100; % int
cfg.modelFit.subspaceItr = 1000; % int
cfg.modelFit.Klist = [1 2];   % array
cfg.modelFit.Ctlist = 10:10:20; % array
cfg.modelFit.Cnlist = inf; % inf for BV indicates: no spatial regularization

```
Further it is possible to set parameters for the switching process optimization. For instance
you can set the number of finite elements used for discretization of Γ(s, t), 
where cfg.gamma.num_fe = NT. If xt is very long, setting NT << T, will reduce the computational time. But, in this case a model switching can happen only every len(xt)/num_fe. If the 

```matlab
%%%%%%%%%
% Gamma options 
cfg.gamma.num_fe   = 600; % number of finite elemnets to use. 
cfg.gamma.discrete = [];
cfg.gamma.distMatrix = distMatrix;
cfg.gamma.positionSwitch = [];

```


#### Choosing the best model
The best model is chosen wrt. Information Criteria. Available criteria are AIC, AICc, BIC. 
In the application we used AICc, to adjust please change the string 'AICc' in line to 'AIC' or 'BIC'

```matlab
optRes = rHandler.getBest('AICc');
```

#### PostAnalysis
The result handler is an FEM-Framework object, so make sure FEM-Framework in your MATLAB path
```matlab
addpath(genpath('PATH_TO_FEM_Framework'));
```

The framework contains some functionality for post-analysis. For details please check
applications/swiss_precipitation/postanalysis/inference_QL98.m. 

For PostAnalysis the inference Handler needs to be initialized: 

```matlab
infHandler = Models.MCMC.GPD.Inference(optRes, iHandler, locNames, timeCell, coord, []);
```
hereby optRes, iHandler are obtained from the running the framework (check main_spatial.m), 
locNames, timeCell and coord are each a cell of length N (number of locations):
* locNames (names of the locations), 
* timeCell (contains the time of each events for every location),
* coord a Nx2 cell contains the latitude/longitude coordinates for each location. 
Once the infHandler is initialized, it can be used for post analysis:

```matlab
%% get confidence intervalls for theta via Fischer Matrix
conf_int = infHandler.getConf4Theta();


%% temporal evalaution of the parameters for a locations
[ksij, sigmaj, Gammaj] = infHandler.getTheta_s(index_loc);

%% estimate the synchronization measures for inferring the spatial dependence structure
[delayES, strengthES] = infHandler.computeES(0);
[delayES_cl1, strengthES_cl1] = infHandler.computeES(1);
[delayES_cl2, strengthES_cl2] = infHandler.computeES(2);

```



#### How to Cite
O. Kaiser, D. Igdalov, O. Martius,  and I. Horenko. On Computationally-Scalable Spatio-Temporal Regression Clustering of Precipitation Threshold Excesses. **arxiv**.


#### Help, it's not working!!
1. if you see "Undefined variable "Core" or class "Core.Datatypes.SGamma.computeS"."
   you did not setup the path in the proper way or restarted MATLAB without saving the path, check Matlab path


2. email us

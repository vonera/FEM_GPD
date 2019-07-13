## FEM-BV-GPD
We present here the FEM-BV-GPD framework for a data-driven spatio-temporal regression based clustering of threshold excesses.



### Functionality: 
Given a series of threshold excesses and covariates for multiple locations the framework provides an optimal sparse solution for the spatio-temporal clustering problem by incorporating information theory based model selection criteria and L1-regularized regression. The framework goes beyond a priori probabilistic and deterministic assumptions typical for standard approaches like Hidden Markov Models (HMM) and Gaussian Mixture Models (GMM). 


The nonstationarity and the non-homogeneity of the data are resolved by describing the underlying behavior by a set of `$K\leq 1$` locally stationary models and a spatial nonstationary switching process `$Γ(s, t)$`, where `$K > 1$` indicates the existence of systematically missing covariates.


The realization of framework allows to consider simultaneously different configurations of the problem referring to hyper-parameter tuning of:
	* the number of clusters,
	* number of switches,
	* different configuration of covariates
	* different L1-regularizations 

The final output of the framework is 



The resulting spatial FEM-BV-GPD can be employed as a robust exploratory sparse regression analysis tool for spatio-temporal extremes, that enables to iden- tify the significant resolved covariates and to account for the influence from the systematically missing ones.

#### Spatial dependence structure 
The resulting switching process `$Γ(s, t)$` provides a posteriori a pragmatic description of the underlying spatial dependence structure by grouping together all locations that exhibit similar behavior in `$Γ(s, t)$`, e.g., by estimating the Event Synchronization (ES) measure matrix [47] which allows to identify contiguous regions that exhibit similar spatio-temporal behavior. The estimation of the ES-matrix is included in the framework.

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
	 and execute these three lines after each other in the MATLAB command line.
	* Go to FEM_GPD/demo and run main_demo.m
	* In case you use Gurobi, uncomment the following lines in main_demo and adjust the PATH_TO_GUROBI (check gurobi instructions)
		 path = pwd;
		 cd PATH_TO_GUROBI/mac64/matlab
		 gurobi_setup
		 cd(path)
	* Run main_demo.m



### Usage

In Brief: Follow main_spatial.m in FEM_GPD/applications/swiss_precipitation.
In more detail: 
	* Create a folder my_application under FEM_GPD/applications
	* copy adjustCFG.m and main_spatial.m from FEM_GPD/applications/swiss_precipitation to FEM_GPD/applications/my_application
	* adjust adjustCFG.m (refer to wiki for more details on the configuration file)
	* adjust main_spatial.m wrt to data loading
	* run main_spatial.m
	* check wiki for a howoto for PostAnalysis
	


#### Data 
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
%%%%%%%%%%%%%%%%%%%%%%%%%
% expected same number of combinations in covariates for ksi, si
cfg.model.ut_combis{1} = {[3 5 7 8 11]}; % cell of arrays: {[1 2 3] [2 3 4]} ...
cfg.model.ut_combis{2} = {[3 5 7 8 11]}; 

%%%%%%%%%%%%%%%
% FEM options
cfg.modelFit.subspaceEps = 1e-03; % float
cfg.modelFit.annealing = 100; % int
cfg.modelFit.subspaceItr = 1000; % int
cfg.modelFit.Klist = [1 2];   % array
cfg.modelFit.Ctlist = 10:10:20; % array
cfg.modelFit.Cnlist = inf; % please note inf for BV indicates: no spatial regularisation
```matlab


#### Choosing the best model
The best model is chosen wrt. Information Criteria. Available criteria are AIC, AICc, BIC. 
In the application we used AICc, to adjust please change the string 'AICc' in line to 'AIC' or 'BIC'
```matlab
optRes = rHandler.getBest('AICc');
```

#### PostAnalysis
The framework contains some functionality for post-analysis.

##### Statistical Significance
* Fischer Matrix
* QQ plots

##### Spatial dependence structure



#### Citation
[FEM-BV-GEV]
[ES- 47]


#### Help, it's not working!!
1. if you see "Undefined variable "Core" or class "Core.Datatypes.CGamma.computeC"."
   you did not setup the path in the proper way or restarted MATLAB without saving the path (check  "1) Matlab path" )


3. email us
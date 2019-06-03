Downloading:
 


Installation/configuration:
	the downloaded FEM_Framework directory contains subdirectories
		examples (simple showcase tests, run them first after configuration of your environment)
		FEM_Framework (framework code, Core components and Models)
		tests (unit tests, and unit test framework, read section "unit test" of this doc)
		workdir (here is a nice place to work, the content of workdir is not under SVN, you can add it to your on repo)	  
	  
	1) Matlab path
		add "FEM_Framework" to the Matlab path (without subfolders)
		optionally add "examples" and "workdir" to the Matlab path (without subfolders) 
		(to make the changes to the Matlab path permanent, call from Matlab console pathtool and cklick the Save button)
	
	2) compiling C++ extensions
		The Framework provides different implementation of Gamma estimators.
		ADAPTIVE  (Fast heuristic only for discrete Gamma, default in most examples and unit test, need to be compiled) 
		LINEAR  (Slow, can both discrete and continues Gamma) 
		LINEARJas  (Slow, can both discrete and continues Gamma, use this one in case the Adaptive is not possible) 
		
		to be able to run the examples just out of the box you should compile the C++ extension of the ADAPTIVE Gamma estimator.
		- go into ..../FEM_Framework/+Core/+Gamma/+Adaptive and run from the Matlab console 
       		mex mexCoarseGrainGamma.cpp
       	 	mex mexFindOptimalGammaImproved.cpp

		
Running the examples:
		After completing of the (Installation/configuration) section you should be able to run the code in examples directory
		In case you was not able to compile the Adaptive Gamma extension you can switch to an another Gamma estimator by changing of 
		one statement in a concrete example the uses Adaptive Gamma. 
		e.g example_VARX.m:
				change the statement "cfg = FEM.getdefaultCfg('VARX','ADAPTIVE');" to cfg = FEM.getdefaultCfg('VARX','LINEARJas');
		and run the code.
		
		
		

Running Unit tests 
	1) add the unittest framework to the matlab path 
		add .../tests/matlab_xunit/xunit to the matlab path  (without subfolders)
    2) add FEM_WORKSPACE/tests/Core to the matlab path (with subfolders)
	3) add FEM_WORKSPACE/tests/Models to the matlab path (with subfolders)
	4) (for now) run manually the test of interest
		
		
Help, it's not working!!
1 if you see "Undefined variable "Core" or class "Core.Datatypes.CGamma.computeC"."
   you did not setup the path in the proper way or restarted matlab without saving the path (check  "1) Matlab path" )

2) if you see "Undefined variable "Core" or class "Core.Gamma.Adaptive.mexFindOptimalGammaImproved"
   you try to use  Adaptive Gamma solver but the C++ extension is not compiled (check the points 
   "2) compiling C++ extensions" or "Running the examples:" )
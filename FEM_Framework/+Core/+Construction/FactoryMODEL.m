%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : FactoryMODEL.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description :
%
% The development was guided by Prof. Dr. Illia Horenko head of the group
%    http://icsweb.inf.unisi.ch/cms/index.php/groups/20-illia-horenko.html
% USI-ICS Horenko group
%    http://icsweb.inf.unisi.ch/cms/index.php/groups/4-group-horenko.html
%--------------------------------------------------------------------------
% LICENSE
% (c) 2015, D.Igdalov, O.Kaiser USI-ICS Horenko Group, All Rights Reserved
%--------------------------------------------------------------------------
% This Software Framework may not be distributed in parts or its entirety
% without prior written agreement by O.Kaiser or D.Igdalov.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
% DEALINGS IN THE SOFTWARE.
%-------------------------------------------------------------------------
classdef FactoryMODEL
    %FACTORYMODEL this is an interface class
    %   to add a new MODEL to FEM Framework you have to implement this
    %   interface..
    
    properties
        mytype
    end
    
    methods
        %------------------------------------------------------------------
        % (String) mytype is the modelname used in config (this should be
        % upper case)
        function self = FactoryMODEL(mytype)
            
            self.mytype = mytype;
        end
        
        %------------------------------------------------------------------
        function checkTYPE(s, type)
            
            if ~strcmp(type, s.mytype)
                error(['expected model type ',s.mytype, ', got ', type]);
            end
            
        end
        
        %------------------------------------------------------------------
        function mfit = createMFit(self, cfg) % class(cfg) -> FEMConfig
            
            self.checkTYPE(cfg.model.type);
            
            mfit = Core.ModelFit(cfg.modelFit);
            
        end
        
        %------------------------------------------------------------------
        function rHandler = createResultHandler(self, cfg)
            
            self.checkTYPE(cfg.model.type);
            
            rHandler = Core.Datatypes.ResultHandler();
            
        end
        
        %------------------------------------------------------------------
        function iHandler = createInputHandler(self, cfg, xt, ut)
            
            self.checkTYPE(cfg.model.type);
            
            iHandler = Core.Datatypes.InputHandler(cfg, xt, ut);
   
        end
        
        %------------------------------------------------------------------
        function input = createInput(self, cfg, xt, ~, meta)
            
            self.checkTYPE(cfg.model.type);
            
            input = Core.Datatypes.Input(xt, meta);
   
        end
        
        %------------------------------------------------------------------
        function gammaEst = createGEstimator(self, cfg, input, K, Ct, Cn) % gamma
            
            self.checkTYPE(cfg.model.type);
            
            if strcmp(cfg.gamma.type, 'SPATIAL_BV_Mat')
                gammaEst = Core.Gamma.Spatial.BV.Mat.GEstimator(cfg.gamma, input, K, Ct, Cn);
            elseif strcmp(cfg.gamma.type, 'SPATIAL_BV_Gurobi')
                gammaEst = Core.Gamma.Spatial.BV.Gurobi.GEstimator(cfg.gamma, input, K, Ct, Cn);
            elseif strcmp(cfg.gamma.type, 'parforSPATIAL_BV_Gurobi')
                gammaEst = Core.Gamma.Spatialparfor.BV.Gurobi.GEstimator(cfg.gamma, input, K, Ct, Cn);
            else 
                error(['something wrong with gammaAlg ' cfg.gamma.type ...
                    ' expected { SPATIAL_BV_Mat,'...
                    'SPATIAL_BV_Gurobi, parforSPATIAL_BV_Gurobi}']);
            end
        end
        
        %------------------------------------------------------------------
        function mfcfg = createExapmleMFitCFG(~)
            
            mfcfg = Core.Construction.MFitConfig(10,false,10,0.00000005,[1,2,3],[5,12],[0]);
        
        end
        
        %------------------------------------------------------------------
        function gcfg = createExapmleGammaCFG(~, Gamma)
            
            if strcmp(Gamma, 'SPATIAL_BV_Mat')
                gcfg = Core.Construction.GammaSConfig('SPATIAL_BV_Mat', 99, 'NO', 'NO');
            elseif strcmp(Gamma, 'SPATIAL_BV_Gurobi')
                gcfg = Core.Construction.GammaSConfig('SPATIAL_BV_Gurobi', 99, 'NO', 'NO');
            elseif strcmp(Gamma, 'parforSPATIAL_BV_Gurobi')
                gcfg = Core.Construction.GammaSConfig('parforSPATIAL_BV_Gurobi', 99, 'NO', 'NO');     
            else
                error(['bad Gamma, expected SPATIAL_BV_Mat ',...
                    'SPATIAL_BV_Gurobi, parforSPATIAL_BV_Gurobi'])
            end
            
        end
           
    end
    
    methods (Abstract)
        
        % this method should return an object of a class that is a
        % specialized implementation of the FEMConfigModel class
        cfg = createModelCFGfromStruct(self, struct_cfg)
        
        % this method should return an object of a class that is a
        % specialized implementation of the FEMConfigModel class
        cfg = createExapmleModelCFG(type)
        
        % this method should return an object of a class that is
        % specialized implementation of the FEMParamEstimator class
        pest = createParamEstimator(self, cfg, input)
        
    end
end


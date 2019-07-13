%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : Factory.m
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
classdef Factory
    %FACTORY for FEM workflow
    %   everything static here 
    
    methods(Static)
        
        %------------------------------------------------------------------
        function mfac = createModelFactory(type)
            try
                mfac = eval(['Models.',type,'.Factory()']);
            catch err
                disp(err.message);
                error(['Can not create Factory for ',type, ', check if Models.',type, '.Factory class exist']);
            end
        end
        
        %------------------------------------------------------------------
        function mfit = createModelFit(cfg) % FEMConfig

            mfac = Core.Construction.Factory.createModelFactory(cfg.model.type);
            
            mfit = mfac.createMFit(cfg);
            
        end
                
        %------------------------------------------------------------------
        function pest = createParamEstimator(cfg, input) % model

            mfac = Core.Construction.Factory.createModelFactory(cfg.type);
            
            pest = mfac.createParamEstimator(cfg, input);
            
        end
        
        %------------------------------------------------------------------
        function gammaEst = createGammaEstimator(cfg, input, K, Ct, Cn)
            
            mfac = Core.Construction.Factory.createModelFactory(cfg.model.type);
                                    gammaEst = mfac.createGEstimator(cfg, input, K, Ct, Cn);
        
        end
               
        %------------------------------------------------------------------
        function iHandler = createInputHandler(cfg, xt, ut) % FEMConfig

            mfac = Core.Construction.Factory.createModelFactory(cfg.model.type);
            
            iHandler = mfac.createInputHandler(cfg, xt, ut);
            
        end
        
        %------------------------------------------------------------------
        function rHandler = createResultHandler(cfg) % FEMConfig

            mfac = Core.Construction.Factory.createModelFactory(cfg.model.type);
            
            rHandler = mfac.createResultHandler(cfg);
            
        end
        
        %------------------------------------------------------------------
        function input = createInput(cfg, xt, ut, meta) % FEMConfig

            mfac = Core.Construction.Factory.createModelFactory(cfg.model.type);
            
            input = mfac.createInput(cfg, xt, ut, meta);
            
        end
        
       
        %------------------------------------------------------------------
        %
        % config creation         
        %------------------------------------------------------------------
        
        %------------------------------------------------------------------
        function cfg = createFEMConfig(cfg_Struct)
            
            if ~isfield(cfg_Struct,'modelFit')
                error('bad cfgStruct, "modelFit" field expected');
            end
            modelFit = Core.Construction.Factory.createMFitCfg(cfg_Struct.modelFit);
            
            if ~isfield(cfg_Struct,'gamma')
                error('bad cfgStruct, "gamma" field expected');
            end 
            gamma = Core.Construction.Factory.createGammaCfg(cfg_Struct.gamma);
            
            if ~isfield(cfg_Struct,'model')
                error('bad cfgStruct, "gamma" field expected');
            end 
            model = Core.Construction.Factory.createModelCfg(cfg_Struct.model);

            cfg = Core.Construction.FEMConfig(modelFit, gamma, model);
        end
   
        %------------------------------------------------------------------
        function cfg = createMFitCfg(cmfit)

            cfg = Core.Construction.MFitConfig(...
                                    cmfit.annealing,...
                                    cmfit.subspaceItr,...
                                    cmfit.subspaceEps,...
                                    cmfit.Klist,...
                                    cmfit.Ctlist,...
                                    cmfit.Cnlist);
        end
        
        %------------------------------------------------------------------
        function cfg = createGammaCfg(cgam)

            if strcmp(strcmp(cgam.type, 'SPATIAL_BV_Mat') ||...
                    strcmp(cgam.type, 'SPATIAL_BV_Gurobi')||...
                    strcmp(cgam.type, 'parforSPATIAL_BV_Gurobi'))
                cfg = Core.Construction.GammaSConfig(cgam.type,...
                    cgam.num_fe,...
                    cgam.options,...
                    cgam.discrete,...
                    cgam.diagnostics);
            else 
                error(['unknown gamma.type ' cgam.type ...
                    ' expected {SPATIAL_BV_Mat, SPATIAL_BV_Gurobi, parforSPATIAL_BV_Gurobi}']);
            end
        end
         
        %------------------------------------------------------------------
        function cfg = createModelCfg(cmod)

            mfac = Core.Construction.Factory.createModelFactory(cmod.type);
            
            cfg = mfac.createModelCFGfromStruct(cmod);

        end
              
        %------------------------------------------------------------------
        function cfg = createExampleFEMConfig(Model, Gamma)
            
            mfac = Core.Construction.Factory.createModelFactory(Model);
            
            mfcfg = mfac.createExapmleMFitCFG();
            
            gcfg = mfac.createExapmleGammaCFG(Gamma);
            
            mcfg = mfac.createExapmleModelCFG(Model);
            
            cfg = Core.Construction.FEMConfig(mfcfg, gcfg, mcfg);
            
        end

    end
end


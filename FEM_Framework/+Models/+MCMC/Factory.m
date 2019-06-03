%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : Factory.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : factory configuration
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
classdef Factory < Core.Construction.FactoryMODEL
    %FACTORY Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        
        %------------------------------------------------------------------
        function self = Factory()
            
            self = self@Core.Construction.FactoryMODEL('MCMC');
        end
        
        %------------------------------------------------------------------
        function cfg = createModelCFGfromStruct(self, struct_cfg) % model
            
            self.checkTYPE(struct_cfg.type);
                        
            cfg = Models.MCMC.Config(struct_cfg.type,...
                struct_cfg.subtype,...
                struct_cfg.proposMaxTimes,...
                struct_cfg.samples,...
                struct_cfg.lassoFLAG,...
                struct_cfg.ridgeFLAG,...
                struct_cfg.gaussFLAG,...
                struct_cfg.eps,...
                struct_cfg.init_beta,...
                struct_cfg.beta_factor,...
                struct_cfg.init_noise,...
                struct_cfg.noise_more,...
                struct_cfg.noise_less,...
                struct_cfg.ut_combis, ...
                struct_cfg.paramTHRESHOLD);
            
        end
        
        %------------------------------------------------------------------
        function pest = createParamEstimator(self, cfg, input) % model
            
            self.checkTYPE(cfg.type);
            
            if strcmp(cfg.subtype, 'GEVmatlab')
                pest = Models.MCMC.GEV.PEstimator(cfg, input);
            elseif strcmp(cfg.subtype, 'GEVcpp')
                pest = Models.MCMC.GEV.MCMCMetFNT.PEstimator(cfg, input);
            elseif strcmp(cfg.subtype, 'GPDmatlab')
                pest = Models.MCMC.GPD.PEstimator(cfg, input);
            elseif strcmp(cfg.subtype, 'GPDcpp')
                pest = Models.MCMC.GPD.GPDcpp.PEstimator(cfg, input);
            else
                error(['something wrong with MCMC Model SubType ' cfg.subtype ...
                    ' expected {GEVmatlab, GEVcpp, GPDmatlab, GPDcpp}']);
            end
        end
        
        %------------------------------------------------------------------
        function iHandler = createInputHandler(self, cfg, xt, ut) % model
            
            self.checkTYPE(cfg.model.type);
            
            if strcmp(cfg.model.subtype, 'GEVmatlab')
                iHandler = Models.MCMC.GEV.InputHandler(cfg, xt, ut);
            elseif strcmp(cfg.model.subtype, 'GEVcpp')
                iHandler = Models.MCMC.GEV.InputHandler(cfg, xt, ut);
            elseif strcmp(cfg.model.subtype, 'GPDmatlab')
                iHandler = Models.MCMC.GPD.InputHandler(cfg, xt, ut);
            elseif strcmp(cfg.model.subtype, 'GPDcpp')
                iHandler = Models.MCMC.GPD.InputHandler(cfg, xt, ut);
            else
                error(['InputHandler for SubType ' cfg.subtype ' not supported']);
            end
        end
        
        %------------------------------------------------------------------
        function input = createInput(self, cfg, meta, xt, ut) % model
            
            self.checkTYPE(cfg.model.type);
            
            if strcmp(cfg.model.subtype, 'GEVmatlab')
                input = Models.MCMC.GEV.Input(meta, xt, ut);
            elseif strcmp(cfg.model.subtype, 'GEVcpp')
                input = Models.MCMC.GEV.Input(meta, xt, ut);
            elseif strcmp(cfg.model.subtype, 'GPDmatlab')
                input = Models.MCMC.GPD.Input(meta, xt, ut);
            elseif strcmp(cfg.model.subtype, 'GPDcpp')
                input = Models.MCMC.GPD.Input(meta, xt, ut);
            else
                error(['Input for SubType ' cfg.subtype ' not supported']);
            end
        end
        
        %------------------------------------------------------------------
        function cfg = createExapmleModelCFG(self, type) % string
            self.checkTYPE(type);
            cfg = Models.MCMC.Config(type,'GEVmatlab',100, 500, 0, 0, false, 0.0005,1,1.5,0.01,1.01,0.1,{{[1,2]}, {[1,2]}, {[1,2]}}, 0);
        end
        
    end
end


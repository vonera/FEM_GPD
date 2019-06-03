%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : FEM.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : get defaults, input, fit
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

classdef FEM
    %FEM Summary of this class goes here
    %   Detailed explanation goes here
    
    
    methods(Static)
        
        %------------------------------------------------------------------
        function listmodels()
            disp('Models: MCMC ');
            disp('Gammas: SBV ');
        end
        
        %------------------------------------------------------------------
        function cfg = getdefaultCfg(Model, Gamma)
            cfg = Core.Construction.Factory.createExampleFEMConfig(Model, Gamma);
        end
        
        %------------------------------------------------------------------
        function inputHandler = createInputHandler(cfg, xt, ut)
            if ~isa(cfg, 'Core.Construction.FEMConfig')
                error('cfg expected to be object of Class FEMConfig, create it using FEM.getdefaultCfg(MODEL)');
            end
            inputHandler = Core.Construction.Factory.createInputHandler(cfg, xt, ut);
        end
        
        %------------------------------------------------------------------
        function resultHandler = fit(inputHandler)
            if ~isa(inputHandler, 'Core.Datatypes.InputHandler')
                error('input expected to be object of Class Core.Datatypes.InputHandler or Subclass, create it using FEM.createInput()');
            end
            mf = Core.Construction.Factory.createModelFit(inputHandler.cfg);
            resultHandler = mf.fit(inputHandler);
        end
    end
end


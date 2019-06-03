%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : GEstimator.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : Gurobi Estimate
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
classdef GEstimator < Core.Gamma.Spatialparfor.BV.GEstimator
    %GESTIMATORSBV Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        model       %Gurobi structure for model, containing Aeq, Aneq ...
        params      %Gurobi structure for params
    end
    
    methods
        
        %------------------------------------------------------------------
        function self = GEstimator(cfg, input, K, Ct, Cn) % gamma cgf
            self = self@Core.Gamma.Spatialparfor.BV.GEstimator(cfg, input, K, Ct, Cn);
            
            
            
            self.params.seed = randi(1,1);
            
            if any(strcmpi(self.cfg.discrete,{'I','B','C'}))
                self.model.vtype = self.cfg.discrete;
            elseif ~isempty(self.cfg.discrete)
                error('Wrong self.cfg.discrete, shold be I or B or C')
            end
            
            self.params.outputflag = self.cfg.diagnostics;
            
        end
        
        
        %------------------------------------------------------------------
        function [gamma] = estimate_local(self, resid, gamma_old, j)
            
            % create the model structure for Gurobi solver
            self.model.A = [self.Aneq{j}; self.Aeq{j}];
            self.model.rhs = [full(self.bneq{j}); full(self.beq{j})];
            self.model.sense = [repmat('<', size(self.Aneq{j},1), 1); repmat('=', size(self.Aeq{j},1), 1)];
            self.model.lb = -inf(size(self.Aneq{j}, 2),1);
            
            self.model.obj = full(self.constructRes_vec(resid, j));
            
%             % set the startvalue as the olg gamma
%             tmp = gamma_old(self.pos_switch(2:end),:);
%             self.model.start = [reshape(tmp, self.tBins*self.K, 1) ; zeros(self.K*(self.tBins-1),1)];
%             
%             % set the cutoff, if objective val < cutoff : return solution
%             self.model.Cutoff = self.model.obj'*self.model.start;
            
            result = gurobi(self.model, self.params);
            
            if ~strcmp(result.status,'OPTIMAL')|| strcmp(result.status,'ITERATION_LIMIT')
                throw(MException('ResultChk:InvalidState','gurobi failed %d ',result.status));
            end
            
            gamma = self.reconstructGamma(result.x, j);
        end
        
        
        
        
    end
end



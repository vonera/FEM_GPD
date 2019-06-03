%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : GEstimator.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : solve BV via Gurobi
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
classdef GEstimator < Core.Gamma.Spatial.BV.GEstimator
    %GESTIMATORSBV Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        model       %Gurobi structure for model, containing Aeq, Aneq ...
        params      %Gurobi structure for params
    end
    
    methods
        
        %------------------------------------------------------------------
        function self = GEstimator(cfg, input, K, Ct, Cn) % gamma cgf
            self = self@Core.Gamma.Spatial.BV.GEstimator(cfg, input, K, Ct, Cn);
            
            % create the model structure for Gurobi solver
            self.model.A = [self.Aneq; self.Aeq];
            self.model.rhs = [full(self.bneq); full(self.beq)];
            self.model.sense = [repmat('<', size(self.Aneq,1), 1); repmat('=', size(self.Aeq,1), 1)];
            self.model.lb = -inf(size(self.Aneq, 2),1);
            
            if any(strcmpi(self.cfg.discrete,{'I','B','C'}))
                self.model.vtype = self.cfg.discrete;
            elseif ~isempty(self.cfg.discrete)
                error('Wrong self.cfg.discrete, shold be I or B or C')
            end
            
            self.params.outputflag = self.cfg.diagnostics;
            
        end
        
        %------------------------------------------------------------------
        function [optRes] = estimate(self, optRes)

            if self.K ==1
                optRes.hidden = optRes.hidden.getStatGamma(self.Tlist);
                return
            end
            
            self.model.obj = full(self.constructRes_vec(optRes.resid));
            
            % set gurobi parameters
%             self.model.cutoff = optRes.acf - 1000*eps;
%             self.model.BarConvTol = 0;
%             
%             % use prev value as start
%             tmp = optRes.hidden.gamma(self.pos_switch(2:end),:);
%             self.model.start = [reshape(tmp,self.tBinsNK,1); zeros(self.tBins_1NK,1) ];
            
            
            result = gurobi(self.model, self.params);
       
            if ~strcmp(result.status,'OPTIMAL')|| strcmp(result.status,'ITERATION_LIMIT')
                throw(MException('ResultChk:InvalidState','gurobi failed %d ',result.status));
            end
            
            gamma = self.reconstructGamma(result.x);
            
            % proofe if gamma is apropriate
            if min(sum(gamma)) == 0
                throw(MException('ResultChk:InvalidState','Resulting Gamma for K = %d Ct = %d Cn = %d has empty cluster for all locations', self.K, self.Ct, self.Cn));
            end
            
            hidden = Core.Datatypes.CGamma(gamma, self.Tlist, self.Ct, self.Cn, self.cfg.discrete);
            
%             tmp_x = result.x(1:self.tBinsNK);
%             tmp_resid = self.model.obj(1:self.tBinsNK);
%             new_acf = tmp_x'*tmp_resid;
%             
            new_acf = optRes.param.computeAcf(hidden, optRes.resid);
            
            % check acf convergency
%             if optRes.acf < new_acf
%                 throw(MException('ResultChk:InvalidState', 'Growing ACF detected, acf %d befor gamma opt., after %d', optRes.acf, new_acf));
%             end
            
            % check acf convergency
            if (optRes.acf - new_acf) > -1e-09
                warning('ResultChk:InvalidState','acf growing slightly')
            else
                throw(MException('ResultChk:InvalidState', 'Growing ACF detected, acf %d befor gamma opt., after %d', optRes.acf, new_acf));
            end
            
            
            
            optRes.acf = new_acf;
            optRes.hidden = hidden;
            
        end
        
        %------------------------------------------------------------------
        function hidden = generateRandomGamma(self)
            
            if self.K == 1
                hidden = Core.Datatypes.CGamma.getStatGamma(self.Tlist);
            else
                max_trys = 50;
                for i = 1:max_trys
                    tmp_vtype = self.model.vtype;
                    self.model.vtype = 'C';
                    
                    self.model.obj = full(self.constructRes_vec(rand(sum(self.Tlist), self.K)));
                    result = gurobi(self.model, self.params);
                    
                     self.model.vtype = tmp_vtype;
                    
                    if ~strcmp(result.status,'OPTIMAL')|| strcmp(result.status,'ITERATION_LIMIT')
                        throw(MException('ResultChk:InvalidState','gurobi failed %d ',result.status));
                    end
                    
                    gamma = self.reconstructGamma(result.x);
                    
                    if all(sum(gamma))
                        break;
                    end
                end
                
                if ~all(sum(gamma))
                    throw(MException('ResultChk:InvalidState','Initial random Gamma for K = %d Ct = %d Cn = %d has empty cluster for all locations', self.K, self.Ct, self.Cn));
                end
                hidden = Core.Datatypes.CGamma(gamma, self.Tlist ,self.Ct, self.Cn, self.cfg.discrete);
            end
            
        end
        
        
    end
end



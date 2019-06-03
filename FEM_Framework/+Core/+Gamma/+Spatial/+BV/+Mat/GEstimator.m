%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : GEstimator.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : solve BV via matlab
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
        
        options
    end
    
    methods
        
        %------------------------------------------------------------------
        function self = GEstimator(cfg, input, K, Ct, Cn) % gamma cgf
            self = self@Core.Gamma.Spatial.BV.GEstimator(cfg, input, K, Ct, Cn);
            
            warning('Solver:GEstimator','Please use the much faster Gurobi if possible!')
            %             if ~isempty(self.cfg.discrete)
            %                 error('Only continues Gamma supporten: -> cfg.gamma.discrete = false')
            %             end
            
            self.options = self.cfg.options;
            self.options = optimoptions('intlinprog','Display','off');
        end
        
        %------------------------------------------------------------------
        function [optRes] = estimate(self, optRes)
            
            if self.K ==1
                optRes.hidden = optRes.hidden.getStatGamma(self.Tlist);
                return
            end
            
            res_vec = self.constructRes_vec(optRes.resid);
            
            if self.cfg.discrete
                [x,~,exitflag] = intlinprog(res_vec, [1:length(res_vec)], self.Aneq, self.bneq, self.Aeq, self.beq,[],[],self.options);
            else
                [x,~,exitflag] = linprog(res_vec, self.Aneq, self.bneq, self.Aeq, self.beq,[],[],[],self.options);
            end
            % prevent over/under-estimation (order of machine precision)
            x(x<0) = 0;
            x(x>1) = 1;
            
            if exitflag ~= 1
                warning('ResultChk:LinpErr','Core.Gamma.GEstimatorSBV: exitflag in linprog %d ', exitflag)
            elseif exitflag < 0
                throw(MException('ResultChk:InvalidState','error in linprog %d ',exitflag));
            end
            
            gamma = self.reconstructGamma(x);
            
            % proofe if gamma is apropriate
            if min(sum(gamma)) == 0
                throw(MException('ResultChk:InvalidState','Resulting Gamma for K = %d Ct = %d Cn = %d has empty cluster for all locations', self.K, self.Ct, self.Cn));
            end
            
            hidden = Core.Datatypes.CGamma(gamma, self.Tlist, self.Ct, self.Cn, self.cfg.discrete);
            new_acf = optRes.param.computeAcf(hidden, optRes.resid);
            
            % check acf convergency
            if optRes.acf < new_acf
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
                    res_vec = self.constructRes_vec(rand(sum(self.Tlist), self.K));
                    if self.cfg.discrete
                        x = intlinprog(res_vec, [1:length(res_vec)], self.Aneq, self.bneq, self.Aeq, self.beq,[],[],self.options);
                    else
                        x = linprog(res_vec, self.Aneq, self.bneq, self.Aeq, self.beq,[],[],[],self.options);
                    end
                    
                    gamma = self.reconstructGamma(x);
                    
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



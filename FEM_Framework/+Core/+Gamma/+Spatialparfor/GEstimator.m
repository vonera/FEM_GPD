%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : GEstimator.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : cfg for spatial parfor gamma
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
classdef GEstimator < Core.Gamma.GEstimator
    %GESTIMATOR this is the initializing class for spatial FEM for estimating
    %the optimal Gamma.
    %   Here we init some parameters, pos_switch, distMatrix
    
    properties
        
        distMatrix      % matrix of size (N x N)
        pos_switch      % vector of size (1, tBins+1 ).
        
        Cn              % double
        N               % integer
        NT              % integer
        tBins           % array of integers
        tBinsK          % array of integer

    end
    
    methods
        
        %------------------------------------------------------------------
        function self = GEstimator(cfg, input, K, Ct, Cn) % gamma cgf
            self = self@Core.Gamma.GEstimator(cfg, input, K, Ct);
            
            self.Cn = Cn;
            self.N = length(self.input(:));
            self.NT = sum(self.Tlist);
            if cfg.num_fe <= min(self.Tlist)
                self.tBins = (cfg.num_fe-1)*ones(1,self.N);
            else
                self.tBins = self.Tlist-1;
            end
            self.tBinsK = self.tBins .* K;
            self.pos_switch = self.getPositionSwitch(cfg);
            
            
            self.distMatrix = cfg.distMatrix;
            if ~isempty(self.distMatrix)
                warning('distMatrix is ignored in the parfor version')
            end
        end
    end
    
    methods
        
        %------------------------------------------------------------------
        function [optRes] = estimate(self,optRes)
            
            
            if self.K ==1
                optRes.hidden = optRes.hidden.getStatGamma(self.Tlist);
                return
            else
                gamma_cell = cell(self.N,1);
                resid_all = self.resid2cell(optRes.resid);
                gamma_old = self.resid2cell(optRes.hidden.gamma);
                
                parfor j = 1:self.N
                %for j = 1:self.N
                    gamma_cell{j} = self.estimate_local( resid_all{j}, gamma_old{j}, j );
                end
                gamma = cell2mat(gamma_cell);
            end
            
            % proofe if gamma is apropriate
            if min(sum(gamma)) == 0
                throw(MException('ResultChk:InvalidState','Resulting Gamma for K = %d Ct = %d Cn = %d has empty cluster for all locations', self.K, self.Ct, self.Cn));
            end
            
            hidden = Core.Datatypes.CGamma(gamma, self.Tlist, self.Ct, self.Cn, self.cfg.discrete);
            new_acf = optRes.param.computeAcf(hidden, optRes.resid);
            
%             % check acf convergency
%             if optRes.acf < new_acf
%                 throw(MException('ResultChk:InvalidState', 'Growing ACF detected, acf %d befor gamma opt., after %d', optRes.acf, new_acf));
%             end
            
            % check acf convergency
            if (optRes.acf - new_acf) < 0 && (optRes.acf - new_acf) > -1e-09
                warning('ResultChk:InvalidState','acf growing slightly')
            elseif (optRes.acf - new_acf) <= -1e-09
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
                    Tlist_all = self.Tlist;
                    Ktmp = self.K;
                    gamma_cell = cell(self.N,1);
                    
                    parfor j = 1:self.N
                        resid = rand(Tlist_all(j), Ktmp);
                        gamma_cell{j} = self.estimate_local(resid, zeros(size(resid)), j);
                    end
                    gamma = cell2mat(gamma_cell);
                    
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
    
    methods (Access = 'protected')
        
        %------------------------------------------------------------------
        function [pos_switch] = getPositionSwitch(self, cfg)
            % POSITIONSWITCH estimates the position switches for the FE
            %If, no pos_switch is provided by the user, then aquidistant
            pos_switch = cell(1,self.N);
            
            if isempty(cfg.positionSwitch)
                for j=1:self.N
                    pos_switch{j} = unique( round( 0:self.Tlist(j)/self.tBins(j):self.Tlist(j) ) );
                end
            else
                if iscell(pos_switch) && length(pos_switch)==self.N
                    pos_switch = cfg.positionSwitch;
                    warning('include check of posSwitch')
                else
                    error('if pos_switch provided by the user, it must be a cell with N elements.');
                end
            end
        end
        
        
        %------------------------------------------------------------------
        function res_vec = constructRes_vec(self, resid, j)
            % from (N Residua in row)' Resid el R^{TNxK}
            fe_resid = zeros(self.tBins(j), self.K);
            
            for b = 1:self.tBins(j)
                fe_resid(b,:) = sum(resid( self.pos_switch{j}(b)+1:self.pos_switch{j}(b+1),:),1);
            end
            
            res_vec = [reshape(fe_resid, self.tBinsK(j), 1) ; zeros(self.K*(self.tBins(j)-1),1)];
            
        end
        
        
        %------------------------------------------------------------------
        function gamma = reconstructGamma(self, x, j)
            % to N Gammas under each other in "classic" orientation
            fe_gamma = reshape(x(1:self.tBinsK(j)), self.tBins(j), self.K);
            
            gamma = zeros(self.Tlist(j), self.K);
            
            for b = 1:self.tBins(j)
                gamma(self.pos_switch{j}(b)+1:self.pos_switch{j}(b+1),:) = ones(self.pos_switch{j}(b+1) - self.pos_switch{j}(b), 1) * fe_gamma(b,:);
            end
            
        end
        
        
        %------------------------------------------------------------------
        function resid_cell = resid2cell(self, resid)
            
            resid_cell = cell(self.N,1);
            T_cumsum = cumsum([0 self.Tlist(1:end-1)]);
            
            for j = 1:self.N
                
                t_start = 1+T_cumsum(j);
                t_end = T_cumsum(j)+self.Tlist(j);
                
                resid_cell{j} = resid(t_start:t_end,:);
                
            end
        end
        
    end
    
    
    methods (Abstract)
        hidden = estimate_local(self, resid, gamma_old)
    end
    
end



%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : GEstimator.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : cfg for spatial gamma
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
        pos_switch      % vector of size (1, tBinsN+1 ).
        
        Cn              % double
        N               % integer
        NT              % integer
        tBins           % integer
        tBinsList       % vector of integers, list of tBins for each location
        
        
        tBinsN          % integer
        tBinsNK         % integer
        
        
    end
    
    methods
        
        %------------------------------------------------------------------
        function self = GEstimator(cfg, input, K, Ct, Cn) % gamma cgf
            self = self@Core.Gamma.GEstimator(cfg, input, K, Ct);
            
            self.Cn = Cn;
            self.N = length(self.input(:));            
            self.NT = sum(self.Tlist);
            
            [self.pos_switch, self.tBinsList] = self.getPositionSwitch(cfg);
            
            self.tBinsN = sum(self.tBinsList);
            self.tBinsNK = self.tBinsN * self.K;
            
            if all( self.tBinsList == self.tBinsList(1) )
                self.tBins = self.tBinsList(1);
            else
                self.tBins = 0;
            end
            
            self.distMatrix = cfg.distMatrix;
            if any(self.N ~= size(self.distMatrix))
                [N1,N2] = size(self.distMatrix);
                error(['Wrong size(distMatrix): ',num2str(N1),'x',num2str(N2),...
                    ' should be: ', num2str(self.N)])
            elseif ~isequal(self.distMatrix,self.distMatrix')
                error('Distance Matrix is not symmetric:')
            end
        end
    end
    
    
    methods (Access = 'protected')
        
        
        %------------------------------------------------------------------
        function [pos_switch, tBinsList] = getPositionSwitch(self, cfg)
            % POSITIONSWITCH estimates the position switches for the finite
            %elements. If, no pos_switch is provided by the user, then the
            %pos_switch is estimated aquidistant. And assigns the vector of
            %tBins to each location.
            T_cumsum = cumsum([0 self.Tlist(1:end-1)]);
            
            if isempty(cfg.positionSwitch)
                tBinsList = (cfg.num_fe-1)*ones(1,self.N);                
                pos_switch = unique(round( 0:self.Tlist(1)/tBinsList(1):self.Tlist(1) ));
                
                for j = 2:self.N
                    pos_switch_j = T_cumsum(j) + unique( round( 0:self.Tlist(j)/tBinsList(j):self.Tlist(j) ) );
                    pos_switch = [pos_switch pos_switch_j(2:end)];
                end
            else
                tBinsList = zeros(1,self.N);
                tBinsList(1) = length(cfg.positionSwitch{1})-1;
                pos_switch = cfg.positionSwitch{1};
     
                for j = 2:self.N
                    pos_switch_j = T_cumsum(j) + cfg.positionSwitch{j};
                    pos_switch = [pos_switch pos_switch_j(2:end)];
                    tBinsList(j) = length(cfg.positionSwitch{j})-1;
                end
            end
        end
        
        
    end
    
    
    methods (Abstract)
        
        hidden = estimate(self, optRes)
        
        hidden = generateRandomGamma(self)
        
    end
    
end



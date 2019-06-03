%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : Params.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : loglikelihood, number of params
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
classdef Params < Core.Datatypes.Params
    %PARAMS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        theta
    end
    
    methods
        %------------------------------------------------------------------
        function self = Params(meta, theta)
            self = self@Core.Datatypes.Params(meta);
            self.theta = theta;
        end
        
        %------------------------------------------------------------------
        function parNum = getNumPar(s)
            
            parNum = numel( s.theta(abs(s.theta) >= s.meta.paramTHRESHOLD) );
        end
        
        %------------------------------------------------------------------
        function simulated = simulate_onestep(~, ~, ~)
            simulated = [];
        end
        
        %------------------------------------------------------------------
        % overwrites
        function [ logLike, numParam ] = compute_loglike(self, hidden, input)
            
            logLike  = -self.computeAcf(hidden, self.compute_resid(input));
            numParam = 0;
        end
        
    end
    
    
    methods (Static)
        %------------------------------------------------------------------
        function param = pcell2matrix(Pcell)
            
            K =  size(Pcell,1);
            dim_param = size(Pcell{1},2);
            
            param = zeros(K,dim_param);
            for i = 1:K
                param(i,:) = Pcell{i}';
            end
        end
        
    end
    
end

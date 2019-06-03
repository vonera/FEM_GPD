%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : Params.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : likelihood, error, info
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
classdef Params
    %Params Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = 'public', SetAccess = 'protected')
        meta
    end
    
    methods
        
        %------------------------------------------------------------------
        function self = Params(meta)
            self.meta = meta;
        end

        
        %------------------------------------------------------------------
        function acf = computeAcf(~, hidden, resid)
            
            acf = sum(sum(resid .* hidden.gamma,1) ,2);
            
            %TODO remove this warning after clarifying the right use of
            %compute Acf
            warning('Params:acf','check acf estimate')
        end
        
       
        %------------------------------------------------------------------
        function [ logLike, numParam ] = compute_loglike(self, hidden, input)
            %COMPUTE_GAUSSIAN_LOGLIKE computes the loglikelihood which
            %is used for AIC
            %   Assuming that the error (one-step prediction) is Gaussian 
            %   distributed, estimate the loglikelihood for the models. 
            %   This assumption leads to additional model parameters. 
            
            error = self.compute_error(hidden, input);
            
            [dimx, T] = size(error);
            Like  = zeros(1,T);
            
            % estimate sigma
            sigma = 1/T*(error*error');
            
            % HACK to cut the bad eigenvalues (work only when the error is given)
            % reduce sigma to full rank 
            rank_sigma = rank(sigma,0.000001);
            e = eigs(sigma,rank_sigma);
            prod_eigs = prod(e);
            
            % estimate pseudo inverse
            p = pinv(sigma);
            
            % estimate Likelihood for each time step
            for t=1:T
                Like(t) = ( (2*pi)^rank_sigma * prod_eigs)^(-0.5) * exp( -0.5 * error(:,t)'*p*error(:,t) );
            end
            % Alternative stuff, maybe at some poit we will combine both
            %for t=1:T
            %    Like(t) = mvnpdf(error(:,t), zeros(dimx,1), sigma);
            %end
            %
            % estimate LogLikeli
            logLike = sum( log(Like) ,2);
            
            % estimate number of parameters for Gaussian LogLikelihood
            numParam = dimx*(dimx+1)/2;
            
        end

        
        %------------------------------------------------------------------
        function error  = compute_error(self, hidden, input)
            %COMPUTE_ERROR returns observation - simulation without start
            N = size(input,2);
            T = input(1).meta.T;
            dimx = size(input(1).xt,1);
            
            start = size(input(1).xt,2) - T + 1;
            error = zeros(dimx,N*T);
            
            for j=1:N
                gamma = hidden.gamma(1+(j-1)*T:j*T, :);
                simulated = self.simulate_onestep(gamma, {input(j)});
                error_j = input(j).xt(:,start:end) - simulated{1}(:,start:end);
                len = size(error_j,2);
                error(:,1+(j-1)*len : j*len) = error_j;
            end

        end
        
        
        %------------------------------------------------------------------
        function infoString = getInfoString(self, delim)
            infoString = sprintf(['This is the default infoString method in FEMParamEstResult\n',...
                                  'implement getInfoString(delim) method in', delim, class(self),'\n']);
        end
       
     end
    
    
    methods (Abstract)
        copy   = copy(self)
        parNum = getNumPar(self)
        resid  = compute_resid(self, input)
        simul  = simulate_onestep(self, hidden, input)
    end
end


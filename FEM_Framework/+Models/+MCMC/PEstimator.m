%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : PEstimator.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : MCMC estimate
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
classdef PEstimator < Core.PEstimator
    %PEstimator  MCMC METROPOLIS
    %   Detailed explanation goes here

    
    methods
        
        %------------------------------------------------------------------
        function self = PEstimator(cfg, input)
            self = self@Core.PEstimator(cfg, input);

        end
        
        
        %------------------------------------------------------------------
        function optRes = estimate(self, optRes)
            % Here is the MCMC algorithm starting, please note in case we
            % use shrinakeg methods, lasso or ridge regression, we need to
            % update the optRes.acf for the MCMC optimization but to
            % return the fem-acf = sum(sum(resid.*gamma))
            
            if isempty(optRes.param)
                aOptRes = self.estimateStartValue(optRes);
            else
                aOptRes = optRes.copy();
            end
            
            if aOptRes.acf < optRes.acf
                optRes = aOptRes;
                return;
            else
                noise = self.cfg.init_noise;
                beta  = self.cfg.init_beta;
                times_accepted = 0;
                [K, par_len] = size(aOptRes.param.theta);
                
                paramCov  = zeros(par_len, par_len, K);
                paramMean = zeros(par_len, K);
                
                for sample_step = 1:self.cfg.samples
                    pOptRes = self.proposeNext(aOptRes, paramCov, times_accepted, noise);
                    
                    if optRes.acf > pOptRes.acf
                        optRes = pOptRes;
                        return;                   
                    elseif self.accept(beta, self.shrinkageAcfUpdate(aOptRes.copy()), self.shrinkageAcfUpdate(pOptRes.copy()) )
                        aOptRes = pOptRes;
                        times_accepted = times_accepted + 1;
                        [paramCov, paramMean] = self.updateMeanCov(paramCov, paramMean, times_accepted, aOptRes);
                    end
                    
                    % adaptivity for MCMC
                    if sample_step >= 50
                        [noise, beta] = self.adapt(noise, beta, times_accepted, sample_step);
                    end
                end
                % we reached the end of the loop without valid result
                throw(MException('ResultChk:InvalidState','Can"t find better params,\n best %f',aOptRes.acf));
            end
        end
    end
    
    
    methods (Access = 'protected')
        %------------------------------------------------------------------
        function accepted = accept(~, beta, acceptedNRG, proposedNRG)
        %ACCEPT returns true if proposedNRG is accepted 
        %  accept in case proposedNRG < acceptedNRG.
        %  in case proposedNRG >= acceptedNRG accept with certain
        %  probability, that enables jumps into an other local minimas  
        %------------------------------------------------------------------
            accepted = false;
            acceptance_prob = min(1., exp(beta * (-proposedNRG + acceptedNRG)) );
            if acceptance_prob >= rand
                accepted = true;
            end
        end

        
        %------------------------------------------------------------------
        function [noise, beta] = adapt(self, noise, beta, counter, number_steps)
            % adjust noise
            if mod(number_steps,50)==0
                accept_rate = counter / number_steps;
                disp(['accept_rate ', num2str(accept_rate)]);
                
                if accept_rate > 0.7
                    noise = noise * self.cfg.noise_more;
                    disp(['noise increased, new noise: ',num2str(noise)])
                
                elseif accept_rate < 0.3
                    noise = noise * self.cfg.noise_less;
                    disp(['noise decreased, new noise: ',num2str(noise)])
                end
            end
            
            % adjust beta
            if mod(number_steps,100)==0
                beta = beta * self.cfg.beta_factor;
                disp(['increased beta: ', num2str(beta)]);
            end
        end
        
        
        %------------------------------------------------------------------
        function [paramCov, paramMean] = updateMeanCov(~, paramCov, paramMean, count, aOptRes)
            % this function is used to update the covMatrix wrt
            % [1] Optimal Proposal Distributions and Adaptive MCMC by
            % Jeffrey S. Rosenthal
            
            if count < 3
                return;
            end
            
            for i=1:aOptRes.hidden.K
                meanP = paramMean(:,i);
                covM = paramCov(:,:,i);
                
                % update mean
                meanP_new = ((count-1)*meanP + aOptRes.param.theta(i,:)') / count;
                
                % update covMatrix
                cov_new = (count-2)/(count-1) * covM + meanP*meanP'...
                    - count * (meanP_new*meanP_new') / (count-1) ...
                    + aOptRes.param.theta(i,:)'*aOptRes.param.theta(i,:) /(count-1);
                
                paramMean(:,i) = meanP_new;
                paramCov(:,:,i) = cov_new;
            end
        end
    end
    

    methods (Abstract, Access = 'protected')
        
        optRes = estimateStartValue(self, optRes)
        
        optRes = proposeNext(self, optRes, paramCov, count, noise)
        
        optRes = shrinkageAcfUpdate(self, optRes)

    end
    
    
end


%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : PEstimator.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : interface for estimate gpd-mcmc-cpp
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
classdef PEstimator < Models.MCMC.GPD.PEstimator
    %MODELS.MCMC.GPD.GPDCPP.PESTIMATOR GPD estimator
    %   This Version includes CPP implementation!
    
    methods
        
        %------------------------------------------------------------------
        function self = PEstimator(cfg, input)
            self@Models.MCMC.GPD.PEstimator(cfg, input);
        end
        
        
        %------------------------------------------------------------------
        function optRes = estimate(self, optRes)
            
            proposMaxTimes = 200;
            samples        = self.cfg.samples;
            noise          = self.cfg.init_noise;
            noise_more     = self.cfg.noise_more;
            noise_less     = self.cfg.noise_less;
            beta           = self.cfg.init_beta;
            beta_factor    = self.cfg.beta_factor;
            seed           = ceil(1 + (65535 -1)*rand);
            gaussFLAG      = self.cfg.gaussFLAG;
            lassoFLAG      = self.cfg.lassoFLAG;
            ridgeFLAG      = self.cfg.ridgeFLAG;
            
            if ~isempty(optRes.param)
                aOptRes = optRes.copy();
            else
                aOptRes = self.estimateStartValue(optRes);
            end
            
            
            try
                ut_ksi = cat(2, self.input.ut_ksi)';
                ut_si = cat(2, self.input.ut_si)';
                xt = cat(2,self.input.xt)';
                
                theta = Models.MCMC.GPD.GPDcpp.MEX_GPDestimate(aOptRes.hidden.gamma, aOptRes.resid,...
                    xt, ut_ksi, ut_si, aOptRes.param.theta, samples, proposMaxTimes, ...
                    noise, noise_more, noise_less, beta, beta_factor, seed,...
                    gaussFLAG, lassoFLAG, ridgeFLAG);
            catch excp
                disp(excp.getReport)
                warning('PEstimator:estimate','error in MEX_GPDestimate')
            end
            
            if ~isempty(theta)
                optRes.param = Models.MCMC.GPD.GPDcpp.Params(self.input(1).meta, theta);
                optRes.updateAcf(self.input);
                if optRes.acf < 0
                    error('PEstimator:estimate','error, negative acf!!!!')
                end
                
                % check acf convergency
                if (aOptRes.acf - optRes.acf) < 0 && (aOptRes.acf - optRes.acf) > -1e-09
                    warning('ResultChk:InvalidState','acf growing slightly')
                elseif (aOptRes.acf - optRes.acf) <= -1e-09
                    throw(MException('ResultChk:InvalidState', 'Growing ACF detected, acf %d before param opt., after %d', optRes.acf, aOptRes.acf));
                end
                
            else
                if ~isempty(optRes.param)
                    warning('PEstimator:estimate','empy parameters, return last valid optRes')
                else
                    warning('PEstimator:estimate','empy parameters, return StartValue Result')
                    optRes = aOptRes.copy();
                end
                
            end
            
            
            
            
        end
    end
    
    
end



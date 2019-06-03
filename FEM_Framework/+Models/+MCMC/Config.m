%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : Config.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : model configuration
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
classdef Config < Core.Construction.ModelConfig
    %CONFIG Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = 'public', SetAccess = 'public')
        subtype
        
        proposMaxTimes
        samples
        gaussFLAG
        eps
        init_beta
        beta_factor
        init_noise
        noise_more
        noise_less
        
        ut_combis
        
        lassoFLAG
        ridgeFLAG
        paramTHRESHOLD
    end
    
    methods
        function self = Config(type, subtype, proposMaxTimes, samples, lassoFLAG, ridgeFLAG,...
                gaussFLAG, eps, init_beta, beta_factor, init_noise, noise_more,...
                noise_less, ut_combis, paramTHRESHOLD)
            
            self = self@Core.Construction.ModelConfig(type);
            
            
            self.subtype = subtype;
            self.proposMaxTimes = proposMaxTimes;
            self.samples = samples;
            self.eps = eps;
            self.init_beta = init_beta;
            self.beta_factor = beta_factor;
            self.init_noise = init_noise;
            self.noise_more = noise_more;
            self.noise_less = noise_less;
            self.ut_combis = ut_combis;
            self.gaussFLAG = gaussFLAG;
            self.lassoFLAG = lassoFLAG;
            self.ridgeFLAG = ridgeFLAG;
            self.paramTHRESHOLD = paramTHRESHOLD;
            if (ridgeFLAG~=0 && lassoFLAG~=0)
                error('Lasso and Shirnkage active, please choose only one')
            end
        end
    end
end


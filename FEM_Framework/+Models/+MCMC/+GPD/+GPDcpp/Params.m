%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : Params.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : interface for acf / resid gpd-mcmc-cpp
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
classdef Params < Models.MCMC.GPD.Params
    %PARAMS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function self = Params(meta, theta)
            self@Models.MCMC.GPD.Params(meta, theta) 
        end
        
        
        %------------------------------------------------------------------
        function acf = computeAcf(s, hidden, resid)
            acf = Models.MCMC.GPD.GPDcpp.MEX_GPDcomputeAcf(hidden.gamma, resid);
            
            % shrinkage update
            ksi = s.theta(:,s.ksi_mask);
            ksi = ksi(:,2:end);
            si  = s.theta(:,s.si_mask);
            si  = si(:,2:end);
            
            if s.meta.lassoFLAG ~= 0 
                acf = acf + s.meta.lassoFLAG * (norm(ksi(:),1) + norm(si(:),1));
            elseif s.meta.ridgeFLAG ~= 0 
                acf = acf + s.meta.ridgeFLAG * (norm(ksi(:)) + norm(si(:)))^2;
            end
        end
    
        %------------------------------------------------------------------
        function resid = compute_resid(self, input)
            
            %TODO check input (length, dim, no ones,..)
			ut_ksi = cat(2, input.ut_ksi)';
            ut_si = cat(2, input.ut_si)';
            xt = cat(2,input.xt)';
                        
            resid = Models.MCMC.GPD.GPDcpp.MEX_GPDresid(self.theta, xt, ut_ksi, ut_si);
        end
       
    end
end

%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : DGamma.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : cfg for discrete gamma
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
classdef DGamma < Core.Datatypes.CGamma
    %DGAMMA Holds Gamma estimation result for the descrite case
    %   The computation of the number of parameters in Gamma for AIC 
    %   computation is different in diccrete case.
    
    methods
        function self = DGamma(gamma, T, K, Ct)
            self@Core.Datatypes.CGamma(gamma, T, K, Ct)
        end
        
        %------------------------------------------------------------------
        function parNum = getNumPar(s)
        %GETNUMPAR overwriten method from CGamma 
        %   in diccrete case have to compute in different way
        
            if s.K == 2
                parNum = s.nbins;
            else
                parNum = 2 * s.nbins - 1;
            end 
        end
        
        %------------------------------------------------------------------
        function copy = copy(self)
            copy = Core.Datatypes.DGamma(self.gamma, self.T, self.K, self.Ct);
            copy.femModelParameter = self.femModelParameter;
        end
        
    end
end


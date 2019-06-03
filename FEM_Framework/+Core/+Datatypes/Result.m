%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : Result.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : estiamte IC (AIC/BIC/AICc) for a result
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
classdef Result < handle
    %RESULT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties 

        param               % data model parameter (<-FEMParamEstimator.estimate())
        hidden              % gamma data and model holder 
        resid               % (TN x K) martix
        acf                 % number
        
        ic                  % struct
        
    end
    
    methods  
        
        %------------------------------------------------------------------
        function self = Result(modelParams, hidden, acf, resid)
            self.param = modelParams;
            self.hidden = hidden;
            self.acf = acf;
            self.resid = resid;
            self.ic = []; 
        end
        
        %------------------------------------------------------------------
        function copy = copy(self)
            
            param_copy = [];
            if ~isempty(self.param)
                param_copy = self.param.copy();
            end
            
            hidden_copy = [];
            if ~isempty(self.hidden)
                hidden_copy = self.hidden.copy();
            end
            
            copy = Core.Datatypes.Result(param_copy, hidden_copy, self.acf, self.resid);
            copy.ic = self.ic;
        end
        
        %------------------------------------------------------------------
        function compute_IC(self, inp)

            % this function is used to compute IC for ||X_t - Theta_t||_2^2 
            % types of problems
            T = sum(self.hidden.Tlist);
            K = self.hidden.K;
            
            % get Gaussian Log Like
            [ logLike, numLLParam ] = self.param.compute_loglike(self.hidden, inp);
            
            % estimate number of parameters for the model
            hidParNum = self.hidden.getNumPar();
            modParNum = self.param.getNumPar();
            
            % estimate all number of params
            number_param = hidParNum + modParNum + numLLParam;
            
            % before computing IC, check if overfit occures in one of clusters
            time_in_cluster = sum(self.hidden.gamma,1);

            if all(time_in_cluster > modParNum/K) && T > (number_param  + 1) 
                
                % compute AIC 
                self.ic.AIC = -2*logLike + 2*number_param; 

                self.ic.AICc = self.ic.AIC + 2*number_param*(number_param+1) / (T - number_param - 1);

                % compute BIC
                self.ic.BIC = -2*logLike + log(T) * number_param;
            else
                disp('Overfitting, force IC = inf')
                self.ic.AIC = Inf ;
                self.ic.AICc = Inf;
                self.ic.BIC = Inf;
            end      
        end
         
        %------------------------------------------------------------------
        function infoString = getInfoString(self)
            for i = 1:1
                AIC  = 'not computed';
                AICc = 'not computed';
                BIC  = 'not computed';
                delim ='\t';

                if ~isempty(self.ic) 
                    if ~isempty(self.ic.AIC)
                        AIC = num2str(self.ic.AIC);
                    end
                    if ~isempty(self.ic.AICc)
                        AICc = num2str(self.ic.AICc);
                    end
                    if ~isempty(self.ic.BIC)
                        BIC = num2str(self.ic.BIC);
                    end
                end
            end
            
            infoString = sprintf(['acf = ', num2str(self.acf), delim,...
                                  'AIC = ', AIC, delim,...
                                  'AICc = ', AICc, delim,...
                                  'BIC = ', BIC, '\n',...
                                  'Hidden info : ',delim, class(self.hidden), delim, ...
                                  self.hidden.getInfoString(delim),'\n',...
                                  'Params info : ', delim, class(self.param), delim,...
                                  self.param.getInfoString(delim), '\n']);
                              
        end
        
        %------------------------------------------------------------------
        function updateAcf(self, inp)
            self.resid = self.param.compute_resid(inp);
            self.acf = self.param.computeAcf(self.hidden, self.resid);
        end
        
    end
    
end


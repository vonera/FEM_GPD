%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : Params.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : acf, resid
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
classdef Params < Models.MCMC.Params
    %PARAMS Summary of this class goes here
    %   Detailed explanation goes here
    properties
        ksi_mask
        si_mask
    end
    
    methods
        %------------------------------------------------------------------
        function self = Params(meta, theta)
            self@Models.MCMC.Params(meta, theta)
            
            ksi_l = size(self.meta.umask_ksi,2)+1;
            si_l = size(self.meta.umask_si,2)+1;
                
            self.ksi_mask = [true(1,ksi_l)  false(1,si_l) ]; 
            self.si_mask =  [false(1,ksi_l) true(1,si_l)  ];  
        end
        
        

        
        %------------------------------------------------------------------
        % overwrite this function defined in Core.Datatypes.Params
        function acf = computeAcf(s, hidden, resid)

            acf = sum(sum(resid .* hidden.gamma,1) ,2);

            % Lasso regression only on the covas coefficients, off set is
            % ignored
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
        function copy = copy(self)
            copy = Models.MCMC.GPD.Params(self.meta, self.theta);
        end
        
        %------------------------------------------------------------------
        function resid = compute_resid(self, input)

            if ~isequal(input(1).meta, self.meta)
                error('Invalid input for this parameters')
            end
                       
            Tlist = input(1).meta.Tlist;
            Tlist_cs = cumsum([0 Tlist]);
            K = size(self.theta, 1);
            
            % initialize negLogLike with big values;
            resid = 3e+10*ones(sum(Tlist), K);
            
            for i = 1:K
                [ksi_i, si_i] = self.getForCluster(i);
                
                for j=1:size(input, 2)
                    
                    for t = 1:Tlist(j)        
                        tj = t + Tlist_cs(j);
                        ksi = ksi_i * input(j).ut_ksi(:,t );
                        si  = si_i * input(j).ut_si(:,t);

                        if ksi ~= 0
                            zt = 1 + ksi * input(j).xt(:,t) / si;
                            if zt > 0
                                resid(tj,i) =  log(si) + (1+1/ksi)*log(zt);
                            end
                        else
                            resid(tj,i) = log(si) + input(j).xt(:,t) / si;
                        end
                    end
                end
            end
        end
        
        
        %------------------------------------------------------------------
        function [ksi_i, si_i] = getForCluster(self, i)
            
            ksi_i = self.theta(i, self.ksi_mask);
            si_i = self.theta(i, self.si_mask);
        end
        
        
        %------------------------------------------------------------------
        function updateAcf(self, inp)
            self.resid = self.compute_resid(inp);
            self.acf = self.computeAcf(self.hidden, self.resid);
        end
    end
end

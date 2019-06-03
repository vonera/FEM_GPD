%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : ModelFit.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : alternating optimization
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
classdef ModelFit

    properties
        cfg
    end
    
    methods
        
        %------------------------------------------------------------------
        function self = ModelFit(cfg)
            self.cfg = cfg;
        end
       
        %------------------------------------------------------------------
        function rHandler = fit(self, iHandl)
            
            rHandler = Core.Construction.Factory.createResultHandler(iHandl.cfg);
            
            for i = 1:iHandl.getNumOfInputs()
                
                input = iHandl.get(i);
                
                disp(['### Start optimization for input: ', input(1).meta.getInfoString()]);
                
                for K = self.cfg.Klist

                    optRes = [];

                    for Ct = self.cfg.Ctlist                                       
                        
                        
                        for Cn = self.cfg.Cnlist 
                            paramEst = Core.Construction.Factory.createParamEstimator(iHandl.cfg.model, input);

                            gammaEst = Core.Construction.Factory.createGammaEstimator(iHandl.cfg, input, K, Ct, Cn);

                            disp(['### Start annealing for K=' num2str(K), ' ,Ct=', num2str(Ct), ' ,Cn=', num2str(Cn)]);

                            anneallings = self.cfg.annealing;
                            if K == 1 && ~self.cfg.annealingForK1
                                anneallings = 1;
                            end

                            optRes = self.parallelAnnealing(optRes, paramEst, gammaEst, anneallings);                    

                            % rHandler stores a copy of optRes  
                            rHandler.add(optRes);

                            if K == 1
                                break
                            end
                        end
                        if K == 1
                            break
                        end
                    end
                end
            end
        end
        %------------------------------------------------------------------
    end
    
    
    methods (Access = 'protected')
        
        %------------------------------------------------------------------
        function bestAnneal = parallelAnnealing(self, prevRes, paramEst, gammaEst, num_annealings)
            
            annealings = cell(1,num_annealings);
            
            for a = 1:num_annealings
                if a == 1 && ~isempty(prevRes)
                    try
                        boot = gammaEst.estimate(prevRes.copy());
                    catch ex
                        if strcmp(regexp(ex.identifier, '(?<=:)\w+$', 'match'), 'InvalidState')
                            fprintf(['>> extended bootstrapping of ',...
                                'gamma failed: invalid solver state ',...
                                'detected (', ex.message,'), continue ',...
                                'with simple bootstrapping\n']);
                        else
                            ex.rethrow();
                        end
                        boot = prevRes.copy(); % only this is also ok! 
                    end
                else
                    boot = Core.Datatypes.Result([], gammaEst.generateRandomGamma(), inf, []);
                end
                annealings{a} = self.subspace(boot, paramEst, gammaEst,a);
            end
            
            %find best annealing model
            bestAnneal = annealings{1};
            for a = 2:num_annealings
                
                optResult = annealings{a};
                
                if isempty(bestAnneal)
                    bestAnneal = optResult;
                    
                elseif ~isempty(optResult) && optResult.acf < bestAnneal.acf
                    bestAnneal = optResult;
                
                end
            end
            disp(['>> Best annealing result with acf ', num2str(bestAnneal.acf)]);
        end
                
        
        %------------------------------------------------------------------
        function optRes = subspace(self, boot, paramEst, gammaEst, a)
            
            optRes = boot; 
            prev_acf = inf;

            try
                for s = 1:self.cfg.subspaceItr

                    optRes = paramEst.estimate(optRes.copy());
                    
                    if prev_acf - optRes.acf < self.cfg.subspaceEps
                        disp(['>> annealing step ', num2str(a),...
                            ', subspace convergence after ' num2str(s),...
                            ' steps by ', num2str(optRes.acf)]);
                        return;
                    end
                    
                    optRes = gammaEst.estimate(optRes.copy());

                    prev_acf = optRes.acf;
                    
                end
                disp(['>> annealing step ', num2str(a),', subspace iterations exceeded by ', num2str(optRes.acf)]);
                warning('FEM:subspace','Number of subspace iterations exceeded, cfg.modelFit.subspaceItr');
            catch ex
                % Get last segment of the error message identifier.
                if strcmp(regexp(ex.identifier, '(?<=:)\w+$', 'match'),...
                        'InvalidState')
                    fprintf(['>> annealing step ', num2str(a),...
                        ', subspace step ',int2str(s),...
                        ', invalid solver state detected (',...
                        ex.message,'), returning with last valid result ',...
                        num2str(optRes.acf),'\n']);
                else
                    ex.rethrow();
                end
            end
        end
        %------------------------------------------------------------------
    end
    
end

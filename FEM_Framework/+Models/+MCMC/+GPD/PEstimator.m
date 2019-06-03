%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : PEstimator.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : MCMC algo
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
classdef PEstimator < Models.MCMC.PEstimator
    %MODELS.MCMC.GPD.PESTIMATOR MCMCMetGEV estimator. EstimateStartvalue
    %computes ksi and si as samples from the corresponding convexhull wrt
    %ut_ksi and ut_si. ProposeNext() uses the covariance matrix for 
    %estimation of the next possible parameter set. This was
    %done wrt [1] Optimal Proposal Distributions and Adaptive MCMC
    %by Jeffrey S. Rosenthal, [2] Adaptive Optimal Scaling of 
    %Metropolis-Hastings Algorithms Using the Robbins-Monro Process, 
    %P. H. GARTHWAITE & Y. FAN
    %#############################   
    
    methods
        
        %------------------------------------------------------------------
        function self = PEstimator(cfg, input)
            self@Models.MCMC.PEstimator(cfg, input);

        end
    end
    
    
    methods (Access = 'protected')
        %------------------------------------------------------------------
        function sOptRes = estimateStartValue(self, optRes)
            % TODO make sigma max value dependent on data
            
            K = optRes.hidden.K;
            theta = zeros(K, self.input(1).d_ksi + self.input(1).d_si);
            maxTrys = self.cfg.proposMaxTimes;
            convhull_ksi = [min(cat(2,self.input.ut_ksi),[],2) max(cat(2,self.input.ut_ksi),[],2)]';
            convhull_si = [min(cat(2,self.input.ut_si),[],2) max(cat(2,self.input.ut_si),[],2)]';
			xt_max = abs( max( cat(2,self.input.xt) )' );
            
            for i=1:K             
                iNext = 0;
                flag = 0;
                
                while ~flag && iNext < maxTrys
                    ksi_i = self.sampleFromConvHull(convhull_ksi, [-0.49, 0.49], [])';
                    si_i = self.sampleFromConvHull(convhull_si, [0.0001, xt_max/4], [])';
                    
                    flag = self.proveConstraints(optRes, ksi_i, si_i, i);
                    iNext = iNext+1;
                end
                
                if flag
                    theta(i,:) = [ksi_i, si_i];
                else
                    throw(MException('ResultChk:InvalidState','Can not find StartValue'));
                end
            end
                       
            sOptRes = optRes.copy();
            sOptRes.param = Models.MCMC.GPD.Params(self.input(1).meta, theta);
            sOptRes.updateAcf(self.input);           
        end
        

        %------------------------------------------------------------------
        function pOptRes = proposeNext(self, optRes, paramCov, count, noise)
            % generate next param, if constraints are not fullfilled
            % try till maxNext is reached
            
            [K, len_param] = size(optRes.param.theta);
            cov_id = 0.01/(len_param)*eye(len_param);
            theta = zeros(K,len_param);
            covMatrix = cov_id;
            flag = 0;
            maxTrys = self.cfg.proposMaxTimes;
            
            for i = 1:K
                if ( count > 2*len_param && self.cfg.gaussFLAG )
                    covMatrix = 5.6644/len_param * paramCov(:,:,i);
                end
                              
                iNext = 0;
                flag = 0;
                
                % propose next for cluster i
                while ~flag && iNext < maxTrys
                    theta(i,:) = (1 - noise) * mvnrnd(optRes.param.theta(i,:), covMatrix)...
                                     + noise * mvnrnd(optRes.param.theta(i,:), cov_id);
                    
                    ksi_i = theta(i, optRes.param.ksi_mask);
                    si_i = theta(i, optRes.param.si_mask);
                    
                    flag = self.proveConstraints(optRes, ksi_i, si_i, i);
                    iNext = iNext + 1;
                end
                
                % no parameters for cluster i found
                if ~flag
                    break
                end
            end
            
            % if no appropriate parameters were proposed return old
            if ~flag
                pOptRes = optRes;
            else
                pOptRes = optRes.copy();
                pOptRes.param = Models.MCMC.GPD.Params(self.input(1).meta, theta);
                pOptRes.updateAcf(self.input);
            end
        end
        
          
        %------------------------------------------------------------------
        function shrinkACF = shrinkageAcfUpdate(self, optRes)
            % this function updates the acf in case shrinkage (lasso/ridge)
            % is used
            theta = optRes.param.theta;
            
            if self.cfg.ridgeFLAG ~= 0
                for i = 1:optRes.hidden.K
                    ksi = theta(i, optRes.param.ksi_mask);
                    si  = theta(i, optRes.param.si_mask);
                    param_norm = sum(sum(ksi(2:end).^2)) + sum(sum(si(2:end).^2));
                    optRes.acf = optRes.acf + self.cfg.ridgeFLAG*param_norm;
                end
            elseif self.cfg.lassoFLAG ~= 0
                for i = 1:optRes.hidden.K
                    ksi = theta(i, optRes.param.ksi_mask);
                    si  = theta(i, optRes.param.si_mask);
                    param_norm = sum(abs(ksi(2:end))) + sum(abs(si(2:end)));
                    optRes.acf = optRes.acf + self.cfg.lassoFLAG*param_norm;
                end
            end
            
            shrinkACF  = optRes.acf;
        end
        
        
        %------------------------------------------------------------------
        function flag = proveConstraints(self, optRes, ksi_i, si_i, cl)

            flag = true;
            Tlist_cs = cumsum([0 optRes.hidden.Tlist]);
            
            for j = 1:optRes.hidden.N
                
                % check if:  -0.5 < ksi(ut) < 0.5
                ksi_vec = ksi_i*self.input(j).ut_ksi;
                if any(abs(ksi_vec) > 0.5)
                    flag = false;
                    return;
                end
                
                % check if: si(ut) > 0
                si_vec  = si_i*self.input(j).ut_si;
                if any(si_vec < 0)
                    flag = false;
                    return;
                end
                
                % check if: ( 1 + ksi(ut)*xt / si(ut) ) > 0 for each j
                indexj = optRes.hidden.gamma( 1+Tlist_cs(j) : Tlist_cs(j+1), cl) > 0;
                zt = 1 + (ksi_vec(indexj) .* self.input(j).xt(indexj)) ./ si_vec(indexj);
                if min(zt) < 0
                    flag = false;
                    return;
                end
            end
            
        end
        
      
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function xsample = sampleFromConvHull(self, convhull, boundary, xstart)
            % sampling from a convex hull, such that constraints on x are
            % fulfilled: lb <= Ax <= ub. Rewritting them we obtain:
            % -Ax <= -lb |
            %            | [-A; A] x <= [-lb; ub]
            %  Ax <=  ub |
            A = self.bichoose( convhull );
            
            lb = boundary(1)*ones(size(A,1),1);
            ub = boundary(2)*ones(size(A,1),1);
            
            if isempty(xstart)
                xstart = zeros(size(A,2),1);
            end
            
            % check if stratvalue is feasible: Ax <= b, else get a random
            rndMax = 500;
            iNext = 1;
            while max( ([-A; A]*xstart > [-lb; ub]) ) && iNext < rndMax
                xstart = A\lb + 0.1*A\(ub-lb).*rand(size(A,2),1) ;
                iNext = iNext + 1;
            end
            
            A = [-A; A];
            b = [-lb; ub];
            
            xsample = Models.MCMC.Tools.ConvexSampler.convexSampler(A, b, xstart);
        end
          
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function [ pNtimesN ] = bichoose(~, twoTimesN )
            %(c) 2012, Dimitri Igdalov
            %BICHOOSE
            %   from    xda creates xda
            %           yeb         xdb
            %                       xea
            %                       xeb
            %                       yda
            %                       ydb
            %                       yea
            %                       yeb
            
            N = size(twoTimesN,2);
            pN = 2^N;
            pNtimesN = zeros(pN,N);
            len = pN;
            
            for i = 1:N
                pair = [twoTimesN(1,i)*ones(len/2,1) ; twoTimesN(2,i)*ones(len/2,1)];
                
                for begin = 1:len:pN
                    pNtimesN(begin:begin+len-1,i) = pair;
                end
                
                len = len / 2;
            end
        end
    end
    
    
end


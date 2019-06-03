%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : Inference_with_GammaInterp.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : generate test-example
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
classdef DataGen
    %DATAGEN
    %
    
    methods (Static)
        %------------------------------------------------------------------
        function  input = generate(Gamma, input, Pcell)
            check = true;
            input = Models.MCMC.GPD.DataGen.generate_checked(Gamma, input, Pcell, threshold, Tlist, check);
        end
        
        
        %------------------------------------------------------------------
        function input = generate_checked(Gamma, input, Pcell, threshold, Tlist, check)

            %   in: Tlist <- list of T for each location
            %   in: Gamma
            %   in: input <- contains for each locations the Ut
            %   in: Pcell
            %   out: input <- contains for each location xt and ut
            
            [~, K] = size(Gamma);
            Tlist_cs = cumsum([0 Tlist]);
            
            for j = 1:length(Tlist)
                if check
                    Gammaj = Gamma(1+Tlist_cs(j):Tlist_cs(j+1),:);
                    Models.MCMC.GPD.DataGen.check_input(Gammaj, input(j).ut, input(j).meta, Pcell);
                end
                
                dimksi = 1 + length(input(j).meta.umask_ksi);
                dimsi = 1 + length(input(j).meta.umask_si);
                
                Tj = Tlist(j);
                Ut_ksi = [ones(1,Tj); input(j).ut(input(j).meta.umask_ksi, :)];
                Ut_si  = [ones(1,Tj); input(j).ut(input(j).meta.umask_si, :)];
                
                for t = 1:Tj
                    xt = 0;
                    
                    for i = 1:K
                        ksi = Pcell{i,1}(1,1:dimksi)*Ut_ksi(:,t);
                        si  = Pcell{i,1}(1,dimksi+1:dimksi+dimsi)*Ut_si(:,t);                       
                        xt = xt + Gammaj(t, i) * gprnd(ksi,si,threshold(j),1,1);
                    end
                    
                    input(j).xt(t) = xt;
                end
            end
        end
        
        
        %------------------------------------------------------------------
        % checkers
        function check_input(Gamma, Ut, meta, Pcell)
            
            
            % check Gamma is correct
            if ~Models.MCMC.GPD.DataGen.check_is_gamma(Gamma)
                error('Gamma is not a correct affiliation vector');
            end
            
            % check Gamma and Ut length
%             if ~Models.MCMC.GPD.DataGen.check_length(Gamma, Ut)
%                 error('Length of Gamma and Ut does not match');
%             end
            
            % warn if Ut has added ones in the first row
            if isequal(Ut(1,:), ones(1,size(Ut,2)))
                warning('Ut(1,:) has only ones'); % TODO better handling
            end
            
            % check number of clusters
            if ~Models.MCMC.GPD.DataGen.check_num_clusters(Pcell, Gamma)
                error('Number of clusters defined by Gamma does not match the number defined by Pcell');
            end
            
            % check number of externals
%             if ~Models.MCMC.GPD.DataGen.check_num_external(Pcell, Ut)
%                 error('The number of influence defined by Pcell does not match the dimension of Ut');
%             end
            
            
            % check constraints for GPD parameters
            if ~Models.MCMC.GPD.DataGen.check_param_const(Pcell, Ut, meta)
                error('GPD parameters do not fulfill constraints');
            end
            
        end
        
        
        %------------------------------------------------------------------
        function [ ret ] = check_param_const(Pcell, Ut, meta)
            ret = true;
                
            T = size(Ut,2);
            K = size(Pcell, 1);
            
            dimksi = 1 + length(meta.umask_ksi);
            dimsi = 1 + length(meta.umask_si);
            
            Ut_ksi = [ones(1,T); Ut(meta.umask_ksi, :)];
            Ut_si  = [ones(1,T); Ut(meta.umask_si, :)];
            
            for i = 1:K
                
                ksi = Pcell{i,1}(1,1:dimksi)*Ut_ksi;
                si  = Pcell{i,1}(1,dimksi+1:dimksi+dimsi)*Ut_si;

                if max( abs( ksi ) ) > 0.5
                    disp('Shape parameter ksi must be element of(-0.5, 0.5)');
                    ret = false;
                end
                
                if min( si ) <= 0
                    disp('Scale parameter si must be > 0');
                    ret = false;
                end
            end
        end
        
        
        %------------------------------------------------------------------
        function [ ret ] = check_length(Gamma, Ut)
            ret = true;
            if size(Gamma,1) ~= size(Ut,2)
                ret = false;
            end
        end
        
        
        %------------------------------------------------------------------
        function [ ret ] = check_num_external(Pcell,Ut)
            ret = true;
            K = size(Pcell,1);
            for i=1:K
                if size(Pcell{i,1},2) - 3 ~= size(Ut,1)
                    ret = false;
                    break
                end
            end
        end
        
        
        %------------------------------------------------------------------
        function [ ret ] = check_num_clusters(Pcell, Gamma)
            ret = true;
            if size(Pcell,1) ~= size(Gamma,2)
                ret = false;
            end
        end
        
        
        %------------------------------------------------------------------
        function [ ret ] = check_is_gamma(Gamma)
            ret = true;
            % check Gamma all elements of Gamma sum up to one
            
            if max(abs(sum(Gamma,2) - 1)) > eps
                disp('column Sum of Gamma is not equal 1');
                ret = false;
            end
            
            % check all elements of Gamma are not negative
            if min(min(Gamma)) < 0
                disp('Gamma has negative elements');
                ret = false;
            end
            
            % check all elements of Gamma are not greater than 1
            if max(max(Gamma)) > 1
                disp('Gamma has elements greater than 1');
                ret = false;
            end
            
        end
        
        
        %------------------------------------------------------------------
        function Pmatrix = pcell2matrix(Pcell)
            K = size(Pcell,1);
            len_matrix = size(Pcell{1,1},2);
            Pmatrix = zeros(K,len_matrix);
            
            for i=1:K
                Pmatrix(i,:) = Pcell{i,1};
            end
            
        end
        
        
        %------------------------------------------------------------------
        function acf_exact = computeAcf(input, GammaN, Pcell)
            
            % by mcmc resid is negLogLike
            N = size(input, 2);
            Ttotal = size(GammaN,1);
            K = size(self.theta, 1);
            dimksi = size(input(1).ut_ksi, 1);
            dimsi = size(input(1).ut_si, 1);
            
            % initialize negLogLike with big values;
            resid = 3e+10*ones(Ttotal, K);
            
            for i = 1:K
                
                ksi_i = Pcell{i,1}(1,1:dimksi);
                si_i  = Pcell{i,1}(1,dimksi+1 : dimski+ dimsi);
                
                for j=1:N
                    xtj = input(j).xt - threshold;
                    
                    for t = 1:T
                        tj = t + (j-1)*T;
                        ksi = ksi_i * input(j).ut_ksi(:,t );
                        si  = si_i * input(j).ut_si(:,t);
                        
                        if ksi ~= 0
                            zt = 1 + ksi * xtj(t) / si;
                            if zt > 0
                                resid(tj,i) =  log(si) + (1+1/ksi)*log(zt);
                            end
                        else
                            resid(tj,i) = log(si) + xtj(t) / si;
                        end
                    end
                end
            end
            
            acf_exact = sum(sum(GammaN.*resid,1),2);    
        end
        
    end
    
end

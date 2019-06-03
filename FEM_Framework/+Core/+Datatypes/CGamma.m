%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : CGamma.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : cfg for continous gamma
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
classdef CGamma < handle
    %CGAMMA Summary of this class goes here
    %  !!!ATTENTION!!! 
    %     if you derive from this class you have to overwrite
    %     the copy() method! If not it will always return objects of type 
    %     CGamma and not of the subclass (No copy constructors in matlab)
    
    properties (GetAccess = 'public', SetAccess = 'protected') % spatial
        discrete     % scalar
        Cn           % scalar
        Ct           % scalar
        N            % scalar
        Tlist        % list of scalars
        K            % scalar

        nbins        % scalar
        gamma        % K*NT 
    end
    
    
    properties 
        femModelParameter  % gamma model holder (are conkret FEMmodelParameterMarkov) 
    end
    
    
    methods
        %------------------------------------------------------------------
        function self = CGamma(GammaN, Tlist, Ct, Cn, discrete) % spatial
            [NT, self.K] = size(GammaN);
            if ~isequal(NT,sum(Tlist))
                error('Bad size of gamma: Modulus after division of len(gamma) by N is not Zero');
            end
            self.discrete = discrete;
            self.Cn = Cn;
            self.Ct = Ct;
            self.N = size(Tlist,2);
            self.gamma = GammaN;
            self.Tlist = Tlist;
            self.nbins = [];
        end
        
        %------------------------------------------------------------------
        function nbins = get.nbins(self)
            % Getter for nbins
            % Getters are very slow, but this is not relevant in this
            % particular case because this one is only called by
            % post processing and also only one or two times. 
            if isempty(self.nbins)
                self.nbins = self.computeNbins(self.gamma, self.Tlist);
            end
            nbins = self.nbins;
        end

        %------------------------------------------------------------------
        function parNum = getNumPar(s) % spatial
            % for more info check the document providet by Olga Kaiser
            % parNum = s.nbins * s.K;
            % parNum = s.nbins * (s.K-1);                                                                                                                                                                                                                    
            % parNum = (s.nbins * s.K) - 1; % same as nbins * (K-1) + (nbins - 1)                                                                                                                                                                                           
            
            if s.K == 1
               parNum = 0; 
               return
            end
            
            if s.discrete
                disp('num of parameters for discrete Gamma')
                if s.K == 2
                    parNum = s.nbins;
                else
                    parNum = 2 * s.nbins - 1*s.N;
                end
            else
                disp('num of parameters for continuous Gamma')
                parNum = (s.nbins * s.K) - 1*s.N; % same as nbins * (K-1) + (nbins - 1)
            end
        end
        
        %------------------------------------------------------------------
        function copy = copy(self)  % spatial
            % !!!ATTENTION!!! if you derive from this class you have to overwrite
            % this method!!! if not it will always return objects of type 
            % CGamma (No copy constructors in matlab)
            copy = Core.Datatypes.CGamma(self.gamma, self.Tlist, self.Ct, self.Cn, self.discrete);
            copy.femModelParameter = self.femModelParameter;
        end
                
        %------------------------------------------------------------------
        function infoString = getInfoString(self,delim) % spatial
            Ct_ = self.computeCt(self.gamma, self.Tlist);
            Cn_ = self.computeCn(self.gamma, self.Tlist);
            infoString = sprintf(['K = ',num2str(self.K), delim,...
                                  'N = ',num2str(self.N), delim,... 
                                  'Ct = ',num2str(self.Ct), delim,... 
                                  'Cn = ',num2str(self.Cn),... 
                                  ', Ct_real(max over N) = ',num2str(Ct_),... 
                                  ', Cn_real(max over T) = ',num2str(Cn_),... 
                                  ', bins = ', num2str(self.nbins)]);
        end
               
    end
    
    
    methods (Static)
    
        %------------------------------------------------------------------
        function hidden = fromGamma(GammaN, Tlist) % spatial
            NT = size(GammaN,1);
            if mod(NT,sum(Tlist))
                error('Bad size of gamma: Modulus after division of len(gamma) by N is not Zero');
            end
            [Ct_,Cn_] = Core.Datatypes.CGamma.computeC(GammaN, Tlist);
            discr = Core.Datatypes.CGamma.checkIsGamma(GammaN);
            hidden = Core.Datatypes.CGamma(GammaN, Tlist, Ct_, Cn_, discr);
        end
        
        %------------------------------------------------------------------
        function bool = checkIsDiscrete(GammaN) % spatial
            % like that GammaN(GammaN<1)>0 but with tolerance my_eps 
            my_eps = 1e-14;
            bool = ~any(GammaN(abs(GammaN-1) > my_eps) > my_eps);
        end
        
        %------------------------------------------------------------------
        function bool = checkIsGamma(GammaN) % spatial
            my_eps = 1e-14;
            bool = ~any(abs(sum(GammaN,2)-1) > my_eps);
        end
        
        %------------------------------------------------------------------
        function hidden = getStatGamma(Tlist) % spatial
            GammaN = ones(sum(Tlist),1);
            hidden = Core.Datatypes.CGamma(GammaN, Tlist, 0, 0, true);
        end
        
        
        %------------------------------------------------------------------
        function nbins = computeNbins(GammaN, Tlist) % spatial
            
            NT = size(GammaN,1);
            if mod(NT,sum(Tlist))
                error('Bad size of gamma: Modulus after division of len(gamma) by N is not Zero');
            end
            
            Tlist_cs = cumsum([0 Tlist]);
            nbins = 0;
            for j = 1:length(Tlist)
                Gamma = GammaN(1+Tlist_cs(j):Tlist_cs(j+1),:);
                nbins_j = 1;
                
                for t = 1:Tlist(j)-1
                    
                    if ~isequal(Gamma(t,:),Gamma(t+1,:))
                        nbins_j = nbins_j + 1;
                    end
                end
                
                nbins = nbins + nbins_j;
            end
        end
        
        %------------------------------------------------------------------
        function Ct= computeCt(GammaN, Tlist) % spatial
            % computes max Ct along all N locations.
            Tlist_cs = cumsum([0 Tlist]);
            Ct = 0;
            
            for j = 1:length(Tlist)
                Gamma = GammaN(1+Tlist_cs(j):Tlist_cs(j+1),:)';
                Ct_max = max(max(sum(abs(diff(Gamma,1,2)),2)));
                
                if Ct < Ct_max
                    Ct = Ct_max;
                end
            end
            
        end
        
        
        %------------------------------------------------------------------
        function Cn = computeCn(GammaN, Tlist) % spatial
            % computes max Cn along all T timesteps.
 
            %warning('Is not completly clear how to compute the Cn')
            Cn = NaN;
        end
        
        
        %------------------------------------------------------------------
        function allG = getPlotableStat(Gamma)
            K_ = size(Gamma,2);
            
            allG = 0;
            for i = 1:K_
                allG = allG + Gamma(:,i) * i;
            end
            allG = allG/K_;
        end
        
    end
    
end


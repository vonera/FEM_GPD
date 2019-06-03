%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : GEstimator.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : matrices for BV
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
classdef GEstimator < Core.Gamma.Spatial.GEstimator
    %GESTIMATORSBV Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        N_1
        NK
        tBinsN_1
        tBinsK
        tBins_1NK
        tBinsN_1K
        tBinsNtmp
        
        Aeq
        beq
        Aneq
        bneq
    end
    
    methods
        
        %------------------------------------------------------------------
        function self = GEstimator(cfg, input, K, Ct, Cn) % gamma cgf
            self = self@Core.Gamma.Spatial.GEstimator(cfg, input, K, Ct, Cn);
            
            self.N_1 = self.N-1;
            self.NK = self.N * self.K;
            self.tBins_1NK = sum(self.tBinsList-1) * self.K;
                        
            if all( self.tBinsList == self.tBinsList(1) ) && ~isinf(Cn)
                self.tBinsNtmp = self.tBins*0.5*(self.N^2 - self.N);
                self.tBinsN_1 = self.tBins * self.N_1;
                self.tBinsK = self.tBins * self.K;
                self.tBins_1NK = (self.tBins-1) * self.NK;
                self.tBinsN_1K = self.tBins * self.N_1 * self.K;
            else
                self.tBinsNtmp = 0;
                self.tBinsN_1 = 0;
                self.tBinsK = 0;
                self.tBinsN_1K = 0;
            end
            
            
            [self.Aeq, self.beq] = self.getAeqBeq();
            self.Aneq = self.getAneq();
            self.bneq = self.getBneq();
        end
        
    end
    
    methods (Abstract)
        
        hidden = estimate(self, optRes)
        
        hidden = generateRandomGamma(self)
        
    end
    
    
    methods (Access = 'protected')
        
        %------------------------------------------------------------------
        function res_vec = constructRes_vec(self, resid)
            % from (N Residua in row)' Resid el R^{TNxK}
            fe_resid = zeros(self.tBinsN, self.K);
            
            for b = 1:self.tBinsN
                fe_resid(b,:) = sum(resid( self.pos_switch(b)+1:self.pos_switch(b+1),:),1);
            end
            
            tBinsList_1 = self.tBinsList-1;
            res_vec = [reshape(fe_resid, self.tBinsNK, 1) ; zeros(self.K*sum(tBinsList_1) + self.K*self.tBinsNtmp,1)];
            
        end
        
        
        %------------------------------------------------------------------
        function gamma = reconstructGamma(self, x)
            % to N Gammas under each other in "classic" orientation
            fe_gamma = reshape(x(1:self.tBinsNK), self.tBinsN, self.K);
            
            gamma = zeros(self.NT, self.K);
            
            for b = 1:self.tBinsN
                gamma(self.pos_switch(b)+1:self.pos_switch(b+1),:) = ones(self.pos_switch(b+1) - self.pos_switch(b), 1) * fe_gamma(b,:);
            end
            
        end
        
        
        %------------------------------------------------------------------
        function [Aeq, beq] = getAeqBeq(self)
            % GETAEQBEQ computes the equality contraints
            %%%%%%%%
            % this code for better understanding of the structure of Aeq
            %%%%%%%%
            %Nrows = self.K*(self.tBins-1)*self.N + self.K*self.tBinsNtmp;
            %Aeq = sparse( zeros(self.tBinsN, self.tBinsNK + Nrows) );
            %beq = sparse(ones(self.tBinsN, 1));
            %
            %for i = 1:self.tBinsN
            %    Aeq(i, i:self.tBinsN:self.tBinsNK) = 1;
            %end
            %%%%%%%
            % this one is 100-times faster than the above formulation
            %%%%%%%
            Nrows = self.tBins_1NK + self.K*self.tBinsNtmp;
            
            aa = cell(1,self.K);
            [aa{:}] = deal( sparse(1:self.tBinsN,1:self.tBinsN,1) );
            
            % here use sparse(zero()) because Nrows sometimes = 0
            Aeq = horzcat( aa{:}, sparse( zeros(self.tBinsN, Nrows) ));
            beq = sparse(ones(self.tBinsN, 1));
            
        end
        
        
        %------------------------------------------------------------------
        function Aneq = getAneq(self)
            % GETANEQ computes the inequality constraints
            
            if self.N > 1 && any(self.tBinsList>1) && ~isinf(self.Cn)
                Aneq1 = self.const_etaCT();
                Aneq2 = self.const_etaGamma();
                Aneq3 = self.const_zetaCJ();
                Aneq4 = self.const_zetaGamma();
                Aneq = vertcat(Aneq1, Aneq2, Aneq3, Aneq4);
                
            elseif self.tBinsN == 1
                Aneq3 = self.const_zetaCJ();
                Aneq4 = self.const_zetaGamma();
                Aneq = vertcat(Aneq3, Aneq4);
                
            elseif self.N == 1 || isinf(self.Cn)
                Aneq1 = self.const_etaCT();
                Aneq2 = self.const_etaGamma();
                Aneq = vertcat(Aneq1, Aneq2);
            end
            
            Aneq = vertcat( sparse(1:size(Aneq,2),1:size(Aneq,2),-1), Aneq);
            
        end
        
        
        %------------------------------------------------------------------
        function bneq = getBneq(self)
            % GETBNEQ
            bneqCT = sparse(1:self.NK, 1, self.Ct);
            bneqCN = sparse(1:self.tBinsK, 1, self.Cn);
            bneq_slack_cn = sparse(2*self.K*self.tBinsNtmp, 1);
            bneq_slack_ct = sparse(2*self.tBins_1NK, 1);
            
            if self.N > 1 && any(self.tBinsList>1) && ~isinf(self.Cn)
                bneq = vertcat(bneqCT, bneq_slack_ct, bneqCN, bneq_slack_cn);
            elseif any(self.tBinsList == 1)
                bneq = vertcat(bneqCN, bneq_slack_cn);
            elseif self.N == 1 || isinf(self.Cn)
                bneq = vertcat(bneqCT, bneq_slack_ct);
            end
            
            nrows = size(self.Aneq,2);
            bneq = vertcat( sparse(nrows, 1, 0 ), bneq );
            
        end
        
        
        %------------------------------------------------------------------
        function Aneq = const_etaCT(self)
            %CONST_ETACT provide the matrix for contrait (10)
            % time jumps of the switching process for each location should
            % be limited by CT.
            % In the paper the vector of unkowns is Omega = (gamma, eta, zeta)
            % time, location, cluster
            % eta belongs to time contraint and thus we need just to sum up
            % for every location the part of the matrix corresponding to Omega:
            % AneqOmega <= CT
            % with
            %       Aneq= |        1...1 0...0                 |
            %             |        0...0 1...1                 |
            %             |   D               .            F   |
            %             |                     .              |
            %             |                       .            |
            %             |                         1...1      |
            %
            %
            % size(Aneq) = K*N x (K*T*N + K*(T-1)*N + K*T*(N-1) )
            % and
            % bneq = CT*ones( N*K, 1 )
            
            D = sparse(self.NK, self.tBinsNK);
            F = sparse(self.NK, self.K*self.tBinsNtmp);
            E = cell(1,self.N);
                        
            for j = 1:self.N
                tBins_j = self.tBinsList(j);
                E{j} = sparse(1,1:(tBins_j-1),1);
            end
            
            Etmp =  blkdiag(E{:});
            EE = cell(1,self.K);
            [EE{:}] = deal(Etmp);
            E = blkdiag(EE{:});
            
            Aneq = horzcat(D,E,F);
            
        end
        
        
        %------------------------------------------------------------------
        function Aneq = const_etaGamma(self)
            %CONST_ETAGAMMA provide the matrix for contraits (11-12)
            % -gamma^(j)(t+1) + gamma^j(t) - eta^j(t) <= 0
            % In the paper the vector of unkowns is Omega = (gamma, eta, zeta).
            % with
            %       Aneq =|  -1 1 0... 0              -1 0...0                       |
            %      %      |   0   . .  0               0 . ..0                       |
            %             |   0 ..0 -1 1               0 ..0-1              C        |
            %             |               -1 1 0... 0         -1 0...0               |
            %             |                0   . .  0          0 . ..0               |
            %             |                0 ..0 -1 1          0 ..0-1               |
            %
            %
            %
            % size(Aneq) = K*(T-1)*N x (K*T*N +  K*(T-1)*N + K*T*(N-1) )
            % and
            % bneq = zeros( K*(T-1)*N, 1 )
            % please note the above formulation is valid only when each
            % location has the same T
            
            
            A = cell(1,self.N);
            B = cell(1,self.N);
            
            for j = 1:self.N
                tBins_j = self.tBinsList(j);
                i_tmp = 1:tBins_j-1;
                
                index_row = [i_tmp i_tmp];
                index_col = [i_tmp i_tmp+1];
                s = [ones(1,tBins_j-1) -1*ones(1,tBins_j-1)];
                
                A{j} = sparse(index_row, index_col, s);
                B{j} = sparse(i_tmp, i_tmp, -1);
            end
            
            Atmp =  blkdiag(A{:});
            
            AA = cell(1,self.K);
            [AA{:}] = deal(Atmp);
            A = blkdiag(AA{:});
            
            
            Btmp =  blkdiag(B{:});
            BB = cell(1,self.K);
            [BB{:}] = deal(Btmp);
            B = blkdiag(BB{:});
            
            % in case Cn is empty -> no space regularisation, C in above
            % form is empty
            if isinf(self.Cn)
                C = [];
            else
                C = sparse(size(A,1), self.K*self.tBinsNtmp);
            end
            
            % to suttiesfy: -gamma^(j)(t+1) + gamma^j(t) - eta^j(t) <= 0
            % and           +gamma^(j)(t+1) - gamma^j(t) - eta^j(t) <= 0
            Aneq = vertcat( horzcat(A,B,C), horzcat(-1*A,B,C) );
            
        end
        
        
        %------------------------------------------------------------------
        function Aneq = const_zetaCJ(self)
            %CONST_ZETACT provide the matrix for contrait (6)
            % space jumps of the switching process for each location are
            % be limited by CJ.
            % In the paper the vector of unkowns is Omega = (gamma, eta, zeta)
            % zeta belongs to space contraint thus for every time step we
            % need to sum up over locations (with a gap of T)
            % for every time step the part of the matrix corresponding to Omega:
            % AneqOmega <= CT
            % with
            %       Aneq= |           1 0...0 1 0...0 ...                          |
            %             |           0 1 0...0 1 0...0 ...                        |
            %             |   J    H                       1 0...0 1 0...0 ...     |
            %             |                                0 1 0...0 1 0...0 ...   |
            %             |                                     .                  |
            %             |                                       .                |
            %
            % thereby the first row sum up zeta for the first time step for the first
            % cluster, the second row for the second time step und so on ...
            % the second block do the same for the second cluster and so on.
            %
            % size(Aneq) = T*K x (K*T*N + K*(T-1)*N + K*T*0.5(N^2-N))
            % and
            % bneq = CJ*ones( T*K, 1 )
            
            if all( self.tBinsList == self.tBinsList(1) )
                J = sparse(self.tBinsK, self.tBinsNK);
                H = sparse(self.tBinsK, self.tBins_1NK);
                
                Ntmp = 0.5*(self.N^2 - self.N);
                
                L = cell(1,Ntmp);
                [L{:}] = deal( sparse(1:self.tBins, 1:self.tBins,1) );
                L = horzcat(L{:});
                
                L_tmp = cell(1,self.K);
                [L_tmp{:}] = deal(L);
                L = blkdiag(L_tmp{:});
                
                Aneq = horzcat(J,H,L);
            else
                error('Space regularization for different time steps is not ready yet')
            end
            
        end
        
        
        %------------------------------------------------------------------
        function Aneq = const_zetaGamma(self)
            %CONST_ZETAGAMMA provide the matrix for contraits (11-12)
            %  gamma^(j+1)(t) + gamma^j(t) - zeta^j(t) <= 0
            % In the paper the vector of unkowns is Omega = (gamma, eta, zeta)
            % with e.g,.:
            %       Aneq= | 1 0...0  0-1... 0          -1 0...0 0.. 0         |
            %             | 1 0...0  0 0 -1 0           0.....0 0-1.0         |
            %             | 1 0...0  0 0  0-1           0.....0 0-1.0         |
            %             | 0 1 0.0  0 0 -1 0     I     0-1...0 0...0         |
            %             | 0 1 ..0  0 0  0-1           0...  0 0 0-1         |
            %             |         .                           .             |
            %             |           .                           .           |
            %
            % size(Aeq) = K*T*0.5(N^2-N) x (K*T*N + K*(T-1)*N + K*T*0.5(N^2-N))
            % and
            % beq = zeros( K*T*0.5(N^2-N), 1 )
            % construct tmp matrices       
            
            if all( self.tBinsList == self.tBinsList(1) )                 
                G = zeros(self.tBinsNtmp, self.tBinsN);
                row_start = 1;
                tBins = self.tBinsList(1);
                
                for t = 1:tBins
                    for j=1:self.N_1;
                        Atmp = zeros(self.N-j, tBins*self.N);
                        Atmp(:, t:tBins:end) = [zeros(self.N-j,j-1),...
                            self.distMatrix(j+1:end,j),...
                            -diag(self.distMatrix(j+1:end,j)) ];
                        
                        G(row_start:row_start+self.N_1-j,:) = Atmp;
                        row_start = 1 + row_start + self.N_1-j;
                    end
                end
                
                GG = cell(1,self.K);
                [GG{:}] = deal(sparse(G));
                G = blkdiag(GG{:});
                
                HH = cell(1,self.K);
                [HH{:}] = deal(sparse(1:self.tBinsNtmp, 1:self.tBinsNtmp, -1));
                H = blkdiag(HH{:});
                
                II = cell(1,self.K);
                [II{:}] = deal(sparse(self.tBinsNtmp, (tBins-1)*self.N));
                I = blkdiag(II{:});
                
                %to suttiesfy: -gamma^(j+1)(t) + gamma^j(t) - zeta^j(t) <= 0
                % and          +gamma^(j+1)(t) - gamma^j(t) - zeta^j(t) <= 0
                Aneq = vertcat(horzcat(-1*G,I,H),horzcat(G,I,H));
            else
                error('Space regularization for different time steps is not ready yet')
            end
            
        end
        
        
    end
end



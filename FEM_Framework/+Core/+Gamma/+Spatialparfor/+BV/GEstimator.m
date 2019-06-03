%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : Factory.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : creation matrices for BV
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
classdef GEstimator < Core.Gamma.Spatialparfor.GEstimator
    %GESTIMATORSBV Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Aeq     % cell of matrices
        beq     % cell of matrices
        Aneq    % cell of matrices
        bneq    % cell of matrices
    end
    
    methods
        
        %------------------------------------------------------------------
        function self = GEstimator(cfg, input, K, Ct, Cn) % gamma cgf
            self = self@Core.Gamma.Spatialparfor.GEstimator(cfg, input, K, Ct, Cn);
            
            for j = 1:self.N
                [self.Aeq{j}, self.beq{j}] = self.getAeqBeq(j);
                self.Aneq{j} = self.getAneq(j);
                self.bneq{j} = self.getBneq(j);
            end 
        end
        
    end
    
    
    methods (Access = 'private')
        
        %------------------------------------------------------------------
        function [Aeq, beq] = getAeqBeq(self, j)
            % GETAEQBEQ computes the equality contraints
            Nrows = (self.tBins(j)-1)*self.K;
            aa = cell(1,self.K);
            [aa{:}] = deal( sparse(1:self.tBins(j),1:self.tBins(j),1) );
            
            % here use sparse(zero()) because Nrows sometimes = 0
            Aeq = horzcat( aa{:}, sparse( zeros(self.tBins(j), Nrows) ));
            beq = sparse(ones(self.tBins(j), 1));
            
        end
        
        
        %------------------------------------------------------------------
        function Aneq = getAneq(self,j)
            % GETANEQ computes the inequality constraints
            Aneq1 = self.const_etaCT(j);
            Aneq2 = self.const_etaGamma(j);
            Aneq = vertcat(Aneq1, Aneq2);
            
            Aneq = vertcat( sparse(1:size(Aneq,2),1:size(Aneq,2),-1), Aneq);
            
        end
        
        
        %------------------------------------------------------------------
        function bneq = getBneq(self, j)
            % GETBNEQ
            bneqCT = sparse(1:self.K, 1, self.Ct);
            bneq_slack_ct = sparse(2*(self.tBins(j)-1)*self.K, 1);
            bneq = vertcat(bneqCT, bneq_slack_ct);
            
            nrows = size(self.Aneq{j},2);
            bneq = vertcat( sparse(nrows, 1, 0 ), bneq );
            
        end
        
        
        %------------------------------------------------------------------
        function Aneq = const_etaCT(self, j)
            %CONST_ETACT provide the matrix for contrait (10)
            D = sparse(self.K, self.tBinsK(j));
            E = sparse(1,1:(self.tBins(j)-1),1);
            
            EE = cell(1,self.K);
            [EE{:}] = deal(E);
            E = blkdiag(EE{:});
            
            Aneq = horzcat(D,E);
            
        end
        
        
        %------------------------------------------------------------------
        function Aneq = const_etaGamma(self, j)
            %CONST_ETAGAMMA provide the matrix for contraits (11-12)
            tBinsj_1 = self.tBins(j)-1;
            i_tmp = 1:tBinsj_1;
            index_row = [i_tmp i_tmp];
            index_col = [i_tmp i_tmp+1];
            s = [ones(1,tBinsj_1) -1*ones(1,tBinsj_1)];
            
            A = sparse(index_row, index_col, s);
            B = sparse(i_tmp, i_tmp, -1);
            
            Atmp =  blkdiag(A);
            
            AA = cell(1,self.K);
            [AA{:}] = deal(Atmp);
            A = blkdiag(AA{:});
            
            Btmp =  blkdiag(B);
            BB = cell(1,self.K);
            [BB{:}] = deal(Btmp);
            B = blkdiag(BB{:});
            
            % to suttiesfy: -gamma^(j)(t+1) + gamma^j(t) - eta^j(t) <= 0
            % and           +gamma^(j)(t+1) - gamma^j(t) - eta^j(t) <= 0
            Aneq = vertcat( horzcat(A,B), horzcat(-1*A,B) );
            
        end
        
    end
end



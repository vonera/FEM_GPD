%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : ResultHandler.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : collect results
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
classdef ResultHandler < handle
    %RESULTHANDLER 
    %   Detailed explanation goes here
    
    properties
        rCell
        next
        criteria_list = {'AIC','AICc','BIC'};
    end
    
    methods
        
        %------------------------------------------------------------------
        function self = ResultHandler()
            self.rCell = cell(1,100);
            self.next = 1;
        end
        
        %------------------------------------------------------------------
        function add(self, result)
            % add result (if not empty, or default) to result cell
            
            % TODO check default result (boot or random init..)
            if ~isempty(result)
                self.rCell{self.next} = result.copy();
                self.next = self.next + 1;
            end
        end
        
        %------------------------------------------------------------------
        function addAll(self, resultHandler)
            % copy all results from another resultHandlet to result cell 
            
            n2add = resultHandler.next - 1;
            
            self.rCell(self.next:self.next+n2add-1) = resultHandler.rCell(1:n2add);
            
            self.next = self.next + n2add;
        end
        
        %------------------------------------------------------------------
        function combineAll(self, resultHandler)
            % replace results by better results for the same optimisation  
            % try to find same results in both sets (own and incomming) and
            % replace in own by incomming if it is better. For that result 
            % has to implement two methods (issimilar, isbetter) if 
            % incomming result is not in own set, it will be added. 
            
        end
        
        %------------------------------------------------------------------
        function printAllInfo(self)
            fprintf('\n');
            for i = 1:self.next - 1
                disp(self.rCell{i}.getInfoString());
            end
        end
        
        %------------------------------------------------------------------
        function computeIC(self, iHandl)
            for i = 1:self.next - 1
                input = iHandl.getFor(self.rCell{i}.param.meta);
                self.rCell{i}.compute_IC(input);
            end
        end
        
        %------------------------------------------------------------------
        function result = getBest(self, criteria)
            
            c = '';
            for cval = self.criteria_list
                if strcmp(criteria,cval{:})
                    c = criteria;
                    break;
                end
            end
            
            result = [];
            
            if ~isempty(c)
                if self.next > 1
                    result = self.rCell{1};
                    for i = 2: self.next-1
                        if eval(['self.rCell{i}.ic.',c,'< result.ic.',c])
                            result = self.rCell{i};
                        end
                    end
                end    
                
            else
                warning(['unknown criteria: ',criteria]);
            end
        end
        
        %------------------------------------------------------------------
        function iterable = getIterable(self)
            iterable = self.rCell(1,1:self.next-1);
        end
    end
end


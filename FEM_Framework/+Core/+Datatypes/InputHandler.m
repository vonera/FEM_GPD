%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : InputHandler.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : create input for different configurations
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
classdef InputHandler < handle
    %INPUTHANDLER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = 'public', SetAccess = 'private')
        cfg
        xt              % a cell of size N (N = number of location) 
        ut              % a cell of size N (N = number of location)
        meta_cell
    end
    
    methods
        
        %------------------------------------------------------------------
        function self = InputHandler(cfg, xt, ut)
            if ~iscell(xt)
                xt = {xt};
            end
            if ~iscell(ut)
                ut = {ut};
            end
            
            self.cfg = cfg;
            self.xt = xt;
            self.ut = ut;
            self.meta_cell = {};
            self.fill_meta_list();
        end
   
        %------------------------------------------------------------------
        function num = getNumOfInputs(self)
            num = length(self.meta_cell);
        end
      
        %------------------------------------------------------------------
        function input = get(self, idx)
            if idx < 0 || idx > length(self.meta_cell)
               error('index exeeded');
            end
            for j = 1:size(self.xt,2)
                % TODO: agly hack but for now its ok
                input(j) = Core.Construction.Factory.createInput(self.cfg, self.xt{j}, self.ut{j}, self.meta_cell{idx});
            end
            
        end
        
        %------------------------------------------------------------------
        function input = getFor(self, meta)
            input = [];
            for idx = 1:self.getNumOfInputs()
                if isequal(self.meta_cell{idx}, meta)
                    input = self.get(idx);
                    break
                end
            end
            
            if isempty(input)
                error('unknown meta');
            end
        end
        
    end
    
    methods (Access = 'protected')
        
        %------------------------------------------------------------------
        function add_meta(self, meta)
            % 
            if ~isa(meta, 'Core.Datatypes.Meta')
                error('Bad type of meta, expected Core.Datatypes.Meta or subclass');
            end
            % TODO: make it nice 
            self.meta_cell{end+1} = meta;
            
        end
        
        %------------------------------------------------------------------
        % to customize overwrite this method!!
        function fill_meta_list(self)
            self.add_meta(Core.Datatypes.Meta(size(self.xt{1},2)));
        end
    end
end


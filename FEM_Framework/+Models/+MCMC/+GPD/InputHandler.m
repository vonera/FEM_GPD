%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : InputHandler.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : create meta
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
classdef InputHandler < Core.Datatypes.InputHandler
    %INPUTHANDLER MCMC.GPD inputHandler
    %   Detailed explanation goes here
    
    methods
        
        %------------------------------------------------------------------
        function self = InputHandler(cfg, xt, ut)
            self = self@Core.Datatypes.InputHandler(cfg, xt, ut);
        end
        
    end
    
    methods (Access = 'protected')
        
        %------------------------------------------------------------------
        % overwriten method!!
        function fill_meta_list(self)
            for index = 1:size(self.cfg.model.ut_combis{1},2)
                umask_ksi = self.cfg.model.ut_combis{1}{index};
                umask_si = self.cfg.model.ut_combis{2}{index};
                
                N = size(self.xt,2);
                Tlist = zeros(1,N);
                for j = 1:N
                    Tlist(j) = size(self.xt{j},2);
                end

                self.add_meta(Models.MCMC.GPD.Meta(Tlist, umask_ksi, umask_si, self.cfg.model.lassoFLAG, self.cfg.model.ridgeFLAG, self.cfg.model.paramTHRESHOLD));
            end
        end
        
    end
end

%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : MFitConfig.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : configuration for alternating optimization
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
classdef MFitConfig
    %MFitConfig Config for Model Fit
    %   Detailed explanation goes here
    
    properties
        annealing       % number
        annealingForK1  % boolean
        subspaceItr     % number 
        subspaceEps     % number
        Klist           % list of numbers
        Ctlist          % list of numbers
        Cnlist          % list of numbers
    end
    
    methods
        function self = MFitConfig(annealing, annealingForK1, subspaceItr,...
                subspaceEps, Klist, Ctlist, Cnlist)
            if annealing == 1 && size(Ctlist,2) > 1 
                warning('DANGER!! for only 1 annealing and len(Ctlist)>1 you get less randomniss!!')
            end
            self.annealing = annealing;
            self.annealingForK1 = annealingForK1;
            self.subspaceItr = subspaceItr;
            self.subspaceEps = subspaceEps;
            self.Klist = Klist;
            self.Ctlist = Ctlist;
            self.Cnlist = Cnlist;
        end
    end
end


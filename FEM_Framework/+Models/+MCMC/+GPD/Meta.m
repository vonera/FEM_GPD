%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : Meta.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : mask for parameters
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
classdef Meta < Core.Datatypes.Meta
    %META Meta for MCMC GPD full non stat 
    %   Detailed explanation goes here
    
    properties (GetAccess = 'public', SetAccess = 'private')
        umask_ksi         % list of numbers index list [1,3,7]
        umask_si          % list of numbers index list [1,3,7]
        Tlist             % list of T for each location
        lassoFLAG         % constant, >= 0
        ridgeFLAG         % constant, >= 0
        paramTHRESHOLD    % is uded for LASSO/Ridge regression, all params below
                          % the threshold are set to 0 and do not count for
                          % num of params
    end

    methods
        function self = Meta(Tlist, umask_ksi, umask_si, lassoFLAG, ridgeFLAG, paramTHRESHOLD)
            
            self = self@Core.Datatypes.Meta(0);
            self.Tlist = Tlist;    
            self.umask_ksi = umask_ksi;
            self.umask_si = umask_si;
            self.lassoFLAG = lassoFLAG; 
            self.ridgeFLAG = ridgeFLAG;
            self.paramTHRESHOLD = paramTHRESHOLD;
        end
        
        
        function infostring = getInfoString(self)
            infostring = [getInfoString@Core.Datatypes.Meta(self), ' umask_ksi = ',mat2str(self.umask_ksi),' umask_si = ', mat2str(self.umask_si)];
        end
    end
end



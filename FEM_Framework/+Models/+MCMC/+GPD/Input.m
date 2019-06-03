%------------------------------------------------------------------------
% Project     : FEM-BV Framework based on the Theoretical development of
%                                                       I.Horenko and group
% Version     : 1.0
%
% Category    : time series analysis
% Filename    : Input.m
% Created by  : O.Kaiser <gulgosha@gmail.com>, D.Igdalov <vonera@gmail.com>
% Description : set covas for params
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
classdef Input < Core.Datatypes.Input
    %INPUT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ut_ksi
        ut_si
        
        d_ksi 
        d_si
    end
    
    methods
        %------------------------------------------------------------------
        function self = Input(xt, ut, meta)
            self = self@Core.Datatypes.Input(xt, meta);
            
            self.ut_ksi = [ ones(1,size(ut,2)) ; ut(meta.umask_ksi,:)];
            self.ut_si  = [ ones(1,size(ut,2)) ; ut(meta.umask_si, :)];
            
            self.d_ksi =  size(self.ut_ksi,1);
            self.d_si = size(self.ut_si,1);
        end
        
        %------------------------------------------------------------------
        function [xti, uti_ksi, uti_si] = getForTimeSteps(self, logic_vec_times)
            xti     = self.xt(:, logic_vec_times);
            uti_ksi = self.ut_ksi(:, logic_vec_times);
            uti_si  = self.ut_si(:, logic_vec_times);
        end
        
        %------------------------------------------------------------------
        function [d_ksi, d_si] = getUdims(self)
            d_ksi = self.d_ksi;
            d_si = self.d_si;
        end
        
    end
end

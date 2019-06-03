classdef ConvexSampler
    %CONVEXSAMPLER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function x = convexSampler(A,b,x)
            % Copyright (c) 2011, Michael Weitzel
            % All rights reserved.
            %
            % Redistribution and use in source and binary forms, with or without
            % modification, are permitted provided that the following conditions are
            % met:
            %
            %     * Redistributions of source code must retain the above copyright
            %       notice, this list of conditions and the following disclaimer.
            %     * Redistributions in binary form must reproduce the above copyright
            %       notice, this list of conditions and the following disclaimer in
            %       the documentation and/or other materials provided with the distribution
            %
            % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
            % AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
            % IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
            % ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
            % LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
            % CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
            % SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
            % INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
            % CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
            % ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
            % POSSIBILITY OF SUCH DAMAGE.
            % including sampleStepHitAndRun(A,b,x), sampleStepGibbs(A,b,x),
            % distToConstraints(A,b,x,r)
            
            % Sampling of uniformly distributed vectors x satisfying the linear system
            % of inequalities A.x <= b. It is assumed that the convex region described
            % by A.x <= b is not unbounded.
            if ~all(A*x<=b)
                error('initial point infeasible');                
            end
            x = Models.MCMC.Tools.ConvexSampler.sampleStepGibbs(A,b,x);
            %x = sampleStepHitAndRun(A,b,x);
        end
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function x = sampleStepHitAndRun(A,b,x)
            % Generate a random direction vector (uniformly covering the surface of a
            % unit sphere)
            r = randn(size(A,2),1);
            r = r/norm(r);
            [d_neg,d_pos] = Models.MCMC.Tools.ConvexSampler.distToConstraints(A,b,x,r);
            % a random number [-d_neg,d_pos)
            t = -d_neg+(d_pos+d_neg)*rand(1);
            % next sample
            x = x + t*r;
        end
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function x = sampleStepGibbs(A,b,x)
            r = zeros(size(A,2),1);
            for j=1:size(A,2)
                % generate a unit vector e_j
                r(j) = 1;
                if j>1, r(j-1) = 0; end
                [d_neg,d_pos] = Models.MCMC.Tools.ConvexSampler.distToConstraints(A,b,x,r);
                % a random number [-d_neg,d_pos)
                t = -d_neg+(d_pos+d_neg)*rand(1);
                % element j of next sample
                x(j) = x(j) + t;
            end
        end
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        function [d_neg,d_pos] = distToConstraints(A,b,x,r)
            % Determine the minimum distances between point x and the (two) intersection
            % points of constraints A.x <= b and line g(t)=x+t*r (in directions -r and +r).
            d_pos = Inf;
            d_neg = Inf;
            for i=1:size(A,1)
                rho = (A(i,:)*x-b(i))/(A(i,:)*r);
                xci = x - r*rho; % intersection point with constraint i
                d_x_xci = norm(xci - x); % distance
                if r'*(xci - x) >= 0.0 % xci is in direction +r
                    if d_x_xci < d_pos
                        d_pos = d_x_xci;
                    end
                else % xci is in direction -r
                    if d_x_xci < d_neg
                        d_neg = d_x_xci;
                    end
                end
            end
        end
    end
end


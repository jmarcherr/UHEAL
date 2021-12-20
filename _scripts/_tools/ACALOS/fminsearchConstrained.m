% Constrained function minimization
%
% Based on fminsearchbnd.m (BSD license)
% Copyright (c) 2006, John D'Errico
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
%
% v0.1 Stephan Ewert, 2014


function [x,fval,exitflag,output]=fminsearchConstrained(fun,x0,lb,ub,options,varargin)

bounded = zeros(size(ub)); % 0 = unconstrained, 1 = constrained
for i=1:length(bounded)
    bounded(i) = any([isfinite(lb(i)) isfinite(ub(i))]);
end

%bounded

% transform starting values using the inverse of sinetransform
x0t = x0;

for i = 1:length(x0)
    if bounded(i)
        % lower and upper bounds
        if x0(i) < lb(i)
            x0t(i) = -pi/2;
        elseif x0(i)> ub(i)
            x0t(i) = pi/2;
        else
            x0t(i) = 2*(x0(i) - lb(i))/(ub(i)-lb(i)) - 1;   % maps to -1 ... 1
            %x0t(i) = asin(max(-1,min(1,x0t(i))));      % maps to -pi/2 ... pi/2, this can be zero which is not good for fminsearch
            x0t(i) = 2*pi + asin(max(-1,min(1,x0t(i))));    % shift by one period to avoid zero start values
        end
        % we could put the additive term here would be more logic
        % x0t(i) = x0t(i) + 2*pi;
    end
end


% our internal constraint function defined inline
constraintfun = @(x, varargin) fun(sinetransform(x,lb,ub,bounded), varargin{:});

% now call fminsearch
% DO: added options as an parameter
[x,fval,exitflag,output] = fminsearch(constraintfun,x0t,options,varargin);

% transform the continuous variable in into the original bounds
x = sinetransform(x, lb, ub, bounded);

% subfunction
function out = sinetransform(x,lb,ub,bounded)

% transform into a constraint range with lower and upper bounds
% using sine transform

% initialize to original values
out = x;

for i = 1:length(x)
    if bounded(i) == 1
        out(i) = (sin(x(i))+1)/2; % periodically ranges from 0 to 1
        out(i) = out(i)*(ub(i) - lb(i)) + lb(i); % ranges from lower bound to upper bound
    end
end

% eof
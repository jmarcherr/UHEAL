% fit = loudness_function_bh2002(x, fitparams, inverse)
%
% function calculates the loudness function according to fitparams
%

%
% PARAMETERS:
%     x:
%       either levels to calculate CU
%       if inverse = true, x contains CU and calculates levels
%     fitparams:
%       fitparams =[lcut, mlow, mhigh]
%       mlow and mhigh have to be positive values
%       if negative values are given, the absolute value will be used.
%       the smallest value for mlow and mhigh is set to 0.001 CU/dB
%     inverse:
%           activates the inverse loudness function
%
% OUTPUT:
%       y:
%           either CU (inverse=false) or levels (inverse=true)
% Authors: Dirk Oetting, SE 16.05.2013 16:27, DO/SE 24.07.2013

% includes DO fixes

% ------------------------------------------------------------------------------
%  Adaptive Categorical Loudness Scaling Procedure for AFC for Mathwork's MATLAB
% 
%  Author(s): Stephan Ewert, Dirk Oetting
% 
%  Copyright (c) 2013-2014, Stephan Ewert, Dirk Oetting
%  All rights reserved.
% 
%  This work is licensed under the 
%  Creative Commons Attribution-NonCommercial-NoDerivs 4.0 International License (CC BY-NC-ND 4.0). 
%  To view a copy of this license, visit
%  http://creativecommons.org/licenses/by-nc-nd/4.0/ or send a letter to Creative
%  Commons, 444 Castro Street, Suite 900, Mountain View, California, 94041, USA.
% ------------------------------------------------------------------------------

function [y, failed] = loudness_function_bh2002(x, fitparams, inverse)

if nargin < 3
    inverse = false;
end

Lcut = fitparams(1);
m_lo = fitparams(2);
m_hi = fitparams(3);

% DO check values of given parameters

if m_lo <= 0 
    if m_lo == 0
        m_lo = 0.001;
        %warning('m_low has to be positive, but is zero! m_low was set to 0.001!');
    elseif m_lo < 0
        m_lo = - m_lo;
        %warning('m_low has to be positive! Absolute value of m_low was used!');
    end
    
end
if m_hi <= 0 
    if m_hi == 0
        m_hi = 0.001;
        %warning('m_high has to be positive, but is zero! m_high was set to 0.001!');
    elseif m_hi < 0
        m_hi = - m_hi;
        %warning('m_high has to be positive! Absolute value of m_high was used!');
    end
end

% calculate point where Bezier function should go through
L15 = y2x_lin(15, Lcut, 25, m_lo);
L35 = y2x_lin(35, Lcut, 25, m_hi);
C = [L15 Lcut L35; 15 25 35];


% SE 16.05.2013 16:28 failed flag
failed = 0;

if ~inverse
    cu = ones(size(x))*NaN;
    if m_lo == m_hi
        cu = x2y_lin(x, Lcut, 25, m_lo);
    else
        % calculate all values below Lcut
        idx = find(x <= Lcut);
        cu(idx) = x2y_lin(x(idx), Lcut, 25, m_lo);

        % find all values above Lcut
        idx = find(x > Lcut);
        cu(idx) = x2y_lin(x(idx), Lcut, 25, m_hi);

        % calculate transition range between 15 and 35 CU
        idx = find(cu>15 & cu<35);
        if any(idx)
            [ cu(idx), failed ] = BezierX2YFor3ControlPoints(x(idx), C);
        end
    end

    if failed
        y =[];
        return;
    end

    % limit to 0 to 50
    cu = (cu<0)*0 + (cu>=0 & cu<=50).*cu + (cu>50)*50;
    y = cu;
else
    % limit x from 0 to 50:
    x = (x<0)*0 + (x>=0 & x<=50).*x + (x>50)*50;

    if m_lo == m_hi
        y = y2x_lin(x, Lcut, 25, m_hi);
    else
        levels = ones(size(x))*NaN;
        % find all values below CU = 15
        idx = find(x <= 15);
        levels(idx) = y2x_lin(x(idx), Lcut, 25, m_lo);

        % find all values above CU = 35
        idx = find(x >= 35);
        levels(idx) = y2x_lin(x(idx), Lcut, 25, m_hi);
        % calculate transition range between 15 and 35 CU
        idx = find(x>15 & x<35);

        if any(idx)
            [ levels(idx), failed ] = BezierX2YFor3ControlPoints(x(idx), C, true);
        end

        if failed
            y =[];
            return;
        end

        y = levels;
    end
end

return

function y = x2y_lin(x, x0, y0, m)
y = y0 + m*(x-x0);
return

function x = y2x_lin(y, x0, y0, m)
x = (y-y0)/m + x0;
return

function [y, failed] = BezierX2YFor3ControlPoints(x, C, inverse)
if nargin < 3
    inverse = false;
end

failed = 0;

% calculate t of x
t = NaN;
y = zeros(size(x));
y0 = C(2,1);
y1 = 2*C(2,2)-2*C(2,1);
y2 = C(2,1)-2*C(2,2)+C(2,3);

x2 = C(1,1)-2*C(1,2)+C(1,3);
x1 = 2*C(1,2)-2*C(1,1);
x0 = C(1,1);

if ~inverse
    % m_low and m_high are not identical
    if x2 ~= 0
        t1 = -x1/(2*x2) + 0.5*sqrt((x1/x2)^2 - 4*(x0-x)/x2);
        t2 = -x1/(2*x2) - 0.5*sqrt((x1/x2)^2 - 4*(x0-x)/x2);
        % check if t values are between zero and one
        if any(t1) && all(imag(t1)==0) && min(t1)>=0 && max(t1)<=1
            t = t1;
        elseif any(t2) && all(imag(t2)==0) && min(t2)>=0 && max(t2)<=1
            t = t2;
        else
            % SE 16.05.2013 16:31 changed error handling
            warning('something strange happens in BezierX2YFor3ControlPoints');
            failed = 1;
            y = [];
            return;
        end
        % L_cut ~= L15 should always be the case
    elseif x1 ~= 0
        t = (x-x0)/x1;
    end
    % if inverse
else
    t = x/y1 - y0/y1;
end
% calculate y of t
idx = find(t >= 0 & t <= 1);
if ~inverse
    y(idx) = y2*t(idx).^2 + y1*t(idx) + y0;
else
    y(idx) = x2*((t+x1/(2*x2)).^2) - (x1^2)/(4*x2) + x0;
end


return


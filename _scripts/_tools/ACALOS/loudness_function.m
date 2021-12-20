% fit = loudness_function(x, fitparams, inverse) 
%
% function calculates the loudness function according to fitparams
%

%     
% PARAMETERS:
%     x:
%       either levels to calculate CU
%       if inverse = true, x contains CU and calculates levels
%     fitparams:
% 		  [m_low, HTL, m_high, (UCL)]
%		  UCL is optional and will override m_high
%         Values describing the loudness function.
%     inverse:
%           activates the inverse loudness function
%     
% OUTPUT:
%       y:
%           either CU (inverse=false) or levels (inverse=true)
% Authors: Dirk Oetting

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

function y = loudness_function(x, fitparams, inverse) 

if nargin < 3
    inverse = false;
end
CP = 25;
m_lo = fitparams(1);
HTL = fitparams(2);
b = 2.5 - m_lo*HTL;
Lcut = (CP - b)/m_lo;

if length(fitparams) == 4
    UCL = fitparams(4);
    if Lcut >= UCL 
        disp('warning: wrong setting for UCL, function converted to linear');
        UCL = (50 - b)/m_lo;
    end
    m_hi = (50-CP)/(UCL - Lcut);
else
    m_hi = fitparams(3);
end

fitparams = [Lcut m_lo m_hi];
y = loudness_function_bh2002(x, fitparams, inverse);


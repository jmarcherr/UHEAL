% fit = fit_loudness_function(measured_data, fit_mode)
% version: 0.92
% This work is based on the publication
% "Optimized loudness-function estimation for categorical loudness scaling
% data." (2014) by Dirk Oetting, Thomas Brand and Stephan D. Ewert
% DOI: 10.1016/j.heares.2014.07.003

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

function [fit,rawData] = fit_loudness_function(measured_data, fit_mode)

% this function creates a new fit structure from loudness scaling data points.
% A fit structure can be used to transform sound levels to categorical
% units and vice versa.
%
%     
% PARAMETERS:
%     measured_data:
%         The data points from the categorial loudness scaling. 
%         [Level1 Response1 Level2 Response2 ...]
%     fit_mode:
%          BY: Brand and Hohmann (2001) fitting function
%          BX: same as BY but minimize error in x-direction (level
%          direction)
%          BTX: BX + threshold estimation, minimize error in x-direction
%          BTUX: BTX + UCL estmation, minimize error in
%          x-direction (recommended for hearing-aid fitting)
%          BTPX: threshold estimation, UCL estmation using data from Pascoe (1998)
%          minimize error in x-direction
%          BTUY: threshold estimation, UCL estmation, minimize error in
%          y-direction
%     
% OUTPUT:
%     fit:
%         Fit structure describing the loudness function.
%         fit = [Lcut m_low m_high];

%  Authors: 
%         Dirk Oetting, 2013, dirk.oetting@idmt.fraunhofer.de
%         Stephan D. Ewert  
%         Thomas Bisitz
%
%

if nargin < 2
    fit_mode = 'BTUX';
end
fit_mode = upper(fit_mode);

if  isequal('BY', fit_mode)
  [Lcut, m_low, m_high] = bh_fit_2001(measured_data);
  HTL = Lcut-22.5*1/m_low;
  fit = [m_low HTL m_high];
elseif  isequal('BX', fit_mode)
  [m_low,HTL,m_high] = bezier_fit(measured_data,'x_v11');
  fit = [m_low HTL m_high];
elseif  isequal('BTX', fit_mode)
  [m_low,HTL,m_high] = bezier_fit(measured_data,'x_v3');
  fit = [m_low HTL m_high];
elseif  isequal('BTUX', fit_mode)
  [m_low,HTL,m_high] = bezier_fit(measured_data,'x_v36');
  fit = [m_low HTL m_high];
elseif  isequal('BTUY', fit_mode)
  [m_low,HTL,m_high] = bezier_fit(measured_data,'y_v36');
  fit = [m_low HTL m_high];
elseif  isequal('BTPX', fit_mode)
  [m_low,HTL,m_high] = bezier_fit(measured_data,'x_v37');
  fit = [m_low HTL m_high];
else
    disp(['The fit mode "' fit_mode '" is not available.']);
end
% convert values according to function given in Brand and Hohmann (2001)
m_lo = fit(1);
HTL = fit(2);
m_hi = fit(3);
Lcut = HTL + 22.5/m_lo;      
fit = [Lcut m_lo m_hi];
rawData = [measured_data(1:2:end)',measured_data(2:2:end)'];


function [Lcut,Mlow,Mhigh] = bh_fit_2001(daten)
% function [Lcut,Mlow,Mhigh] = bh_fit_2001(daten)
%
% 2: modifizierte Kostenfunktion
% Berechnet [Lcut,Mlow,Mhigh] nach "Original"-Verfahren
%  (Summe der Fehlerquadrate minimieren)
% Messwerte:  daten = [level1, cu1, level2, cu2, ... ]
%
% Thomas Bisitz, Mai/Juni 2005
%

% Startwerte für Optimierung:
%fitparams =[Lcut,Mlow,Mhigh];
fitparams_start =[75,0.35,0.65];

[fitparams, exval, exitflag, output] = ...
    fminsearch( @(fitparams)costfc2(fitparams,daten),fitparams_start);

Lcut = fitparams(1);
Mlow = fitparams(2);
Mhigh = fitparams(3);



function [m_low,HTL,m_high] = bezier_fit(daten,mode)
%
% calculate loudness function (straigh lines with Bezier transition)
% Dirk Oetting / Juni 2012
%
% Startwerte für Optimierung:
% fitparams =[HTL,m_low,UCL];


pascoe_fit = false;
if isequal('x_v37', mode)
    pascoe_fit = true;
    mode = 'x_v3';
end

cu = daten(2:2:end);
cu = cu(:);
cu_all = cu;
levels_dB = daten(1:2:end);
levels_dB=levels_dB(:);

% init randomized timer
if exist('rng','builtin')==5
    rng(sum(100*clock))
else
    rand('twister',sum(100*clock))
end
idx = find(cu<=25);

% check if cu and level values are not equal
if any(diff(cu(idx))>0) && any(levels_dB(idx)>0) && length(idx) > 5
    a1 = polyfit(cu(idx),levels_dB(idx),1);
    m1 = 1/a1(1);
    b1 = -1/a1(1)*a1(2);

    a2 = polyfit(levels_dB(idx),cu(idx),1);
    m2 = a2(1);
    b2 = a2(2);

    m_low = mean([m1 m2]);
    htl = mean([-b1/m1 -b2/m2]);
else
    m_low = 1;
    htl = mean(levels_dB(idx));
end

mins            = [0.2  -20  0.2];
maxs            = [5     100   5];


if isequal('x_v11', mode)
    % no weighting
    optFn = @(fitparams,daten) costfcn_bezier_x(fitparams,daten,false,true);
elseif isequal('y_v11', mode)
    % no weighting
     optFn = @(fitparams,daten) costfcn_bezier_y(fitparams,daten);
elseif isequal('_v3', mode(2:4))
    % find a good estimate for the hearing threshold
    optFn = @(fitparams,daten) costfc_psychometric(fitparams,daten);
    searchOptions = struct('Display','none');
    fitparams_start = [0.5 0];
    mins_htl = [0.4 -30];
    maxs_htl = [1 100];

    % use data points up to responses of soft
    htl_levels = levels_dB;
    htl_probability = cu~=0;
    range_of_htl = 0;
    if sum(htl_probability==0) == 0
        % if no responses "not heard" are given, select 5 dB below the
        % lowest response as a "not heard" response
        fitparams = [0.4 min(htl_levels)-5];
        estimated_htl = fitparams(2);
        mins(2) = estimated_htl-range_of_htl;
        maxs(2) = estimated_htl+range_of_htl;
    else
        daten_htl = reshape([htl_levels htl_probability]',length(htl_levels)*2,1);
        [fitparams,fval,exitflag,output] = fminsearchConstrained(optFn,fitparams_start,mins_htl,maxs_htl,searchOptions,daten_htl);
        estimated_htl = fitparams(2);
        if estimated_htl > min(htl_levels) && estimated_htl < max(htl_levels)
 
            % mins max of hearing threshold estimation
            mins(2) = estimated_htl-range_of_htl;
            maxs(2) = estimated_htl+range_of_htl;
        else
            % don't limit the min and max values
            disp('Fitting method BTUX: HTL estimation failed for this dataset. Threshold estimation was not applied.');
        end
    end

    % remove data points with CU = 0
    idx_remove = find(cu==0);
    cu(idx_remove) = [];
    levels_dB(idx_remove) = [];
    daten = reshape([levels_dB cu]',2*length(cu),1);

    
    if isequal('_v36', mode(2:end))
        % set upper slope if less than 4 data points are available above 35
        % CU
        if  sum(cu_all>=35) < 4 
            % assume median slope
            mins(3) = 1.53;
            maxs(3) = mins(3);
            disp(['Fitting method BTUX: not enough data points in the upper loudness range. m_high was fixed to ' num2str(mins(3)) ' CU/dB']);
        end
    end

    if strcmp(mode(1),'x')
        optFn = @(fitparams,daten) costfcn_bezier_x(fitparams,daten,true,true);
    elseif strcmp(mode(1),'y')
         optFn = @(fitparams,daten) costfcn_bezier_y(fitparams,daten);
    end    
end



number_of_runs = 10;
error_x = nan(1,number_of_runs);
fitparams = [];

for ik = 1:number_of_runs
    % slope of loudness function
    fitparams_start(1) = m_low + (rand(1)-0.5)*0.05;
    % level of htl
    if isequal('_v11', mode(2:end))
        fitparams_start(2) = htl + (rand(1)-0.5)*5;
        fitparams_start(3) = mins(3) + (maxs(3) - mins(3))*rand(1);
    elseif isequal('_v3', mode(2:end)) || isequal('_v36', mode(2:end))
        fitparams_start(2) = mins(2) + (maxs(2) - mins(2))*rand(1);
        fitparams_start(3) = mins(3) + (maxs(3) - mins(3))*rand(1);
    end
    searchOptions = struct('Display','none');
    [fitparams(ik,:),error_x(ik),exitflag,output] = fminsearchConstrained(optFn,fitparams_start,mins,maxs,searchOptions,daten);
end
[mx,idx] = min(error_x);
fit_params = fitparams(idx,:);

if pascoe_fit % v_37
    re_fit = false;
    MaxUCL = 140;
     if loudness_function(50,fit_params,true)>MaxUCL
         UCL = MaxUCL;
         re_fit = true;
     end  
    if fit_params(3) < 0.25 || sum(cu_all>=35) < 4
      UCL = PascoeUCL(fit_params(2));
      re_fit = true;
    end
    if re_fit
        mins(3) = UCL;
        maxs(3) = mins(3);
        max_Lcut = UCL - 5/25;
        mins(1) = 22.5/(max_Lcut - mins(2));

        [fit_params,fval,exitflag,output] = fminsearchConstrained(optFn,fitparams_start,mins,maxs,searchOptions,daten);
        disp('Fitting method BTUX: re-fit using UCL estimatin of Pascoe (1988)');
        % convert from UCL to m_high
        b = 2.5 - fit_params(1)*fit_params(2);
        Lcut = (25 - b)/fit_params(1);   
        fit_params(3) = 25/(fit_params(3) - Lcut);
        % limit slope to 5 CU/dB
        fit_params(3) = min(fit_params(3),5);
    end
end
 
m_low = fit_params(1);
HTL = fit_params(2);
m_high = fit_params(3);







function [x] = costfc2(fitparams,daten)
% function [x]=sumerr2(fitparams,daten)
% (zu minimierende Kostenfunktion)
% gibt Summe der quadratischen Fehler des Fit-Modells zu den Daten aus
% verwendet linear fortgeführte klsku.dll
% fitparams =[Lcut,Mlow,Mhigh]
% daten = [level1, cu1, level2, cu2, ... ]

level = daten(1:2:end);
cu = daten(2:2:end);

cu_fit = loudness_function_bh2002(level,fitparams);

%% lineare Fortführung von klsku:
x2= loudness_function_bh2002(50,fitparams,true); %get level at CU = 50
x0=loudness_function_bh2002(0,fitparams,true); %get level at CU = 0

% linear extrapolation of cu values
cu_fit(level<x0) = fitparams(2)*(level(level<x0) - x0) + 0;
cu_fit(level>x2) = fitparams(3)*(level(level>x2) - x2) + 50;

% set points, where the measured data is 50 and the fit function is bigger
% than 50 to the constant value 50, so that the difference cu-cu_fit = 0
cu_fit(cu==50 & (cu_fit > cu))   = 50;

% set points, where the measured data is 0 and the fit function is lower
% than 0 to the constant value 0, so that the difference cu-cu_fit = 0
cu_fit(cu==0 & (cu_fit < cu))    = 0;

% calculate the distance between cu and cu_fit
delta_x = cu_fit - cu;
x = sum((delta_x).^2);


function p = psychometric_function(fitparams,x)
slope = fitparams(1);
threshold = fitparams(2);
p = 1./(1+exp(-slope.*(x-threshold)));

function [x] = costfc_psychometric(fitparams,x)
% calculates the error of the psychometric function
daten = x{1};
level = daten(1:2:end);
probability = daten(2:2:end);


psychometric_fit = psychometric_function(fitparams, level);
x = sum((psychometric_fit - probability).^2);



function [x] = costfcn_bezier_y(fitparams,x)
% function [x]=sumerr2(fitparams,daten)
% (zu minimierende Kostenfunktion)
% gibt Summe der quadratischen Fehler des Fit-Modells zu den Daten aus
% verwendet linear fortgeführte klsku.dll
% fitparams =[Lcut,Mlow,Mhigh]
% daten = [level1, cu1, level2, cu2, ... ]
daten = x{1};

level = daten(1:2:end);
cu = daten(2:2:end);

cu_fit =  loudness_function(level,fitparams);

% lineare Fortführung von klsku:
x2 = loudness_function(50,fitparams,true); % get level at CU = 50
x0 = loudness_function(0,fitparams,true);  % get level at CU = 0

% linear extrapolation of cu values
cu_fit(level<x0) = fitparams(1)*(level(level<x0) - x0) + 0;
cu_fit(level>x2) = fitparams(3)*(level(level>x2) - x2) + 50;

% set points, where the measured data is 50 and the fit function is bigger
% than 50 to the constant value 50, so that the difference cu-cu_fit = 0
cu_fit(cu==50 & (cu_fit > cu))   = 50;

% set points, where the measured data is 0 and the fit function is lower
% than 0 to the constant value 0, so that the difference cu-cu_fit = 0
cu_fit(cu==0 & (cu_fit < cu))    = 0;

% calculate the distance between cu and cu_fit
delta_x = cu_fit - cu;
x = sum((delta_x).^2);


function [x] = costfcn_bezier_x(fitparams,x,weighting,limit_50,outliner_range)
% calculates the error of the broken stick function
if nargin < 5
    outliner_range = 40;
end
if nargin < 4
    limit_50 = true;
end
if nargin < 3
    weighting = false;
end

daten = x{1};

level = daten(1:2:end);
cu = daten(2:2:end);

level_fit = loudness_function(cu,fitparams,true);
delta_x = level_fit - level;

cus = loudness_function([0 50],fitparams,true);
cu0_level = cus(1);
UCL = cus(2);
% set contribution from error function to zeros for levels which are
% outside the fit-function range.

if limit_50==true
    idx_cu = find(cu==50);
    idx_level = find(level>UCL);
    idx = intersect(idx_cu,idx_level);
    delta_x(idx) = 0;
end

% limit contribution of data point below loudness 0 CU as 0
idx_cu = find(cu==0);
idx_level = find(level<cu0_level);
idx = intersect(idx_cu,idx_level);
delta_x(idx) = 0;

% limit contribution to error of 40 dB
delta_x(find(abs(delta_x)>outliner_range)) = outliner_range;

if weighting
        cus = 0:5:50;
        
        std_kls = zeros(11,1);
        for m = 1:11
            idx = (cu==cus(m));
            if ~any(idx)
                continue;
            end
            variante = 2;
            if variante == 1
                std_kls(m) = std(level(idx));
                if std_kls(m) < 3
                    std_kls(m) = 3;
                end
            else
                std_mean = [NaN 8.5	11.0 11.8 9.1 6.9 5.2 4.3 3.9 3.4 3.6];
                std_kls(m) = std_mean(m);
            end
            delta_x(idx) = delta_x(idx)./(std_kls(m));
        end
end

x = sum(delta_x.^2);

function UCL=PascoeUCL(HTL,mode)
%
%	UCL=PascoeUCL(HTL,mode)
%
%	schätzt die Unbehaglichkeitsschwelle UCL aus der Hörschwelle HTL ab
%	nach:
%	Pascoe, D. P., (1988), Clinical measurements of the auditory dynamic range
%       and their relation to formulas for hearing aid gain, In: J.H. Jensen (Editor),
%       Hearing Aid fitting, 13th Danavox Symposium, 129-152.
%
%	mode = 'smoothed': einfache Abschätzung nach Tabelle 4 (pooled data)

if nargin<2
	mode = 'smoothed';
end

if strcmp(mode,'smoothed')
	UCLverHTL=[...
		-100	100;
		40	100;
		120	140];
end
UCL=interp1(UCLverHTL(:,1),UCLverHTL(:,2),HTL);

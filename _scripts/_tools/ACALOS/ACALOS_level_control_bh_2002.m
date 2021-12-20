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

function [stACALOS] = ACALOS_level_control_bh_2002(stACALOS)

% attach patient answer to stACALOS.CU_response in calling function
% initiate function by using init_ACALOS_level_control
% the current presentation level is stored in stACALOS.new_level

% incorporates fixes by Dirk Oetting (DO)
% Authors: René Asendorf, Dirk Oetting, Stephan D. Ewert
% Version: 1.0rc

% fail safe: intialize stACALOS if not initialized
% DOES NOT HAPPEN IN AFC, stACALOS is initialized in acls_init
if ~isfield(stACALOS,'phase')
    stACALOS.new_level = 65;
    stACALOS.level_history = [stACALOS.new_level];
    stACALOS.CU_response =[];
    stACALOS.up_flag = [0];
    stACALOS.down_flag = [0];
    stACALOS.phase = 1;
    stACALOS.finished = false;
    stACALOS.iterations = 2;
    % maximal presentation level, only use values that can be reach from
    % the first level in 5 dB steps
    stACALOS.maximal_presentation_level = 105;
    stACALOS.L_max = stACALOS.maximal_presentation_level;
    % minimal presentation level.
    stACALOS.minimal_presentation_level = -10;
    stACALOS.levels_phase3 = [];
    
    % alternate between high and low track (0) or randomly switch between
    % tracks (1), Brand and Hohmann (2002) used (0)
    stACALOS.randomizeRangeEstimation = 0;
    stACALOS.levelRounding = 1;


    % T Wittkop 07.05.2013 09:55, changes in OMA 2013 version
    % At the end of phase 1, the actual CU's of the last responses
    %      are used (instead fix 5 and 50, approved by T.Brand).
    %
    % set variable to true if wanted
    stACALOS.useMeasuredCU_phase3 = false;
    stACALOS.failed = false;

    return
else
    % get some variables from struct
    phase = stACALOS.phase;
    level_history = stACALOS.level_history;
    CU_response = stACALOS.CU_response;
    up_flag  = stACALOS.up_flag;
    down_flag = stACALOS.down_flag;

    % SE 07.05.2013 11:51
    useMeasuredCU_phase3 = stACALOS.useMeasuredCU_phase3;
end

% Level limitations (i.e. level range compression) due to response 50
%      in phases 2+ apply immediately and in next phase(s) (approved by T.Brand).
%      The resulting running limit increases by 5 dB after each phase
%      (like in literature).

% SE 07.05.2013 10:01 second point unclear, Level range is immediately applied (l 229)

% get current phase
if phase == 1
    % if first presentation level was audible proceed to phase 2
    if CU_response(end) ~= 50 && CU_response(end) ~= 0
        phase = 2;
        upper_limit_found = 0;
        lower_limit_found = 0;
    end
elseif phase == 2
    upper_limit_found = 0;
    lower_limit_found = 0;

    % if upper limit was already found proceed with the lower limit
    idx = find(up_flag == 1);
    for ii = 1:length(idx)
        if CU_response(idx(ii)) == 50 || level_history(idx(ii)) >= stACALOS.maximal_presentation_level;
            upper_limit_found = 1;
        end
    end

    % if lower limit was already found proceed with the upper limit
    idx = find( down_flag == 1);
    for ii = 1:length(idx)-1
        % DO 2013-07-17 Two options for finding the lower limit
        % 1) response CU = 0 and next response is CU > 0
        % 2) response CU > 0 and minimal_presentation_level reached
        if (CU_response(idx(ii)) == 0 && CU_response(idx(ii+1)) > 0) || (CU_response(idx(ii)) >  0 && level_history(idx(ii)) <= stACALOS.minimal_presentation_level)
            lower_limit_found = 1;
        end
    end

    % if lower and upper limit were found in phase 2, proceed with phase 3
    if upper_limit_found == 1 && lower_limit_found == 1
        phase = 3;
    end
end

% start phase - find audible level
if phase == 1
    level = level_history(end);
    % DO: 30.7.2014, increase level only by 5 dB of we are above 90 dB HL
    if CU_response(end) == 0; % if level wasn't audible
        if level >= 90
            new_level = level + 5;
        else
            new_level = level + 15;
        end
        if new_level > stACALOS.maximal_presentation_level
            stACALOS.new_level = stACALOS.maximal_presentation_level; % []
            stACALOS.failed = true;
            warning('could not find an audible level for this subject');
        end
        up_flag = [up_flag 0];
        down_flag = [down_flag 0];

    elseif CU_response(end) == 50; % if level was to loud

        new_level = level - 15;
        if new_level < stACALOS.minimal_presentation_level
            stACALOS.new_level = stACALOS.minimal_presentation_level; % []
            stACALOS.failed = true;
            warning('could not find an audible level for this subject');

        end
        up_flag = [up_flag 0];
        down_flag = [down_flag 0];

    end

end

% second phase - estimation of dynamic range
if phase == 2
    % DO: 2013-07-15 enable random range estimation
    if stACALOS.randomizeRangeEstimation > 0
        % randomize upper and lower track
        if rand(1) > 0.5
            try_upper_limit = true;
        else
            try_upper_limit = false;
        end

        % if up or down track is finished only present not-finished track
        if lower_limit_found
            try_upper_limit = true;
        end
        if upper_limit_found
            try_upper_limit = false;
        end
    end

    % find upper limit when randomizeRangeEstimation is true and ACALOS
    % should try upper limit OR if the randomizedRangeEstimation is false
    % do an alternating between the upper and the lower track
    if ( stACALOS.randomizeRangeEstimation > 0 && try_upper_limit ) || ( stACALOS.randomizeRangeEstimation == 0 && ((up_flag(end) == 0 && upper_limit_found==0) || lower_limit_found) )
        idx = find(up_flag == 1);
        if isempty(idx)
            % DO: 2013-07-15, fixed: wrong starting level handling
            idx_startlevel = max(find(~(down_flag|up_flag)));
            level = level_history(idx_startlevel);
        else
            level = level_history(idx(end));
        end
        
        
        if level + 10 >= 90 
            % for levels > 90 dB increase by 5 dB            
            new_level = level + 5;
        else
            % increase by 10 dB for leves below 90 dB
            new_level = level + 10;
        end

        % check if new level is below maximum presentation level
        if new_level >= stACALOS.maximal_presentation_level
            stACALOS.new_level = stACALOS.maximal_presentation_level;
            disp('maximal presentation level reached');
        end

        up_flag = [up_flag 1];
        down_flag = [down_flag 0];

    else % find lower limit
        idx = find(down_flag == 1);

        if isempty(idx) % if only starting level was tested until now
            % DO: 2013-07-15, fixed: wrong starting level handling
            idx_startlevel = max(find(~(down_flag|up_flag)));
            level = level_history(idx_startlevel);
        else
            level = level_history(idx(end));
        end

        % if there is no first occurrence of "not heard" in the lower track
        % decrease level
        if isempty(find(CU_response(idx)==0,1))
            new_level = level - 15;
        else
            % increase level to detect an audible response again
            new_level = level + 5;
        end

        up_flag = [up_flag 0];
        down_flag = [down_flag 1];
    end
end


% now we found a lower and upper limit and can present levels between them
if phase == 3
    idx = find(down_flag==1);
    % assume to find the level for 5 CU
    L5 = level_history(idx(end));
    % SE 07.05.2013 09:49 read out true CU value for that level
    CU_low = CU_response(idx(end));

    idx = find(up_flag == 1);
    % assume to find the level for 50 CU
    L50 = level_history(idx(end));
    % SE 07.05.2013 09:49 read out true CU value for that level
    CU_high = CU_response(idx(end));

    % SE changes in OMA 2013 according to Thomas Wittkop:
    % if we want to use measured CU values don't fall back to the
    % assumption that they were 50 and 5
    if ( ~useMeasuredCU_phase3 )
        % Assume that extreme responses were 50 and 5
        CU_high = 50;
        CU_low = 5;
    end

    % calculate the CU range
    CU_range = CU_high - CU_low;
    stACALOS.dynamic_range = L50 - L5;
    stACALOS.L5 = L5;
    stACALOS.L50 = L50;
    
    % estimate L15, L25, L35, L45 at beginning of phase 3
    if isempty( stACALOS.levels_phase3 )
        % linear interpolation of L15, L25, L35 and L45
        % L5 is neglected in the first iteration, see Brand, Hohmann (2002)
        x = [15 25 35 45];
        m = ((L50-L5)/CU_range);
        b = L5-m*CU_low;
        levels_phase3 = m*x + b;
        
        % counter for number of iterations
        iteration_phase3 = 1;
        % counter for number of level to present (four levels for
        % iteration_phase3 = 1 and five levels for iteration_phase3
        %
        % start with first level
        cur_presentation = 1;
        % maximal presentation level is set to level where 50 CU are
        % assumed
        L_max = L50;
        stACALOS.levels_phase3 = levels_phase3;
        stACALOS.L_max = L50;
    % check if we still have to generate levels
    elseif stACALOS.iteration_phase3 - 1 <= stACALOS.iterations
        cur_presentation = stACALOS.cur_presentation;
        iteration_phase3 = stACALOS.iteration_phase3;
        levels_phase3 = stACALOS.levels_phase3;
        L_max = stACALOS.L_max;
    % procedure has finished
    else
        stACALOS.finished = true;
        stACALOS.new_level = [];
        return
    end

    % apply new L_max limits from previous response immediately
    if iteration_phase3 + cur_presentation > 2
        % check if last response was 50 CU ("extremely loud")
        if CU_response(end) == 50
            % set new maximal presentation level
            L_max = level_history(end);
        elseif level_history(end) >= L_max && CU_response(end) < 50
            % increase L_max if it wasn't judged "extremely loud" again
            L_max = L_max + 5;
            % but do not increase L_max beyond given limits
            L_max = min(L_max,stACALOS.maximal_presentation_level);
        end
    end    
    
    if iteration_phase3 == 1 && cur_presentation == 1
        %% randomize order (try to prevent level differences > 0.5* dynamic range)
        % DO, 26.8.2014, use pseudo randomize order
        % two conditions should be met for the pseudo randomized order:
        % 1) the index distance to the next index should be 1 (up or down)
        % or 2 (up or down) for all steps
        % 2) the index distance should not be 1 for all entries (avoid
        % monotonous level increase or decrease)
        idx = randperm(length(levels_phase3));
        iteration = 0;
        while ~(sum(abs(diff(idx))<=2) == length(levels_phase3)-1) || ~(sum(abs(diff(idx)) == 1) < length(levels_phase3)-1)
            iteration = iteration + 1;
            idx = randperm(length(levels_phase3));
            if iteration == 10000
                % if no pseudo randomized order could be found after 10000
                % iterations, use the increasing order and issue a waring
                idx = 1:length(levels_phase3);
                warning('Levels in phase 3 cannot be arranged in pseudo randomized order.');
                break;
            end
        end   
        levels_phase3 = levels_phase3(idx);
    elseif iteration_phase3 > 1 && cur_presentation == 1;
        % estimate new levels after each new iteration
        % use all available responses to estimate new levels
        fitparams_start =[0.5,30];

        % estimate linear loudness function
        options = optimset('MaxFunEvals', 1000);
        fitparams = fminsearch(@(fitparams)costfc2(fitparams,stACALOS),fitparams_start,options);
        
        % estimate levels for next iteration
        m = fitparams(1);
        b = -fitparams(1)*fitparams(2);
        levels_phase3 = [5 15 25 35 45]./m - b/m;

        
        % limiting level to L_max before rearranging levels to meet
        % distance criterion
        % Dirk Oetting, 12.2.2014

        % if the maximum estimated level is above the maximum allowed level
        % L_max (adaptive during the measurement, at the beginning 105 dB HL)
        % we limit the presentation level to L_max and distribute
        % the other levels in the remaining dynamic range by preserving the
        % minimum estimated level
        if max(levels_phase3) > L_max 
            m = 40/(L_max  - min(levels_phase3));
            b = 5 - m*min(levels_phase3);
            levels_phase3 = [5 15 25 35 45]./m - b/m;
        end
        
        % if the lowest estimated level is below the minimal
        % presentation level, we limit the levels to the minimal
        % presentation level and rearrange the levels so that the maximum
        % level is preserved
        if min(levels_phase3) < stACALOS.minimal_presentation_level
            m = 40/(max(levels_phase3) - stACALOS.minimal_presentation_level);
            b = 5 - m*stACALOS.minimal_presentation_level;
            levels_phase3 = [5 15 25 35 45]./m - b/m;
        end
        
        % addional check for levels above L_max
        idx = levels_phase3 > L_max;
        levels_phase3(idx) = L_max;

        
        %% randomize order (try to prevent level differences > 0.5* dynamic range)
        % DO, 26.8.2014, use pseudo randomize order
        % two conditions should be met for the pseudo randomized order:
        % 1) the index distance to the next index should be 1 (up or down)
        % or 2 (up or down) for all steps
        % 2) the index distance should not be 1 for all entries (avoid
        % monotonous level increase or decrease)
        idx = randperm(length(levels_phase3));
        iteration = 0;
        while ~(sum(abs(diff(idx))<=2) == length(levels_phase3)-1) || ~(sum(abs(diff(idx)) == 1) < length(levels_phase3)-1)
            iteration = iteration + 1;
            idx = randperm(length(levels_phase3));
            if iteration == 10000
                % if no pseudo randomized order could be found after 10000
                % iterations, use the increasing order and issue a warning
                idx = 1:length(levels_phase3);
                warning('Levels in phase 3 cannot be arranged in pseudo randomized order.');
                break;
            end
        end   
        levels_phase3 = levels_phase3(idx);
    end
    
    % select new level
    new_level = levels_phase3(cur_presentation);

    % if current presentation was the last presentation of this iteration
    % step forward to next iteration
    cur_presentation = cur_presentation + 1;
    if cur_presentation > length(levels_phase3)
        % SE 02.05.2013 14:50
        % end if we are in the last iteration phase
        if iteration_phase3 - 1 == stACALOS.iterations
            % but return with new_level for the very last measurement
            stACALOS.finished = true;
        end
        % switch to next iteration
        cur_presentation = 1;
        iteration_phase3 = iteration_phase3 + 1;
    end
    % write phase specific variables to struct;
    stACALOS.levels_phase3 = levels_phase3;
    stACALOS.cur_presentation = cur_presentation;
    stACALOS.iteration_phase3 = iteration_phase3;
    stACALOS.L_max = L_max;
    %fprintf('L_max: %.1f CU\n',L_max);
end

% further check to limit level to allowed range
new_level = min(max(stACALOS.minimal_presentation_level,new_level),stACALOS.L_max);

% save variables to struct;
stACALOS.up_flag = up_flag;
stACALOS.down_flag = down_flag;
stACALOS.phase = phase;
% DO 2013-07-15 round new level in 1 dB steps
% SE 22.07.2013 11:05 customized rounding
switch ( stACALOS.levelRounding )
    case 1
        stACALOS.new_level = round(new_level);
    case 2
        stACALOS.new_level = round(new_level*10)/10;
    otherwise
        % do not round
        stACALOS.new_level = new_level;
end
stACALOS.level_history = [level_history stACALOS.new_level];


function erry = costfc2(fitparams,stACALOS)
% function [x]=costfc2(fitparams,stACALOS)
% use function to minimize error between data points and fitted function


level = stACALOS.level_history;
cu = stACALOS.CU_response;

% remove all data points with CU = 0, these points influence the
% linear fit to end with a too shallow slope
% Dirk Oetting, 30.1.2014
idx_cu0 = find(cu==0);
cu(idx_cu0) = [];
level(idx_cu0) = [];
% end of modification, DO 30.1.2014


m = fitparams(1);
b = -fitparams(1)*fitparams(2);

cu_fit = m*level+b;

% set points, where the measured data is 50 and the fit function is bigger
% than 50 to the constant value 50, so that the difference cu-cu_fit = 0
cu_fit(cu==50 & (cu_fit > cu))   = 50;

% set points, where the measured data is 0 and the fit function is lower
% than 0 to the constant value 0, so that the difference cu-cu_fit = 0
cu_fit(cu==0 & (cu_fit < cu))    = 0;


% calculate the distance between cu and cu_fit
delta_y = cu_fit - cu;
erry = sum((delta_y).^2);


% eof
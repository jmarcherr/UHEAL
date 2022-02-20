function [results] = get_data(d,s)
%load relevant fil
load([d(s).name]);
% subject name
sub_id = d(s).name;

%% audiograms
if ~isempty(dataalm.aud)
    %get audiogram and age
    [age,aud_L,aud_R,aud_freq,gender] = get_aud(dataalm);
    age = dataalm.subinfo.age;
    gender = dataalm.subinfo.gender;
    CP = dataalm.subinfo.CP;
    HV = dataalm.subinfo.HV;
    
    %get stim ear
    if strcmp(dataalm.id,'UH099')
        stimear = 1;
    elseif isempty(dataalm.stim.ffr.ear)
        stimear = 1;
    else
        stimear = dataalm.stim.ffr.ear(1);
    end
    %plot stim ear
    if stimear ==1 %left ear
        aud = aud_L;
    else % right ear
        aud = aud_R;
    end
    
    age_sub = age; %log age
    gender_sub = gender;
    CP_sub = CP;
    HV_sub = HV;
else
    age_sub = nan;
    gender_sub = nan;
    CP_sub = nan;
    aud = nan(12,1);
    aud_freq = nan(12,1);
    stimear = 1;
    HV_sub = nan;
    
end
results.sub_id = sub_id;
results.age_sub = age_sub;
results.gender_sub = gender_sub;
results.CP_sub = CP_sub;
results.HV_sub = HV_sub;
results.aud = aud;
results.aud_freq = aud_freq;
results.rds = dataalm.rds;
results.nesi = dataalm.nesi;
results.tts = dataalm.tts;
results.ssq = dataalm.ssq;
results.acalos = dataalm.acalos;

%% MEMR

    if ~isempty(dataalm.memr)
        % MEMR
        if stimear == 1 % left stim ear
            % does memr exist for stim ear?
            if isfield(dataalm.memr,'L')
                reflex = dataalm.memr.L.reflex_ipsi.response;
                levels = dataalm.memr.L.reflex_ipsi.labels;
                f_center = dataalm.memr.L.reflex_ipsi.f_center;
                freq = dataalm.memr.L.reflex_ipsi.freq;
            else
                reflex = dataalm.memr.R.reflex_ipsi.response;
                levels = dataalm.memr.R.reflex_ipsi.labels;
                f_center = dataalm.memr.R.reflex_ipsi.f_center;
                freq = dataalm.memr.R.reflex_ipsi.freq;
                
            end
        else            % right stim ear
            % does memr exist for stim ear?
            if isfield(dataalm.memr,'R')
                reflex = dataalm.memr.R.reflex_ipsi.response;
                levels = dataalm.memr.R.reflex_ipsi.labels;
                f_center = dataalm.memr.R.reflex_ipsi.f_center;
                freq = dataalm.memr.R.reflex_ipsi.freq;
            else
                reflex = dataalm.memr.L.reflex_ipsi.response;
                levels = dataalm.memr.L.reflex_ipsi.labels;
                f_center = dataalm.memr.L.reflex_ipsi.f_center;
                freq = dataalm.memr.L.reflex_ipsi.freq;
            end
        end
        f_lim = [6:15];
        [m_fc,I] = max((reflex(f_lim,end)'));
        growth_sub = reflex(I+f_lim(1)-1,:);
        growth_sub_alt = sum(abs(reflex));
        reflex_sub = reflex;
        MEM_slope = polyfit(levels,growth_sub_alt,1);        
    else
        warning('No MEMR data for subject')
        reflex = nan(16,7);
        levels = nan;
        MEM_slope = nan;
        growth_sub_alt = nan;
        reflex_sub = reflex;
        f_center = nan;
        freq = nan;
    end

    
    % get values
    results.memr.reflex_sub = reflex;
    results.memr.growth_sub = growth_sub_alt;
    results.memr.MEM_slope = MEM_slope;
    results.memr.levels = levels;
    results.memr.f_center = f_center;
    results.memr.freq = freq;
    %% TEOAE
    % Load scaped data
    if ~isempty(dataalm.teoae)
        teoae_amp = dataalm.teoae{1};
        teoae_resp = [dataalm.teoae{3},dataalm.teoae{2}];
        if isnan(teoae_amp(:,1,stimear))
            if stimear ==1
                stimear = 2;
            else
                stimear = 1;
            end
        end
        
        % save
        teoae_sub_amp = squeeze(teoae_amp(:,:,stimear));
        teoae_sub_resp = teoae_resp;
        results.teoae.amp = teoae_sub_amp;
        results.teoae.resp = teoae_sub_resp;
    else
        results.teoae.amp = nan(5,3);
        results.teoae.resp = nan(93,3);
    end
    clc
    disp(['Subject ' num2str(s) ' processed.'])
    

end





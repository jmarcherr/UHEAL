clear all;close all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
cd ..
cd ..
UHEAL_startup
ft_defaults

cd('_EEG/_preprocdata_AEP')
% subject
subdir = dir('UH*')

rates = [0.5,1,1.5,2];
%%
tic;
for s=21:length(subdir)
    
    % go to EEG folder
    cd(rootdir)
    cd('_EEG/_preprocdata_AEP')
    try
        load(subdir(s).name)
    catch
        load(subdir(s).name)
    end
    
    fs = data.fsample;
    %% ------------Comments --------------------------------------
    % AEP analysis
    % load stim file
    cd(datadir)
    cd(subdir(s).name(1:4))
    stimname = dir('AEP_stim_*');
    try
        load(stimname.name)
    catch
        load(stimname.name)
    end
    % missing trials
    if ~isempty(data.missing_trials)
        ids = stim.id(setdiff(1:length(stim.id),data.missing_trials));
    else
        ids = stim.id;
    end
    
    isi = stim.isi;
    
    
    
    %chan oi
    if strcmp(subdir(s).name(1:4),'UH20') % noisy FC4(16 channel) use
        chaoi = [2,4,7,8,9,10,12,14];
    else
        chaoi = [2,4,7,8,9,10,12,14,16];
    end
    chans =chaoi;
    
    
    for kk=[1:4] % condition loop
        
        trials_oi =find(ids==kk);
        
        cfg = [];
        cfg.channel = {data.label{[chans]}};%,data.label{[2]}, data.label{5}};
        data_cond = ft_selectdata(cfg,data);
        time = data.time{1};
        
        %tidx = find(time>=-.05 & time<.5); % chirps -1 to 8
        tidx = find(time>=0 & time<=.5);
        epoched_data=[];
        epoched_data =epoch_data_3(data_cond,1,trials_oi); %epoch x chan x time
        
        % artifact rejection and weighting
        reject_epoch = 0;
        reject_idx = 0;
        
        rjt_trials = [];
        for ii=1:size(epoched_data,1) % loop over trials
            for cc= 1:size(epoched_data,2) %loop over channels
                data_art = epoched_data(ii,cc,tidx);
                thr = 80;
                if max(abs(data_art)) > thr
                    rjt_trials = [rjt_trials ii];
                    %plot(squeeze(data_art))
                    %pause
                end
            end
        end
        rjt = unique(rjt_trials);
        trials_oi = [trials_oi];
        valid_trials = setxor(rjt,1:length(trials_oi));
        
        data_cc = epoched_data(valid_trials,:,:);
        data_cc = permute(data_cc,[2,3,1]);
        
        clc
        nr_reject = length(rjt)/(length(trials_oi));
        fprintf('%.2f %% of trials rejected! \n',nr_reject*100)
        %% weighted average
        fs = data.fsample;
        time = data.time{1};
        var_tmp =[];data_weighted=[];epoch_var=[];
        for cc=1:size(data_cc,1)
            for ii=1:size(data_cc,3) %chan x time x trials
                var_tmp(cc,ii) =var(data_cc(cc,tidx,ii));
                data_weighted(cc,:,ii) = data_cc(cc,tidx,ii)./var_tmp(cc,ii);
                epoch_var(cc,ii) = var_tmp(cc,ii).^(-1);
                
            end
        end
        
        catdata = [];
        for cc=1:size(data_cc,1)
            
            catdata =[catdata squeeze(data_weighted(cc,:,:))]; %cat chans
        end
        
        epoch_var = epoch_var(:);
        summed_trials =sum(catdata,2);
        summed_weights = sum(epoch_var);
        data_w = summed_trials./summed_weights;
        
        data_prefilt(s,kk,:) = data_w;
        data_trials{kk} = catdata./epoch_var';
        
        %     % HP filtering for EFR extraction
        filt_coef = [80 1000]; %4 Hz
        filt_def = designfilt('bandpassiir','FilterOrder',4, ...
            'HalfPowerFrequency1',filt_coef(1),'HalfPowerFrequency2',filt_coef(2), ...
            'SampleRate',fs,'designmethod','butter');
        
        data_filt(s,kk,:) =filtfilt(filt_def,data_w);%data_w;
        vts{s,kk} = trials_oi(valid_trials);
    end
    
    % timelock analysis
    subtrials = vts(s,:);
    % ids
    [fid,sid] = sort(cell2mat(subtrials));
    cfg = [];
    cfg.trials = fid;
    cfg.keeptrials  = 'yes';
    cfg.lpfilttype  = 'firws';
    cfg.lpfilter    = 'yes';
    cfg.lpfreq      = 30;
    
    
    
    data_cond = ft_preprocessing(cfg,data_cond);
    cfg = [];
    cfg.keeptrials  = 'yes';
    timelock{s} = ft_timelockanalysis(cfg,data_cond)
    timelock{s}.trialids = ids(fid);
    
    
    
    
    
    %% save processed AEP
    data_aep = struct;
    data_aep.aep = squeeze(data_filt(s,:,:));
    data_aep.trials = data_trials;
    data_aep.time = data.time{1};
    data_aep.tidx = tidx;
    data_aep.fs = fs;
    data_aep.subid = subdir(s).name(1:4);
    data_aep.nr_reject = nr_reject*100;
    data_aep.timelock = timelock{s}
    data_aep.subtrials = subtrials
    
    
    cd(rootdir)
    cd('_EEG/_AEP_results')
    %cd(subdir(s).name)
    savefile = [subdir(s).name(1:4) '_AEP_processed.mat'];
    
    save(savefile,'-struct','data_aep','-v7.3');
    
    
    
end
%%
cd('_timelock')
savefile = ['timelock_aep2.mat'];
    
save(savefile,'timelock','-v7.3');





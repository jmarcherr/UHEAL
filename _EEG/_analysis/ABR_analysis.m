clear all;close all
addpath('/work1/jonmarc/UHEAL_master/UHEAL')
UHEAL_startup
ft_defaults
eeg_dir = '/work1/jonmarc/UHEAL_master/UHEAL/_EEG/_preprocdata_ABR';
cd(eeg_dir)
% subject
subdir = dir('UH*')
%CP = [1,1];
%%
rates = [9,40];
for s=21%:length(subdir)
    % Select preprocessed .mat file to process
    data_filt = [];
    cd(datadir);
    cd(subdir(s).name);
    %stim 
    stimname = dir('click_abr_stim_*');
    try
       load(stimname.name);
    catch
       load(stimname.name);
    end
    stimear = stim.ear;
    % go to EEG folder
    cd(rootdir)
    cd('_EEG/_preprocdata_ABR')
    cd(subdir(s).name)
    try
        load([subdir(s).name,'_ABR.mat'])
    catch
        load([subdir(s).name,'_ABR.mat'])
    end


    cd(rootdir)
    %labels = {(1:18)};{'lMast','rMast','Cz','Fz','Fp1','FCz','ltip','rtip'};
    %data.label = labels(1:size(data.trial{1},1));
    fs = data.fsample;
    %% ------------Comments --------------------------------------
    %P_ABR:
    %1 = click, 80 nHL, 9.1Hz
    %2 = click, 80 nHL, 40Hz
    
    %chan oi
    if stimear ==1
        chans =[find(strcmp(data.label,'EXG1'))]; % right tiptrode
    else
        chans =[find(strcmp(data.label,'EXG2'))];
    end
    
    ccc = 1:2;
    for kk=[ccc] % condition loop
        trials_oi =find(data.trialinfo==kk); % condition trials
        
        cfg = [];
        cfg.channel = {data.label{[chans]}};
        data_cond = ft_selectdata(cfg,data); % select data
        time = data.time{1};
        
        tidx = find(time>=-5e-3 & time<15e-3); % chirps -1 to 15 ms
        
        epoched_data =epoch_data_3(data_cond,1,trials_oi); %epoch x chan x time
        
        % artifact rejection and weighting
        reject_epoch = 0;
        reject_idx = 0;
        
        rjt_trials = [];
        for ii=1:size(epoched_data,1) % loop over trials
            for cc= 1:size(epoched_data,2) %loop over channels
                data_art = epoched_data(ii,cc,tidx);
                thr = 20;
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
        if length(valid_trials)<3000
            data_cc = nan;
        else
        %data_cc{tt} = zscore(epoched_data(valid_trials,:)')';
        data_cc = epoched_data(valid_trials,:,:);
        data_cc = permute(data_cc,[2,3,1]);
        end
        data.valid_trials = valid_trials;
        
        clc
        nr_reject = length(rjt)/(length(trials_oi));
        fprintf('%.2f %% of trials rejected! \n',nr_reject*100)
        %% weighted average
        fs = data.fsample
        time = data.time{1};
        %tidx = find(time>=0 & time<.5);
        var_tmp =[];data_weighted=[];epoch_var=[];
        if length(valid_trials)<3000
            var_tmp = nan;
            data_weighted=nan;
            epoch_var=nan;
        else
            for cc=1:size(data_cc,1)
                for ii=1:size(data_cc,3) %chan x time x trials
                    var_tmp(cc,ii) =var(data_cc(cc,tidx,ii));
                    data_weighted(cc,:,ii) = data_cc(cc,tidx,ii)./var_tmp(cc,ii);
                    epoch_var(cc,ii) = var_tmp(cc,ii).^(-1);
                    
                end
            end
        end
        %catdata = squeeze(zeros(size(data_weighted(1,:,:))));
        catdata = [];
        for cc=1:size(data_cc,1)

                catdata =[catdata -squeeze(data_weighted(cc,:,:))]; %left tiptrode (inverting channel)
        end
        
        epoch_var = epoch_var(:);
        summed_trials =sum(catdata,2);
        summed_weights = sum(epoch_var);
        data_w = summed_trials./summed_weights;
        data_trials = catdata./epoch_var';
        
        %     % HP filtering
        
        %if strcmp('UH31',subdir(s).name)
            data_filt(s,kk,:) = data_w;
        %else
        %filt_coef = [80 3000]; %4 Hz
        %filt_def = designfilt('bandpassiir','FilterOrder',4, ...
        %    'HalfPowerFrequency1',filt_coef(1),'HalfPowerFrequency2',filt_coef(2), ...
        %    'SampleRate',fs,'designmethod','butter');
        %data_filt(s,kk,:) = filtfilt(filt_def,data_w);
        %end
        
        data_trials_sub{s,kk,:,:} = data_trials;
    end
    %% save processed ABR
    data_abr = struct;
    data_abr.abr = squeeze(data_filt(s,:,:));
    data_abr.trials = squeeze(data_trials_sub(s,:,:));
    data_abr.time = data.time{1};
    data_abr.tidx = tidx;
    data_abr.fs = fs;
    data_abr.subid = subdir(s).name; 
    data_abr.nr_reject = nr_reject*100;
    data_abr.stimear = stimear(1);
    

    cd('_EEG/_ABR_results')
    %cd(subdir(s).name)
    savefile = [subdir(s).name '_ABR_processed.mat'];

    save(savefile,'data_abr','-v7.3');
end


%% plotting
close all
figure(10)
delay = 0.60e-3%1.1e-3%
blIDX = round(delay/(1/fs));
for s=1:length(subdir)

        ccc = 1:2;
    
    for kk=ccc
        data_w = squeeze(data_filt(s,kk,:));
        baseline = data_w(find(time(tidx)==0)+blIDX);%
        
        
        tn = time(tidx)-delay;
        tnIDX = find(tn>-1e-3);
        subplot(6,5,s)

        sub_corrected(s,kk,:) =data_w(tnIDX)-baseline;
        p1(kk)=plot(tn(tnIDX),data_w(tnIDX)-baseline);
        hold on
        plot(tn(tnIDX),zeros(size(tn(tnIDX))),'k--');
        ylabel('[\mu V]');
        xlabel('Time [s]')
        box off
        set(gca,'Fontsize',10)
        
        clear colormap
        title(subdir(s).name(3:4))
        %cm=cbrewer('seq','Reds',200);
        %cm(1:10,:)=repmat([1 1 1],10,1);
        %colormap(cm)
        hold on
        
        
        xlim([-.5e-3 6e-3])
        ylim([-.5 .7])
        
    end
    
end
hleg = legend([p1(1) p1(2)],'9 Hz','40 Hz')
hleg.Box = 'off'
set(gcf,'position',[26   220   527   585]);
hold off


%% group average

figure(1000)
plot(tn(tnIDX),squeeze(mean(sub_corrected,1)));
hold on
plot(tn(tnIDX),zeros(size(tn(tnIDX))),'k--');
ylabel('ABR Amplitude [\mu V]');
xlabel('Time [s]')
box off
set(gca,'Fontsize',16)
hleg = legend([p1(1) p1(2)],'9 Hz','40 Hz')
hleg.Box = 'off'
set(gcf,'position',[457   527   269   214]);
hold off
        xlim([-.5e-3 6e-3])
        ylim([-.2 .4])




function [data_out] = epoch_data_3(data_in,trial_id,nr_epochs)
%Epoch data

%More info following

if ~isstruct(data_in)
    error('Data not in right format');
else
    clc
    disp('Epoching data');
end
if length(trial_id)>1
    idx = find(data_in.trialinfo==trial_id(1) | data_in.trialinfo==trial_id(2));
else
    idx=find(data_in.trialinfo==trial_id);
end
if length(nr_epochs)==1
    idx_tmp=idx(1:nr_epochs);
else
    idx_tmp=nr_epochs;
end

% generate epochs
epoched_data = [];
chans = size(data_in.trial{1},1);
times = size(data_in.trial{1},2);
epoched_data = reshape(cell2mat(data_in.trial(idx_tmp)),chans,times,length(idx_tmp)); % epoch x chan x time
clc
disp(['trials done.'])


%if size(epoched_data,1)>1 %mean over channels
%    epoched_data = mean(epoched_data,2);
%end
data_out = permute(epoched_data,[3,1,2]);

end



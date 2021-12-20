clear all;close all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
cd ..
cd ..
UHEAL_startup
ft_defaults

% subjects
subjects =setdiff(4:46,[20,31]);%
cd('_EEG/_preprocdata_FFR');
subdir = dir('UH*');
  
%%
for s=setdiff(31:length(subdir),[2,3,30])
    % Select preprocessed .mat file to process
    cd(datadir);
    cd(subdir(s).name(1:4));
    %stim
    stimname = dir('ffr_SW_stim_*');
    try
        load(stimname.name);
    catch
        load(stimname.name);
    end
    
    stimear = stim.ear;
    % go to EEG folder
    cd(rootdir)
    cd('_EEG/_preprocdata_FFR')

    try
        load([subdir(s).name])
    catch
        load([subdir(s).name])
    end


    cd(rootdir)
    %labels = {(1:18)};{'lMast','rMast','Cz','Fz','Fp1','FCz','ltip','rtip'};
    %data.label = labels(1:size(data.trial{1},1));
    fs = data.fsample;

    %% ------------------------ Preprocessing------------------------
    cd(rootdir)
    fs = data.fsample;
    
    %conditions:
    %1) 326 @ 4Hz
    %2) 326 @ 4Hz inv

    %%
    %chan oi
    if s==20  % noisy FC4 channel
        %chans =1:15;
    else
    chans =1:16%[find(strcmp(data.label,'Cz')) find(strcmp(data.label,'Fz'))]; % right tiptrode
    end
    cco = [1:2];
    for kk=1 % condition loop
        
        trials_oi =find(data.trialinfo==cco(kk));
        trials_inv = find(data.trialinfo==cco(kk)+1); %inverted
        % are there differences in trial nr?
        if length(trials_oi)~=length(trials_inv)
            [a,b]=min([length(trials_oi) length(trials_inv)]);
            trials_oi = trials_oi(1:a);
            trials_inv = trials_inv(1:a);
        end
        
        
        cfg = [];
        cfg.channel = {data.label{[chans]}};
        data_cond = ft_selectdata(cfg,data);
        time = data.time{1};
        
        tidx = find(time>=0 & time<0.5); % pt
        
        % epoch data
        epoched_data =epoch_data_3(data_cond,1,trials_oi); %epoch x chan x time
        epoched_inv = -epoch_data_3(data_cond,1,trials_inv); %inverted
        epoched_data = [epoched_data+epoched_inv]/2;
        
        % artifact rejection and weighting
        reject_epoch = 0;
        reject_idx = 0;
        
        rjt_trials = [];
        for ii=1:size(epoched_data,1) % loop over trials
            for cc= 1:size(epoched_data,2) %loop over channels
                data_art = epoched_data(ii,cc,tidx);
                thr = 25;
                if max(abs(data_art)) > thr
                    rjt_trials = [rjt_trials ii];
                    %plot(squeeze(data_art))
                    %pause
                end
            end
        end
        rjt = unique(rjt_trials)
        valid_trials = setxor(rjt,1:length(trials_oi)-1);
        data_cc = epoched_data(valid_trials,:,:);
        data_cc = permute(data_cc,[2,3,1]);
        data.valid_trials = valid_trials;
        
        clc
        nr_reject = length(rjt)/(length(trials_oi));
        fprintf('%.2f %% of trials rejected! \n',nr_reject*100)
        %% weighted average
        fs = data.fsample
        time = data.time{1};
        tidx = find(time>=0 & time<.5);
        var_tmp =[];data_weighted=[];epoch_var=[];
        for cc=1:size(data_cc,1)
            for ii=1:size(data_cc,3) %chan x time x trials
                var_tmp(cc,ii) =var(data_cc(cc,tidx,ii));
                data_weighted(cc,:,ii) = data_cc(cc,tidx,ii)./var_tmp(cc,ii);
                epoch_var(cc,ii) = var_tmp(cc,ii).^(-1);
                
            end
        end
        
        catdata = [];%squeeze(zeros(size(data_weighted(1,:,:))));
        for cc=1:size(data_cc,1)
            catdata =[catdata squeeze(data_weighted(cc,:,:))];
        end
        epochs_x_trial = 4;
        % Create a Trial linking epochs_x_trial (default: 16) epochs
        num_full_trials = floor(size(catdata,2)/epochs_x_trial);
        sample_x_epoch = size(catdata,1);
        num_chans = length(chans);
        catdata = reshape(catdata(:,...
            1:epochs_x_trial*num_full_trials), sample_x_epoch*epochs_x_trial,...
            num_full_trials);
        
        % Vector of the summed weights for each trial
        %epoch_var= reshape(epoch_var(1,1:num_full_trials*num_chans),...
        %    2, num_full_trials);
        
        epoch_var = epoch_var(:);
        summed_trials =sum(catdata,2);
        summed_weights = sum(epoch_var);
        data_w = summed_trials./summed_weights;
        
        % BP filtering
        filt_coef = [250 1000]; %4 Hz
        filt_def = designfilt('bandpassiir','FilterOrder',4, ...
            'HalfPowerFrequency1',filt_coef(1),'HalfPowerFrequency2',filt_coef(2), ...
            'SampleRate',fs,'designmethod','butter');
        
        data_filt(s,kk,:) = data_w; %weighted data per condition
        if isnan(data_filt(s,kk,1))
            disp('nan')
        end
    end
%end

%%
% ------------------------ Analysis ------------------------
foi = [326];
time = 0:1/fs:length(data_filt(1,1,:))/fs-1/fs;
%for s=1:length(subdir)
    for kk=1
        fid = foi(kk);
        tmp_data = squeeze(data_filt(s,kk,:));
        M=tmp_data;
        f_fft = [];
        %FFT
        f_fft = fft(M)/(length(M)/2);
        %Convert to power
        f_fft_pow(s,kk,:) = abs(f_fft.^2); %
        %Truncate negative freqencies
        fft_sub = squeeze((f_fft_pow(s,kk,1:end/2+1)));
        %Frequency vector
        f = fs/2*linspace(0,1,length(fft_sub));
        % get complex value at signal bin
        f_fft_cmplx(kk) = f_fft(find(f==fid));
        %get powerbin
        f_fft_sub_pow(kk) = fft_sub(find(f==fid));
        
        % get noise
        noisebw = 20;
        nonfreqs = [(fid+4) (fid-4)]';
        linenoise = [300]'
        f_fft_noise(kk) = mean(fft_sub(find([f>=(fid-noisebw) & ...
            f<=(fid+noisebw) & f~=fid & f~=nonfreqs(1) & f~=nonfreqs(2) & f~=linenoise])));
        
        %F-statistic
        F(s,kk)=f_fft_sub_pow(kk)/f_fft_noise(kk);
        bg_freq = find([f>=fid-noisebw & f<=fid+noisebw & f~=fid]);
        F_crit(s,kk) = finv(0.99,2,2*length(bg_freq));
        
        %SNR
        SNR(s,kk) = db(f_fft_sub_pow(kk))-db(f_fft_noise(kk));
        FFR(s,kk) =f_fft_sub_pow(kk);
        %plotting
        figure(kk)
        subplot(2,1,1)
        semilogx(f,db(fft_sub),'k');
        xlim([200 1050])
        hold on
        set(gca,'fontsize',18)
        xlabel('Freq Hz')
        ylabel('FFR amplitude (dB rel. 1uV)')
        set(gcf,'position',[680   828   353   270])
        
        subplot(2,1,2)
        hold on
        plot(time,tmp_data)
        set(gca,'fontsize',18)
        
        
    end
    
    %% Save FFR
    data_ffr = struct;
    data_ffr.f_fft     = fft_sub;
    data_ffr.fft_freq  = f;
    data_ffr.FFR       = FFR(s,:);
    data_ffr.FFR_SNR   = SNR(s,:);
    data_ffr.TS        = squeeze(data_filt(s,:,:));
    data_ffr.time      = data.time{1};
    data_ffr.F         = F(s,:);
    data_ffr.F_crit    = F_crit(s,:);
    data_ffr.tidx      = tidx;
    data_ffr.fs        = fs;
    
    data_ffr.subid = subdir(s).name(1:4); 
    data_ffr.nr_reject = nr_reject*100;
    data_ffr.stimear = stimear(1);
    

    cd('_EEG/_FFR_results/all_chans')
    %cd(subdir(s).name)
    savefile = [subdir(s).name(1:4) '_FFR_processed_allchans.mat'];

    save(savefile,'data_ffr','-v7.3');
    cd(rootdir)
end

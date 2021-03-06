
  %% EFR analysis script for HPC
function ABR_40hz_analysis(s,subdir,datadir,rootdir,clin_dir)
try
    % Select preprocessed .mat file to process
    cd(datadir);
    %cd(subdir(s).name(1:4));
    %stim

    cd(subdir(s).name(1:5));
    stimname = dir('ffr_SW_stim_*');
    
    % get stim ear
    this_stim=load(stimname.name);
    stimear = this_stim.stim.ear;

    % load subject info
    cd(clin_dir)
    load([subdir(s).name(1:5) '.mat']); % dataalm    
    
    % go to EEG folder
    cd(rootdir)
    cd('_EEG/_preprocdata_ABR/_40Hz')
    % load preprocessed data
    load([subdir(s).name])
    cd(rootdir)
    fs = data.fsample;

    %% ------------------------ Preprocessing------------------------
    cd(rootdir)

    
    %conditions:
    %1) 40Hz click EFR


    %%
    %chan oi (all channels)
    chans =1:length(data.label);%[find(strcmp(data.label,'Cz')) find(strcmp(data.label,'Fz'))]; % right tiptrode
    chan_labels = data.label(chans);
    cco = [1:2];
    for kk=1 % condition loop
        
        trials_oi =find(data.trialinfo==cco(kk));

        
        
        cfg = [];
        cfg.channel = {data.label{[chans]}};
        data_cond = ft_selectdata(cfg,data);
        time = data.time{1};
        
        % define time idx
        tidx_FFR = find(time>=0 & time<1);
        tidx = tidx_FFR;
        
        % epoch data
        epoched_data =epoch_data_3(data_cond,1,trials_oi); %epoch x chan x time

        
        % artifact rejection and weighting
        reject_epoch = 0;
        reject_idx = 0;
        
        rjt_trials = {};
        % artefact rejection on all data (with noise)
        for ii=1:size(epoched_data,1) % loop over trials
            for cc= 1:size(epoched_data,2) %loop over channels
                rjt_chan = [];
                data_art = epoched_data(ii,cc,tidx);
                thr = 25;
                if max(abs(data_art)) > thr
                    rjt_chan = [rjt_chan ii];
                    %plot(squeeze(data_art))
                    %pause
                end
            end
            rjt_trials{cc} = rjt_chan;
        end
        rjt = [unique([rjt_trials{:}])];
        valid_trials = setxor(rjt,1:length(trials_oi));
        data_cc = epoched_data(valid_trials,:,:);
        data_cc = permute(data_cc,[2,3,1]);
        data.valid_trials = valid_trials;
        
        clc
        nr_reject = length(rjt)/(length(trials_oi));
        fprintf('%.2f %% of trials rejected! \n',nr_reject*100)
        %% weighted average for FFR + noise +concat
        % FFR and noise
            fs = data.fsample
            time = data.time{1};
            %tidx = find(time>=0 & time<.5);
            var_tmp =[];data_weighted=[];epoch_var=[];
            for cc=1:size(data_cc,1)

                for ii=1:size(data_cc,3) %chan x time x trials
                    % EFR
                    var_tmp(cc,ii) =var(data_cc(cc,tidx_FFR,ii));
                    data_weighted(cc,:,ii) = data_cc(cc,tidx_FFR,ii)./var_tmp(cc,ii);
                    epoch_var_chan(cc,ii) = var_tmp(cc,ii).^(-1);           
                end  
                            % EFR single trial
                epoch_var_trial = epoch_var_chan(cc,:);
                summed_trial =sum(data_weighted(cc,:,:),3);
                summed_weights = sum(epoch_var_trial);
                data_w_trial(cc,:) = summed_trial./summed_weights;

                % for concatenated trials
                % Group epochs of FFR and noise together to make 1.5s trials
                catdata = [];
                catnoise = [];
                % How many epochs to group
                epochs_x_trial = 1;

                % Create a Trial linking epochs_x_trial
                num_full_trials = floor(size(squeeze(data_weighted(cc,:,:)),2)/epochs_x_trial);
                
                sample_x_epoch = size(data_weighted(cc,:,:),2);
                
                % Reshape
                tmp_data = squeeze(data_weighted(cc,:,:));
                % EFR
                catdata = reshape(tmp_data(:,1:epochs_x_trial*num_full_trials),...
                    sample_x_epoch*epochs_x_trial,...
                    num_full_trials);        
                
                % Vector of the summed weights for each trial
                epoch_var= reshape(epoch_var_chan(cc,1:num_full_trials*epochs_x_trial),...
                    epochs_x_trial, num_full_trials);

                
                % EFR
                epoch_var = epoch_var(:);
                summed_trials =sum(catdata(:,:),2);
                summed_weights = sum(epoch_var);
                data_w(cc,:) = summed_trials./summed_weights;
                

            
            end
            
            % BP filtering
            %         filt_coef = [250 2000];
            %         filt_def = designfilt('bandpassiir','FilterOrder',4, ...
            %             'HalfPowerFrequency1',filt_coef(1),'HalfPowerFrequency2',filt_coef(2), ...
            %             'SampleRate',fs,'designmethod','butter');

            data_filt = data_w; %weighted data per condition
            data_trial = data_w_trial;

            if isnan(data_filt)
                disp('nan')
            end

    end

%%
% ------------------------ Analysis ------------------------
foi = [38];
time = 0:1/fs:length(data_filt(1,:))/fs-1/fs;
%get fft EFR

[f,fft_sub,f_fft_noise,FFR,F,SNR,F_crit]=get_fft(data_filt,foi,fs);


    %% Save FFR
    data_efr = struct;
    data_efr.f_fft     = fft_sub;
    data_efr.fft_freq  = f;
    data_efr.FFR       = FFR;
    data_efr.FFR_SNR   = SNR;
    data_efr.TS_cat    = data_filt;
    data_efr.TS        = data_trial;
    %data_ffr.trials    = data_cc;
    data_efr.time      = data.time{1};
    data_efr.F         = F;
    data_efr.F_crit    = F_crit;
    data_efr.tidx      = tidx;
    data_efr.fs        = fs;
    data_efr.fft_noise = f_fft_noise;
    %data_efr.noise.f   = f_noise;
    %data_efr.noise.fft_sub_noise = fft_sub_noise;
    
    data_efr.subid = subdir(s).name(1:5); 
    data_efr.nr_reject = nr_reject*100;
    data_efr.rjt_trials = rjt_trials;
    data_efr.stimear = stimear(1);
    data_efr.channels = chans;
    data_efr.chan_labels = chan_labels;
    data_efr.subinfo = dataalm.subinfo;
catch % no eeg
    data_efr = {};
    data_efr.subid = subdir(s).name(1:5);
    data_efr.subinfo = dataalm.subinfo;
end
    cd(rootdir)
    cd('/work1/jonmarc/UHEAL_master/UHEAL/_EEG/_derivatives/_40hz_results')
    %cd(subdir(s).name)
    savefile = [subdir(s).name(1:5) '_40hz_processed.mat'];

    save(savefile,'data_efr','-v7.3');
    cd(rootdir)
end


function [f,fft_sub,f_fft_noise,FFR,F,SNR,F_crit]=get_fft(data,foi,fs)
chans = size(data,1);
for cc=1:chans  
    fid = foi;
    tmp_data = squeeze(data(cc,:));
    M=tmp_data;
    f_fft = [];
    %FFT
    f_fft = fft(M)/(length(M)/2);
    %Convert to power
    f_fft_pow(cc,:) = abs(f_fft.^2); %
    %Truncate negative freqencies
    fft_sub(cc,:) = squeeze((f_fft_pow(cc,1:end/2+1)));
    %Frequency vector
    f = fs/2*linspace(0,1,length(fft_sub(cc,:)));
    % get complex value at signal bin
    f_fft_cmplx = f_fft(find(f==fid));
    %get powerbin
    f_fft_sub_pow = fft_sub(cc,find(f==fid));
    
    %get noise
    noisebw = 10;
    nonfreqs = [];%[fid+2 fid-2]';%[(fid+4) (fid-4)]';
    linenoise = [100]';
    nbins = find([f>=(fid-noisebw) & ...
        f<=(fid+noisebw) & f~=fid...
                & f~=linenoise]);     

    f_fft_noise(cc,:) = mean(fft_sub(cc,nbins));
    
    %F-statistic
    F(cc)=f_fft_sub_pow/f_fft_noise(cc,:);
    bg_freq = nbins;%find([f>=fid-noisebw & f<=fid+noisebw & f~=fid]);
    F_crit(cc) = finv(0.99,2,2*length(bg_freq));
    
    %SNR
    SNR(cc) = db(f_fft_sub_pow)-db(f_fft_noise(cc,:));
    FFR(cc) =f_fft_sub_pow;

    
    
end

end
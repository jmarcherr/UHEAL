
%% FFR analysis script for HPC
function FFR_4Hz_analysis_par(s,subdir,datadir,rootdir,clin_dir)
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
    cd('_EEG/_preprocdata_FFR_4Hz')
    % load preprocessed data
    load([subdir(s).name])
    cd(rootdir)
    fs = data.fsample;
    
    %% ------------------------ Preprocessing------------------------
    cd(rootdir)
    
    
    %conditions:
    %1) 326 @ 4Hz
    %2) 326 @ 4Hz inv
    
    %%
    %chan oi (all channels)
    chans =1:length(data.label);%[find(strcmp(data.label,'Cz')) find(strcmp(data.label,'Fz'))]; % right tiptrode
    chan_labels = data.label(chans);


        
        trials_oi =find(data.trialinfo);
        
        
        cfg = [];
        cfg.channel = {data.label{[chans]}};
        data_cond = ft_selectdata(cfg,data);
        time = data.time{1};
        
        tidx = find(time>=0 & time<3.5);

        
        % epoch data
        epoched_data =epoch_data_3(data_cond,1,trials_oi); %epoch x chan x time

        
        % artifact rejection and weighting
        reject_epoch = 0;
        reject_idx = 0;
        

            epoched = epoched_data;
            rjt_trials = {};
            tidx = find(time>=0 & time<3.5);
            for ii=1:size(epoched,1) % loop over trials

                for cc= 1:size(epoched,2) %loop over channels
                    %       epoched = filtfilt(filt_def,squeeze(epoched(:,cc,:))');
                    
                    rjt_chan = [];
                    data_art = epoched(ii,cc,:);
                    thr = 110;
                    if max(abs(data_art)) > thr
                        rjt_chan = [rjt_chan ii];
                        %plot(squeeze(data_art))
                        %pause
                    end
                end
                rjt_trials{cc} = rjt_chan;
            end
            rjt = [unique([rjt_trials{:}])];
            valid_trials = setxor(rjt,1:length(trials_oi)-1);
            data_cc = epoched(valid_trials,:,:);
            data_cc = permute(data_cc,[2,3,1]);
            data.valid_trials = valid_trials;
            
            clc
            nr_reject = length(rjt)/(length(trials_oi));
            fprintf('%.2f %% of trials rejected! \n',nr_reject*100)
            %% weighted average
            fs = data.fsample
            time = data.time{1};
            tidx = find(time>=0 & time<3);
            var_tmp =[];data_weighted=[];epoch_var=[];
            for cc=1:size(data_cc,1)
                for ii=1:size(data_cc,3) %chan x time x trials
                    var_tmp(cc,ii) =var(data_cc(cc,:,ii));
                    data_weighted(cc,:,ii) = data_cc(cc,:,ii)./var_tmp(cc,ii);
                    epoch_var(cc,ii) = var_tmp(cc,ii).^(-1);
                    
                end
                epoch_var = epoch_var(cc,:);
                summed_trials =sum(data_weighted(cc,:,:),3);
                summed_weights = sum(epoch_var);
                data_w(cc,:) = summed_trials./summed_weights;
            end
            
            data_filt = data_w; %weighted data per condition
            if isnan(data_filt)
                disp('nan')
            end
            %data_filt = data_filt;
            
            % ITPC
        N=size(data_cc,3);
        fid = [4 8];
        for cc=1:size(data_cc,1) % channels
            for it = 1:N % trials
                
                M=data_cc(cc,tidx,it);
                %FFT
                f_fft = fft(M)/(length(M)/2);
                %Convert to power
                pow(it,:) = abs(f_fft.^2); %
                %Truncate negative freqencies
                ft_sub(it,:) = (pow(it,1:end/2+1));
                %Frequency vector
                f = fs/2*linspace(0,1,length(ft_sub(it,:)));
                % get complex value at signal bin
                itpc(it,:) = f_fft(1:end/2+1);
            end
            
            itpc = itpc./abs(itpc);
            itpc      = sum(itpc);   % sum angles
            itpc      = abs(itpc)/N;   % take the absolute value and normalize
            itpc_sub(cc,:)      = squeeze(itpc);
        end
            
            %%
            % ------------------------ Analysis ------------------------
       

            time = 0:1/fs:length(data_filt(1,:))/fs-1/fs;
            for cc=1:length(chans)
                
                tmp_data = squeeze(data_filt(cc,:));
                M=tmp_data;
                f_fft = [];
                %FFT
                f_fft = fft(M)/(length(M)/2);
                %Convert to power
                f_fft_pow(cc,:) = abs(f_fft.^2); %
                %Truncate negative freqencies
                fft_sub(cc,:) = squeeze((f_fft_pow(cc,1:end/2+1)));
                
            end


       

    %% Save FFR
    data_ffr = struct;
    data_ffr.f_fft     = fft_sub;
    data_ffr.itpc      = itpc_sub;
    data_ffr.fft_freq  = f;
    data_ffr.TS        = data_filt;
    %data_ffr.trials    = data_cc;
    data_ffr.time      = data.time{1};
    data_ffr.tidx      = tidx;
    data_ffr.fs        = fs;
    
    data_ffr.subid = subdir(s).name(1:5); 
    data_ffr.nr_reject = nr_reject*100;
    data_ffr.rjt_trials = rjt_trials;
    data_ffr.valid_trials = data.valid_trials
    data_ffr.stimear = stimear(1);
    data_ffr.channels = chans;
    data_ffr.chan_labels = chan_labels;
    data_ffr.subinfo = dataalm.subinfo;
catch % no eeg
    data_ffr = {};
    data_ffr.subid = subdir(s).name(1:5);
    data_ffr.subinfo = dataalm.subinfo;
end
    cd(rootdir)
    cd('/work1/jonmarc/UHEAL_master/UHEAL/_EEG/_derivatives/_4Hz_results')
    %cd(subdir(s).name)
    savefile = [subdir(s).name(1:5) '_4Hz_processed.mat'];

    save(savefile,'data_ffr','-v7.3');
    cd(rootdir)
end

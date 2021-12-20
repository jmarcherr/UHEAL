clear all;close all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
cd .. 
cd ..
UHEAL_startup
ft_defaults

cd('_EEG/_preprocdata_FFR_4Hz')
% subject
subdir = dir('UH*')

%%
for s=1:length(subdir)
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
    cd('_EEG/_preprocdata_FFR_4Hz')

    try
        load([subdir(s).name])
    catch
        load([subdir(s).name])
    end


    cd(rootdir)

    fs = data.fsample;
    
    chaoi = [2,3,4,7,8,9,10,12,14,16];
    chans =chaoi;%3,4,5,7,8,9,12,14,15,16]; %left tiptrode
    cco=1;
    for kk=1 % condition loop
        
        trials_oi =find(data.trialinfo==1);
        trials_inv = find(data.trialinfo==2); %inverted
        
        if length(trials_oi)~=length(trials_inv)
            [a,b]=min([length(trials_oi) length(trials_inv)]);
            trials_oi = trials_oi(1:a);
            trials_inv = trials_inv(1:a);
        end
        cfg = [];
        cfg.channel = {data.label{[chans]}};%,data.label{[2]}, data.label{5}};
        data_cond = ft_selectdata(cfg,data);
        time = data.time{1};
        
        tidx = find(time>=0 & time<0.5); % pt
        
        % epoch data
        epoched_data = epoch_data_3(data_cond,1,trials_oi); %epoch x chan x time
        epoched_inv = epoch_data_3(data_cond,1,trials_inv); %inverted
        epoched_data = [epoched_data+epoched_inv]/2;
        
        % artifact rejection and weighting
        reject_epoch = 0;
        reject_idx = 0;
        
        rjt_trials = [];
        for ii=1:size(epoched_data,1) % loop over trials
            for cc= 1:size(epoched_data,2) %loop over channels
                data_art = epoched_data(ii,cc,tidx);
                thr = 110;
                if max(abs(data_art)) > thr
                    rjt_trials = [rjt_trials ii];
                    %plot(squeeze(data_art))
                    %pause
                end
            end
        end
        rjt = unique(rjt_trials);
        valid_trials = setxor(rjt,1:length(trials_oi));
        data_cc = epoched_data(valid_trials,:,:);
        data_cc = permute(data_cc,[2,3,1]);
        data.valid_trials = valid_trials;
        
        clc
        nr_reject = length(rjt)/(length(trials_oi));
        fprintf('%.2f %% of trials rejected! \n',nr_reject*100)
        %% weighted average
        fs = data.fsample;
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
                % itpc
        N=size(catdata,2);
        fid = [4 8];
        for ff=1:2
            for it = 1:N
                
                M=catdata(:,it);
                %FFT
                f_fft = fft(M)/(length(M)/2);
                %Convert to power
                pow(it,:) = abs(f_fft.^2); %
                %Truncate negative freqencies
                ft_sub(it,:) = (pow(it,1:end/2+1));
                %Frequency vector
                f = fs/2*linspace(0,1,length(ft_sub(it,:)));
                % get complex value at signal bin
                itpc(it,:) = f_fft(find(f==fid(ff)));
            end
            
            itpc = itpc./abs(itpc);
            itpc      = sum(itpc);   % sum angles
            itpc      = abs(itpc)/N;   % take the absolute value and normalize
            itpc_sub(s,ff)      = squeeze(itpc);
        end
        
        
        epochs_x_trial = 16;
        % Create a Trial linking epochs_x_trial (default: 16) epochs
        num_full_trials = floor(size(catdata,2)/epochs_x_trial);
        sample_x_epoch = size(catdata,1);
        catdata = reshape(catdata(:,...
            1:epochs_x_trial*num_full_trials), sample_x_epoch*epochs_x_trial,...
            num_full_trials);
        

        
        % weighted average
        % Vector of the summed weights for each trial
        %epoch_var= reshape(epoch_var(1,1:num_full_trials*2),...
        %    2, num_full_trials);
        epoch_var = mean(epoch_var,1);
        epoch_var = epoch_var(:);
        summed_trials =sum(catdata,2);
        summed_weights = sum(epoch_var);
        data_w = summed_trials./summed_weights;
        
        
        data_filt(s,kk,:) = data_w; %weighted data per condition
        
    end




%% ------------------------ Analysis ------------------------
close all
foi = [4,8];
%itpc_sub=[];
for ff =1:2
    
    time = 0:1/fs:length(data_filt(1,1,:))/fs-1/fs;
    %for s=1:length(subjects)
        for kk=1
            fid = foi(ff);
            tmp_data = squeeze(data_filt(s,kk,:));
            M=tmp_data;
            N=size(tmp_data,1);
            f_fft = [];
            %FFT
            f_fft = fft(M)/(length(M)/2);
            %Convert to power
            f_fft_pow(s,ff,:) = abs(f_fft.^2); %
            %Truncate negative freqencies
            fft_sub = (squeeze(f_fft_pow(s,ff,1:end/2+1)));
            %Frequency vector
            f = fs/2*linspace(0,1,length(fft_sub));
            % get complex value at signal bin
            f_fft_cmplx = f_fft(find(f==fid));
            %get powerbin
            f_fft_sub_pow = fft_sub(find(f==fid));
            % spectrum
            fft_spec(s,ff,:) = fft_sub;
            
            % get noise
            noisebw =1.2;
            nonfreqs = [fid+4 fid-4]';
            f_fft_noise = mean(fft_sub(find([f>=fid-noisebw & f<=fid+noisebw & f~=fid & f~=nonfreqs])));
            
            %F-statistic
            F(s,ff)=f_fft_sub_pow/f_fft_noise;
            bg_freq = find([f>=fid-noisebw & f<=fid+noisebw & f~=fid]);
            F_crit(s,ff) = finv(0.99,2,2*length(bg_freq));
            
            %SNR
            SNR(s) = db(f_fft_sub_pow)-db(f_fft_noise)
            EFR(s) =f_fft_sub_pow;
            %plotting
            %figure(s)
            
%             %subplot(2,1,1)
%             if ff==1
%             %title(subjects{s})
%             semilogx(f,fft_sub)%,'k');
%             xlim([0 800])
%             hold on
%             set(gca,'fontsize',18)
%             xlabel('Freq Hz')
%             ylabel('EFR amplitude')
%             set(gcf,'position',[680   870   353   228]) 
%             
%             end
            %subplot(2,1,2)
            %hold on
            %plot(time,tmp_data)
            %set(gca,'fontsize',18)
            %xlabel('time(s)')
                   
            
        end
%        title(subjects{s})
%        hleg = legend('326 Hz');
%        hleg.Box = 'Off';
        % SNR and EFR
        SNR_f(s,ff) = SNR(s);
        EFR_f(s,ff) = EFR(s)
        %SNR_f(s,ff,2) = mean(SNR(s,[2]),2);
        
end
    
        %% Save FFR
    data_ffr = struct;
    data_ffr.f_fft     = fft_spec(s,:,:);
    data_ffr.fft_freq  = f;
    data_ffr.EFR = EFR_f(s,:);
    data_ffr.SNR = SNR_f(s,:);
    data_ffr.F   = F(s,:);
    data_ffr.F_crit = F_crit(s,:)
    data_ffr.foi = foi;
    data_ffr.data_w = squeeze(data_filt(s,kk,:));
    data_ffr.iptc_sub = itpc_sub(s,:);
    

    data_ffr.subid = subdir(s).name(1:4); 
    data_ffr.nr_reject = nr_reject*100;
    data_ffr.stimear = stimear(1);
    

    cd('_EEG/_FFR_4Hz_results')
    %cd(subdir(s).name)
    savefile = [subdir(s).name(1:4) '_FFR_4Hz_processed.mat'];

    save(savefile,'data_ffr','-v7.3');
    cd(rootdir)
end
%%

% %% mean spectrogram
% close all
% for s=1:4
%     for kk=1
%         figure(1)
%         subplot(2,2,s)
% y_mean = squeeze(mean(mean(f_fft_pow(s,1,kk,1:end/2+1),1),3));
% %o_mean = squeeze(mean(mean(f_fft_pow(find(~CP),2,[1 3],1:end/2+1),1),3));
% %semilogx(f,o_mean,'k')
% xlim([1 20])
% ylim([0 0.4])
% hold on
% semilogx(f,y_mean)
% xlabel('Frequency (Hz)')
% ylabel('EFR amplitude')
% set(gca,'fontsize',18,'xtick',[2:2:10])
% 
% 
% set(gcf,'position',[680   243   521   562])
% %hleg = legend('326 Hz')
% %hleg.Box = 'Off'
% box off
% title(subjects{s})
%     end
% end

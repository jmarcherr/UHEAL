
function FFR_par_full(datadir,d,dd)

%rng('default'); parcreaterandstream(10, randi(10000)+dd)
rootdir = cd;
cd(datadir)
cd(d(dd).name)
try
% get BDF name
bdf = dir('*.bdf')
if ~isempty(bdf) %|| strcmp(d(dd).name,'UH091') || strcmp(d(dd).name,'UH067') %no ABR
    dataset = bdf.name;
    % load stim file
    stim_file = [];
    stim_file = dir('ffr_SW_stim*');
    if ~isempty(stim_file)
        load(stim_file.name);
        % Get stim ear
        stimear = stim.ear(1);
        %cd(bdfdir)
    end
    
    % load dataalm for this participant
    cd ..
    cd('scraped')
    load([d(dd).name '.mat'])
    eeg_lab = dataalm.subinfo.lab;
    
    % back to datadir
    cd(datadir)
    cd(d(dd).name)
    %% ------------Event extraction --------------------------------------
    if strcmp(eeg_lab,'PHY2')
        triggers = [10,22]; % 9.1/s, 9.1/s+noise, 40/s
    else
        triggers = [10,20];
    end
    abort=0;
    try
        hdr = ft_read_header(dataset);
        cfg=[];
        cfg.layout =  'biosemi64.lay';
        cfg.continuous = 'yes';
        cfg.channel     = 'all';%{'EXG1','EXG2','EXG3','EXG4','EXG5','EXG6','EXG7','EXG8', 'Status'};
        cfg.dataset = dataset;
        cfg.trialdef.eventtype    = 'STATUS';
        cfg.trialdef.eventvalue   = triggers;
        cfg.trialdef.prestim      = .1; %5 ms
        cfg.trialdef.poststim     = 3.5; %20 ms
        cfg = ft_definetrial(cfg);
    catch
        abort=1;
    end
    % TFS FFR dataset
    if ~abort
        for tt=1:length(triggers)
            cfg.trl(cfg.trl(:,4)==triggers(tt),4) = tt;
            trl_org = cfg.trl;
        end
        if size(cfg.trl,1)==972
        % take every 6 only
            %trl_tmp = cfg.trl(1:6:end,:);
            %cfg.trl = trl_tmp;
            missing_trials = 0;

        else
            warning('Subject has missing trials')

            missing_trials = 1;
        end
            % find 6 trial sequence
            if dd==31 || dd==3
                fs=2048;
            else
                fs=4096*4;
            end
          
            for ii=1:length(trl_org)-1
                this_stamp = trl_org(ii,1) % this zero point
                next_stamp(ii) = trl_org(ii+1,1)-this_stamp;
            end
            break_trials = find(next_stamp/fs>18000/fs);
            break_trials = [1 break_trials+1];
            trl_tmp = cfg.trl(break_trials,:);
            cfg.trl = trl_tmp;
            
        

        
        % initial preprocessing
        data_int = ft_preprocessing(cfg);
        fs = data_int.fsample;
        for ii=1:length(trl_tmp)-1           
            these_trials = break_trials(ii):break_trials(ii+1)-1
            this_break = break_trials(ii+1)-1
            time_stamps{ii}= (trl_org([these_trials],1)-trl_org(these_trials(1),1))/fs;
        end
        these_trials = break_trials(ii+1):length(trl_org)
        this_break = length(trl_org)
        time_stamps{ii+1}=(trl_org([these_trials],1)-trl_org(these_trials(1),1))/fs;
       
        % Resampling
        
        cfgres = [];
        cfgres.resamplefs = 4096;
        cfgres.detrend    = 'no';
        data = ft_resampledata(cfgres,data_int);
        %%
        %Rereferencing (l/r mastoid)
        %cfg = [];

        cfg.dataset = dataset;
        cfg.channel     = {'eeg','EXG1','EXG2' '-Status'};%chaoi;;%chaoi;
        cfg.reref       = 'yes';
        if stimear ==1
            cfg.refchannel  = {'EXG1'}; %vertex electrodes%linked mastoids
        else
            cfg.refchannel  = {'EXG2'};
        end
        %cfg.refchannel = {'Cz'};
        if dd == 2 % UH02 wrong recording labels
            if stimear ==1
                cfg.refchannel  = {'EXG1'}; %vertex electrodes%linked mastoids
            else
                cfg.refchannel  = {'EXG2'};
            end
            
            cfg.channel     = {'A1','A2','A3','A4','A5','A6','A7','A8','A9','A10',...
                'A11' 'A12' 'A13' 'A14' 'A15' 'A16','EXG1','EXG2' '-Status'};
        end
        cfg.layout      =  'biosemi64.lay';
        cfg.continuous  = 'yes';
        cfg.dftfilter   = 'yes';
        cfg.dftfreq     = [50 100  150];
        cfg.lpfilttype  = 'firws';
        cfg.lpfilter    = 'yes';
        if dd==31 || dd==3 % recorded with 2048 fs
            cfg.lpfreq      = 1000;%3000;
        else
            cfg.lpfreq      = 2000; 
        end
        cfg.hpfilter    = 'yes';
        cfg.hpfilttype  = 'firws';
        cfg.hpfreq      = 1;%80; % changed from 80

        
        
        % rereferenced data struct
        data = ft_preprocessing(cfg,data);
        data.missing_trials = missing_trials;
        data.trl_org = trl_org;
        data.time_stamps = time_stamps;
        data.break_trials = break_trials;
        if dd==2
            data.label = {'Fp1','F3','AFz','Fz','P9','T7','C3','FC3','C4','Cz',...
                'P10','FCz','T8','F4','Fp2','FC4','EXG1','EXG2'};
        end
        %
        %     cfgd          = [];
        %     cfgd.method   = 'channel';
        %     cfgd.channel = 'all'
        %     cfgd.viewmode = 'butterfly';
        %     ft_databrowser(cfgd,data)
        %     pause
        %clear data_int
        


        
        cd(datadir)
        cd ..
        cd(['_EEG' filesep '_preprocdata_FFR/full'])
        %%  Save mat
        %if ~exist(d(dd).name, 'dir')
        %     mkdir(d(dd).name)
        % %end
        %
        % cd(d(dd).name)
        % if stimear ==1
        %
        %     savefile = [dataset(1:2) '_ABR_ltip.mat']
        % else
        %     savefile = [dataset(1:2) '_ABR_rtip.mat']
        % end
        savefile = [d(dd).name '_FFR.mat'];
        
        save(savefile,'data','-v7.3');

            
        clear data_DG
    end


end
catch
end

cd(rootdir)




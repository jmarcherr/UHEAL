
function ABR_par(datadir,d,dd)

%rng('default'); parcreaterandstream(10, randi(10000)+dd)
rootdir = cd;
cd(datadir)
cd(d(dd).name)

% get BDF name
bdf = dir('*.bdf')
if ~isempty(bdf) || strcmp(d(dd).name,'UH091') || strcmp(d(dd).name,'UH067') %no ABR
    dataset = bdf.name;
    % load stim file
    stim_file = [];
    stim_file = dir('click_abr_stim*');
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
        triggers = [62]; % 9.1/s, 9.1/s+noise, 40/s
    else
        triggers = [60];
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
        cfg.trialdef.poststim     = 1.1; %20 ms
        cfg = ft_definetrial(cfg);
    catch
        abort=1;
    end
    if ~abort
        for tt=1:length(triggers)
            cfg.trl(cfg.trl(:,4)==triggers(tt),4) = tt;
        end
        cfg.trl = cfg.trl(1:40:end,:,:,:);
        
        %%
        %Rereferencing (l/r mastoid)
        %cfg = [];
        cfg.dataset = dataset;
        cfg.channel     = {'eeg','EXG1','EXG2' '-Status'};%chaoi;;%chaoi;
        cfg.reref       = 'yes';
        % if stimear ==1
        %     cfg.refchannel  = {'EXG1'}; %vertex electrodes%linked mastoids
        % else
        %     cfg.refchannel  = {'EXG2'};
        % end
        cfg.refchannel = {'P9','P10'};
        if dd == 2 % UH02 wrong recording labels
            cfg.refchannel = {'A5','A11'};
            cfg.channel     = {'A1','A2','A3','A4','A5','A6','A7','A8','A9','A10',...
                'A11' 'A12' 'A13' 'A14' 'A15' 'A16','EXG1','EXG2' '-Status'};
        end
        cfg.layout      =  'biosemi64.lay';
        cfg.continuous  = 'yes';
        cfg.dftfilter   = 'yes';
        cfg.dftfreq     = [50 100  150];
        cfg.lpfilttype  = 'but';
        cfg.lpfilter    = 'yes';
        cfg.lpfiltord   = 2;
        if dd==31 || dd==3 % recorded with 2048 fs
            cfg.lpfreq      = 1000;%3000;
        else
            cfg.lpfreq      = 3000;
        end
        cfg.hpfilter    = 'yes';
        cfg.hpfilttype  = 'but';
        cfg.hpfreq      = 10; % changed from 5
       %cfg.demean = 'yes';
        %cfg.detrend = 'yes';
        %cfg.method      ='channel'
        cfg.hpfiltord   = 4;
        
        
        % rereferenced data struct
        data = ft_preprocessing(cfg);
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
        
        % Resample to 16 kHz
        %if dd==31 || dd==3 % recorded with 2048 fs
            cfgres = [];
            cfgres.resamplefs = 2048;
            cfgres.detrend    = 'no';
            data = ft_resampledata(cfgres, data);
        %end
        
        cd(datadir)
        cd ..
        cd(['_EEG' filesep '_preprocdata_ABR' filesep '_40Hz'])
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
        savefile = [d(dd).name '_ABR_40Hz.mat'];
        
        save(savefile,'data','-v7.3');
        
        clear data_DG
    end
end

cd(rootdir)




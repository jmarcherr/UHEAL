
function ABR_par(datadir,d,dd)

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
        triggers = [50,62]; % 9.1/s, 9.1/s+noise, 40/s
    else
        triggers = [50,60];
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
        cfg.trialdef.prestim      = 10e-3; %5 ms
        cfg.trialdef.poststim     = 20e-3; %20 ms
        cfg = ft_definetrial(cfg);
    catch
        abort=1;
    end
    if ~abort
        for tt=1:length(triggers)
            cfg.trl(cfg.trl(:,4)==triggers(tt),4) = tt;
        end
        
        %%
        %Rereferencing (l/r mastoid)
        %cfg = [];
        cfg.dataset = dataset;
        cfg.channel     = {'eeg','EXG1','EXG2' '-Status'};%chaoi;;%chaoi;
        cfg.reref       = 'yes';
        cfg.refchannel = {'Cz','Fz','FCz'};
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
        cfg.hpfreq      = 100; % changed from 5
        cfg.demean = 'yes';
        %cfg.detrend = 'yes';
        %cfg.method      ='channel'
        cfg.hpfiltord   = 4;
        
        
        % rereferenced data struct
        data = ft_preprocessing(cfg);


        %     cfgd          = [];
        %     cfgd.method   = 'channel';
        %     cfgd.channel = 'all'
        %     cfgd.viewmode = 'butterfly';
        %     ft_databrowser(cfgd,data)
        %     pause
        %clear data_int
        
        % Resample to 16 kHz
        if dd==31 || dd==3 % recorded with 2048 fs
            cfgres = [];
            cfgres.resamplefs = 16384;
            cfgres.detrend    = 'no';
            data = ft_resampledata(cfgres, data);
        end
        
        cd(datadir)
        cd ..
        cd(['_EEG' filesep '_preprocdata_ABR'])
        %%  Save mat
        savefile = [d(dd).name '_ABR.mat'];
        
        save(savefile,'data','-v7.3');
        
        clear data_DG
    end
end

cd(rootdir)




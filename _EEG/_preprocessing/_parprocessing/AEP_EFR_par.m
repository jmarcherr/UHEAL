
function AEP_EFR_par(datadir,d,dd)

%rng('default'); parcreaterandstream(10, randi(10000)+dd)
rootdir = cd;
cd(datadir)
cd(d(dd).name)

% get BDF name
bdf = dir('*.bdf')
try
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
        triggers = [102]; % 9.1/s, 9.1/s+noise, 40/s
    else
        triggers = [100];
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
        cfg.trialdef.prestim      = .1; %100 ms
        cfg.trialdef.poststim     = .6; %ss ms
        cfg = ft_definetrial(cfg);
    catch
        abort=1;
    end
    if ~abort
        for tt=1:length(triggers)
            cfg.trl(cfg.trl(:,4)==triggers(tt),4) = tt;
        end

        %%
        data_int = ft_preprocessing(cfg);
        
        
        %%
        
        cfgres = [];
        cfgres.resamplefs = 4096;
        cfgres.detrend    = 'no';
        data = ft_resampledata(cfgres,data_int);
        %Rereferencing (l/r mastoid)
        %cfg = [];
        cfg.dataset = dataset;
        cfg.channel     = {'eeg','EXG1','EXG2' '-Status'};%chaoi;;%chaoi;
        cfg.reref       = 'yes';
         %if stimear ==1
         %    cfg.refchannel  = {'EXG1'}; %vertex electrodes%linked mastoids
         %else
         %    cfg.refchannel  = {'EXG2'};
         %end
        if stimear ==1
            cfg.refchannel  = {'P9','P10'}; %vertex electrodes%linked mastoids
        else
            cfg.refchannel  = {'P9','P10'};
        end
        
        %if dd == 2 % UH02 wrong recording labels
        %    
        %    cfg.channel     = {'A1','A2','A3','A4','A5','A6','A7','A8','A9','A10',...
        %        'A11' 'A12' 'A13' 'A14' 'A15' 'A16','EXG1','EXG2' '-Status'};
        %end
        cfg.layout      =  'biosemi64.lay';
        cfg.continuous  = 'yes';
        cfg.dftfilter   = 'yes';
        cfg.dftfreq     = [50 100  150];
        cfg.lpfilttype  = 'firws';
        cfg.lpfilter    = 'yes';
        cfg.lpfreq      = 2000;
        cfg.hpfilter    = 'yes';
        cfg.hpfilttype  = 'firws';
        cfg.hpfreq      = 80; % changed from 80

        
        
        % rereferenced data struct
        data = ft_preprocessing(cfg,data);
%         if dd==2
%             data.label = {'Fp1','F3','AFz','Fz','P9','T7','C3','FC3','C4','Cz',...
%                 'P10','FCz','T8','F4','Fp2','FC4','EXG1','EXG2'};
%         end
        %
        %     cfgd          = [];
        %     cfgd.method   = 'channel';
        %     cfgd.channel = 'all'
        %     cfgd.viewmode = 'butterfly';
        %     ft_databrowser(cfgd,data)
        %     pause
        %clear data_int
        
        % Resample to 4 kHz


        
        cd(datadir)
        cd ..
        cd(['_EEG' filesep '_preprocdata_AEP' filesep 'EFR' filesep 'mastoid'])
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
        savefile = [d(dd).name '_EFR.mat'];
        
        save(savefile,'data','-v7.3');
        
        clear data_DG
    end
    
end
catch
    warning(['Subject ' num2str(dd) 'not processed'])
end

cd(rootdir)




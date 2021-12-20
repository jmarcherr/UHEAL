clear all
ft_defaults
UHEAL_startup

%cd(bdfdir)
cd(datadir)
d = dir('UH*');
numsub = length(d);
 

%% ------------Import data ----------------------------------------
for dd=setdiff(38:length(d),[21,31])
cd(datadir)
cd(d(dd).name)

% get BDF name
bdf = dir('*.bdf')
dataset = bdf.name;
% load stim file
stim_file = dir('ffr_SW_stim*');
load(stim_file.name);
% Get stim ear
stimear = stim.ear(1);
%cd(bdfdir)

%comments:

%% ------------Event extraction --------------------------------------

triggers = [102];
hdr = ft_read_header(dataset);
cfg=[];
cfg.layout =  'biosemi64.lay';
cfg.continuous = 'yes';
    cfg.channel     = {'eeg','EXG1','EXG2', 'Status'};
cfg.dataset = dataset;
cfg.trialdef.eventtype    = 'STATUS';
cfg.trialdef.eventvalue   = triggers;
cfg.trialdef.prestim      = 0.05;
cfg.trialdef.poststim     = 0.5;
cfg = ft_definetrial(cfg);


for tt=1:length(triggers)
    cfg.trl(cfg.trl(:,4)==triggers(tt),4) = tt;
end

%% check trials
close all
fs = 16384;
tmp = (cfg.trl(2:end,:)-cfg.trl(1:end-1,:));
missing_trials = find(tmp(:,1)>=3.8e4);
plot(tmp(:,1)/fs)

% 
% stim.id(172) = [];
% stim.id(366) = [];
% hold on
% plot(stim.id)

%%
%inital preprocessing(all channels no reref)
data_int = ft_preprocessing(cfg);


%%

cfgres = [];
cfgres.resamplefs = 4096;
cfgres.detrend    = 'no';
data = ft_resampledata(cfgres,data_int);

%Rereferencing (l/r mastoid)
%cfg = [];
cfg.dataset = dataset;
cfg.channel     = {'eeg','EXG1','EXG2', '-Status'};%chaoi;
cfg.reref       = 'yes';
cfg.refchannel  = {'EXG1','EXG2'}; %vertex electrodes%linked mastoids
cfg.layout      =  'biosemi64.lay';
cfg.continuous  = 'yes';
cfg.dftfilter   = 'yes';
cfg.dftfreq     = [50:50:1000];
%cfg.dftbandwidth = [4*ones(size(50:50:1000))];
cfg.lpfilttype  = 'firws';
cfg.lpfilter    = 'yes';
%cfg.lpfiltord   = 2;
cfg.lpfreq      = 2000;
cfg.hpfilter    = 'yes';
cfg.hpfilttype  = 'firws';
cfg.hpfreq      = 0.5;
cfg.demean = 'yes';
%cfg.baselinewindow = [-0.1 0];
%cfg.method      ='channel'
%cfg.hpfiltord   = 2;


% rereferenced data struct
data = ft_preprocessing(cfg,data);

% 
%     cfgd          = [];
%     cfgd.method   = 'channel';
%     cfgd.channel = 'all'
%     cfgd.viewmode = 'butterfly';
%     ft_databrowser(cfgd,data_int)
    %%
%     pause
%clear data_int
% Resample to 4096Hz
% cfgres = [];
% cfgres.resamplefs = 4096;
% cfgres.detrend    = 'no';
% data = ft_resampledata(cfgres, data_int);
data.missing_trials = missing_trials;

cd(datadir)
cd ..
cd(['_EEG' filesep '_preprocdata_AEP'])
%%  Save mat


savefile = [d(dd).name '_AEPtips.mat'];
%cd('./chirp')
save(savefile,'data','-v7.3');

clear data_DG
cd(rootdir)
%end
end

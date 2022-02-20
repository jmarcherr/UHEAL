clear all
ft_defaults
UHEAL_startup

%cd(bdfdir)
cd(datadir)
d = dir('UH*');
numsub = length(d);
%d = dir('*.bdf');
 
%% ------------Import data ----------------------------------------
for dd=1:length(d)
cd(datadir)
cd(d(dd).name)

% get BDF name
bdf = dir('*.bdf')
dataset = bdf.name;
% load stim file
stim_file = dir('ffr_SW_stim*');
try
    load(stim_file.name);
catch
    load(stim_file.name);
end
% Get stim ear
stimear = stim.ear(1);
%cd(bdfdir)

%% ------------Event extraction --------------------------------------
triggers = [10,22]; % cmplx tone, up-down-chirp, up-down-chirp SNR 5d
%10 22 = ffr
%50 62 = ABR
%100 = BBAEP


hdr = ft_read_header(dataset);
cfg=[];
cfg.layout =  'biosemi64.lay';
cfg.continuous = 'yes';
cfg.channel     = 'all'%{'EXG1','EXG2','EXG3','EXG4','EXG5','EXG6','EXG7','EXG8', 'Status'};
cfg.dataset = dataset;
cfg.trialdef.eventtype    = 'STATUS';
cfg.trialdef.eventvalue   = triggers;
cfg.trialdef.prestim      = 0.1;
cfg.trialdef.poststim     = .5;
cfg = ft_definetrial(cfg);

for tt=1:length(triggers)
    cfg.trl(cfg.trl(:,4)==triggers(tt),4) = tt;
end



%%
%cfg = [];
cfg.dataset = dataset;
%p9,p10,cz,fz,fp1,fcz,tipL,tipR
cfg.channel     = {'eeg','EXG1','EXG2' '-Status'};%chaoi;
cfg.reref       = 'yes';
if stimear ==1
    cfg.refchannel  = {'EXG1'};%'EXG1'}; %Rtip
else 
    cfg.refchannel  = {'EXG2'};%'EXG1'}; %Rtip
end
if dd==3 || dd==31
    cfg.refchannel = {'P9'}; %no exg channels recorded
end
cfg.layout      =  'biosemi64.lay';
cfg.continuous  = 'yes';
cfg.dftfilter   = 'yes';
cfg.dftfreq     = [50:50:1000];
%cfg.dftbandwidth = [4*ones(size(50:50:1000))];
cfg.lpfilttype  = 'firws';
cfg.lpfilter    = 'yes';
%cfg.lpfiltord   = 2;
cfg.lpfreq      = 1000;
cfg.hpfilter    = 'yes';
cfg.hpfilttype  = 'firws';
cfg.hpfreq      = 80;
%cfg.method      ='channel'
%cfg.hpfiltord   = 2;


% rereferenced data struct
data_int = ft_preprocessing(cfg);

% Resample to 2048Hz
cfgres = [];
cfgres.resamplefs = 2048;
cfgres.detrend    = 'no';
data = ft_resampledata(cfgres, data_int);

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
cd(['_EEG' filesep '_preprocdata_FFR'])
%%  Save mat
% if ~exist(dataset(1:2), 'dir')
%     mkdir(dataset(1:2))
% end
% 
% cd(d(dd).name(1:2))

savefile = [d(dd).name '_FFR.mat']

save(savefile,'data','-v7.3');

clear data_DG
cd(rootdir)
%end
end

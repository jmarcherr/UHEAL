clear all

%addpath('O:\Public\Hearing-Systems-group\cahr\Temporary_ftp\UHEAL')
cd('/work1/jonmarc/UHEAL_master/UHEAL')
UHEAL_startup
ft_defaults

%cd(bdfdir)
cd(datadir)
d = dir('UH*');
numsub = length(d);
%d = dir('*.bdf');
tic
%% ------------Import data ----------------------------------------
for dd=2:10%1:10%setdiff(2:numsub,[21])

cd(datadir)
cd(d(dd).name)

% get BDF name
bdf = dir('*.bdf')
dataset = bdf.name;
% load stim file
stim_file = dir('click_abr_stim*');
load(stim_file.name);
% Get stim ear
stimear = stim.ear(1);
%cd(bdfdir)

%% ------------Event extraction --------------------------------------

triggers = [50,62]; % 9.1/s, 9.1/s+noise, 40/s

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
for tt=1:length(triggers)
    cfg.trl(cfg.trl(:,4)==triggers(tt),4) = tt;
end

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
cfg.refchannel = {'Cz','Fz','FCz'};
if dd == 2 % UH02 wrong recording labels
    cfg.refchannel = {'A10','A4','A12'};
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
if dd==31 || 3 % recorded with 2048 fs
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
if dd==31 || 3 % recorded with 2048 fs
cfgres = [];
cfgres.resamplefs = 16384;
cfgres.detrend    = 'no';
data = ft_resampledata(cfgres, data);
end

cd(datadir)
cd ..
cd(['_EEG' filesep '_preprocdata_ABR'])
%%  Save mat
%if ~exist(d(dd).name, 'dir')
    mkdir(d(dd).name)
%end

cd(d(dd).name)
% if stimear ==1
%     
%     savefile = [dataset(1:2) '_ABR_ltip.mat']
% else
%     savefile = [dataset(1:2) '_ABR_rtip.mat']
% end
savefile = [d(dd).name '_ABR.mat'];

save(savefile,'data','-v7.3');

clear data_DG
cd(rootdir)
%end
end
clc
toc
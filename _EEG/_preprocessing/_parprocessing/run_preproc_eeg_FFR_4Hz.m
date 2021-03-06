clear all

%addpath('O:\Public\Hearing-Systems-group\cahr\Temporary_ftp\UHEAL')
cd('/work1/jonmarc/UHEAL_master/UHEAL')
UHEAL_startup
ft_defaults
addpath('/work1/jonmarc/UHEAL_master/UHEAL/_EEG/_preprocessing/_parprocessing')

%cd(bdfdir)
cd(datadir)
d = dir('UH*');
numsub = length(d);

%% HPC cluster parameters
clust=parcluster('dcc');    % load the MDCS cluster profile
clust.AdditionalProperties.MemUsage = '8GB';
clust.AdditionalProperties.WallTime = '20:00';
clust.saveProfile;
%%
parpool(clust, 20);
% %%
%addpath('~/fieldtrip'),ft_defaults;

%%
cd('/work1/jonmarc/UHEAL_master/UHEAL/_EEG/_preprocessing/_parprocessing')
parfor kk = 1:length(d)
    FFR_par_4Hz(datadir,d,kk);
end

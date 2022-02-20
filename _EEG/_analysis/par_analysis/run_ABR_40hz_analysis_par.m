cd('/work1/jonmarc/UHEAL_master/UHEAL')
UHEAL_startup
ft_defaults
% dir with preproc data
eeg_dir = '/work1/jonmarc/UHEAL_master/UHEAL/_EEG/_preprocdata_ABR/_40Hz';
% dir with subject info
clin_dir = '/work1/jonmarc/UHEAL_master/UHEAL/UHEAL_data/scraped';
cd(eeg_dir)
% subject
subdir = dir('UH*')
%CP = [1,1];

%% HPC cluster parameters
clust=parcluster('dcc');    % load the MDCS cluster profile
clust.AdditionalProperties.MemUsage = '8GB';
clust.AdditionalProperties.WallTime = '20:00';
clust.saveProfile;
%%
parpool(clust, 20);
%%
addpath('/work1/jonmarc/UHEAL_master/UHEAL/_EEG/_analysis/par_analysis')
parfor kk = 1:length(subdir)
    ABR_40hz_analysis(kk,subdir,datadir,rootdir,clin_dir);
end
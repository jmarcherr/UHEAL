cd('/work1/jonmarc/UHEAL_master/UHEAL')
UHEAL_startup
ft_defaults
% dir with preproc data
eeg_dir = '/work1/jonmarc/UHEAL_master/UHEAL/_EEG/_preprocdata_ABR/mid_latencies';
% dir with subject info
clin_dir = '/work1/jonmarc/UHEAL_master/UHEAL/UHEAL_data/scraped';
cd(eeg_dir)
% subject
subdir = dir('UH*')
%CP = [1,1];

%% HPC cluster parameters
clust=parcluster('dcc');    % load the MDCS cluster profile
clust.AdditionalProperties.MemUsage = '10GB';
clust.AdditionalProperties.WallTime = '20:00';
clust.saveProfile;
%%
parpool(clust, 20);
%%
addpath('/work1/jonmarc/UHEAL_master/UHEAL/_EEG/_analysis/par_analysis')
for kk = 20%1:length(subdir)
    midlate_analysis_par(kk,subdir,datadir,rootdir,clin_dir);
end
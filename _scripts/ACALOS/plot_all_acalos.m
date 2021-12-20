%% plot all acalos

close all
addpath(fileparts(matlab.desktop.editor.getActiveFilename))
cd('O:\Public\Hearing-Systems-group\cahr\Temporary_ftp\UHEAL\UHEAL_data')
%addpath('O:\Public\Hearing-Systems-group\cahr\Temporary_ftp\UHEAL')
UHEAL_startup
cd(datadir)
cd([datadir '/scraped']);
freq_aud = [250 500 1000 2000 4000 8000 9000 10000 11200 12500 14000 16000];
d=dir('UH*.mat')
%% get acalos data
for s=1:length(d)
    
    %load relevant file
    %load(['UH' num2str(subid(s))])
    
    load([d(s).name]);
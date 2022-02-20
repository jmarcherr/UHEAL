%cd(fileparts(matlab.desktop.editor.getActiveFilename))
clear
cd /work1/jonmarc/UHEAL_master/UHEAL/
datadir = [cd filesep 'UHEAL_data'];
addpath([cd '/UHEAL_Scraper'])
addpath([cd '/_scripts/ACALOS'])
addpath([cd '/_scripts/_tools/ACALOS'])
datafol_remote = '/work1/jonmarc/UHEAL_master/UHEAL/UHEAL_data';
%%
% run scraper
UHEALscraper_par(datafol_remote)


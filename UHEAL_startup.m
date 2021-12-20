% startup
rootdir = fileparts(which('UHEAL_startup.m'));
datadir = ([rootdir, filesep 'UHEAL_data']);
%bdfdir = ([rootdir,'/_data/raw_bdf'])


try %#ok
    rng(1); 
end  


fprintf('\n project directory now added to the current path \n')

if ~exist(fileparts(which('ft_defaults.m')))
    fprintf('remember to add fieldtrip to you path! \n')
end


addpath(fullfile('_EEG/'))
%addpath(fullfile('_func/'))

addpath (fullfile('_EEG/_preprocessing'))
addpath (fullfile('_EEG/_func'))
addpath (fullfile('_EEG/_analysis'))
addpath (fullfile('_scripts/_tools/cbrewer/cbrewer'))
addpath (fullfile('_scripts'))



%addpath _data
%addpath _data/raw_bdf





fprintf('\n directory addded to the path')




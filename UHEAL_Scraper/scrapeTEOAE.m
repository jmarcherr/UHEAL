%This function is part of the UHEALscraper. It accepts the struct created
%by reading in the Interacoustics Titan TEOAE module .xml file by the
%function 'xml2struct.m' (Copyright (c) 2010, Wouter Falkena) and outputs a
%vastly simplified cell of the key TEOAE data only.

%if several attempts are recorded (multible cells under the 'RecordedData')
%field, then the last one is taken.

function TEOAE = scrapeTEOAE(structRAW)

TEOAEwf = nan(93,2); %initate storage for the avergage waveform ([waveform] x ear)
TEOAEstats = nan(5,3,2); %initiate storage freq x [sig, noise, passed/failed] x ear

structRAW = structRAW.SaData.Session; %strip off unrequired fields

LRperm = [1,2]; %default assumed order of storage of Left and Right ear data 
%check if assumption is true
if strcmp(structRAW.Test.Data.RecordedData{1,1}.Measured.EarSide.Text, 'Right')
    LRperm = [2,1]; %if not; switch order
end

for i = LRperm %loop across ears, starting with left
    
    %check if the data exists before trying to extract it
    if isfield(structRAW.Test.Data.RecordedData{1,i}.Measured.BandResults,'Band')
    
    %extract OAE waveform per ear x and y (chain of transformations required to extract properly)  
  %%!!! There appears here to be two choices, 'BufferA' and 'BufferB' for time series data. They both have Xaxis scales which run from 4 to 12.345.
  % OaeData, NoiseData, and ProbeData don't appear to be time series information.
        OAE = cell2mat(structRAW.Test.Data.RecordedData{1,i}.Measured.BufferA.Point); %or 'BufferB'
        
        OAEx = squeeze(struct2cell([OAE.X])); OAEy = squeeze(struct2cell([OAE.Y]));
        OAEx = cell2mat(cellfun(@str2num,OAEx,'un',0)); OAEy = cell2mat(cellfun(@str2num,OAEy,'un',0));
        
        TEOAEwf(:,i) = OAEy; %y-axis (pressure?) values of the OAE
    
        for j = 1:size(structRAW.Test.Data.RecordedData{1,i}.Measured.BandResults.Band,2) %loop across frequency bands
    
        %extracts key stats for saving
        sig = str2num(structRAW.Test.Data.RecordedData{1,i}.Measured.BandResults.Band{1,j}.Attributes.Signal);
        noise = str2num(structRAW.Test.Data.RecordedData{1,i}.Measured.BandResults.Band{1,j}.Attributes.Noise);
        PF = strcmp(structRAW.Test.Data.RecordedData{1,i}.Measured.BandResults.Band{1,j}.Attributes.Passed,'true');
        TEOAEstats(j,:,i) = [sig, noise, PF];       
                
        end
     
    
    
    end
end
    
TEOAE = {TEOAEstats, TEOAEwf, OAEx}; %save data into cell array (OAEx is the waveform time vector (assumed to be the same for L and R))

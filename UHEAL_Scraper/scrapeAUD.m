%This function is part of the UHEALscraper. It accepts the struct created
%by reading in the Interacoustics diagnostic suite .xml file by the
%function 'xml2struct.m' (Copyright (c) 2010, Wouter Falkena) and outputs a
%vastly simplified cell of the key audiogram data only.

%if several attempts are recorded (multible cells under the 'RecordedData'
%field), then the last one is taken.

function [AUD, error] = scrapeAUD(structRAW)

AUD = zeros(12,2,2); %initiate storage freq x [freq, HL thresh] x ear
freqlookup = [1,2,3,4,5,6,7,8,9,10,11,12 ; 250,500,1000,2000,4000,8000,9000,10000,11200,12500,14000,16000];

error = {}; %initiate flag to notify of an issue

structRAW = structRAW.SaData.Session; %strip off unrequired fields

for i = 1:size(structRAW.Test.Data.RecordedData.Measured,2)%loop across measurement cells (divides L/R and low/high freq if done with different headphones)
    for j = 1:size(structRAW.Test.Data.RecordedData.Measured{1,i}.Tone.TonePoint,2)%loop across contained frequencies
        
        freq = str2num(structRAW.Test.Data.RecordedData.Measured{1,i}.Tone.TonePoint{1,j}.Frequency.Text);
        [~,k] = find(freqlookup == freq); %check where that frequency should go (in case some are missing)
        
        %check if it is left or right ear data
        if strcmp(structRAW.Test.Data.RecordedData.Measured{1,i}.Tone.Earside.Text, 'Left')
            m = 1;
        elseif strcmp(structRAW.Test.Data.RecordedData.Measured{1,i}.Tone.Earside.Text, 'Right')
            m=2;
        end
        AUD(k,1,m) = freq; % save frequency        
        AUD(k,2,m) = str2num(structRAW.Test.Data.RecordedData.Measured{1,i}.Tone.TonePoint{1,j}.IntensityUT.Text); %extract and save HL threshold
               
    end    
end

 


if any(AUD(:,1,:)==0)  %any of AUD's frequency entries have remained as '0' generate an error message and report
    warning('This audiogram file has less than the full expected range of frequencies')
    fprintf('The scrapeAUD function will return this incomplete matrix:\n')
    AUD
      
   error = [error; [' was missing some audiogram frequencies.']]; 
   
   
end

if size(structRAW.Test.Data.RecordedData.Measured,2) < 4
    warning('It appears as this audiogram was not performed with different headphones for the low & high frequencies')
    
    error = [error; [' was not performed with different headphones for the low & high frequencies.']]; 
    
end



    


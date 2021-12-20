function [acalos_results] = get_acalos_raw(acalos_raw)
% get raw acalos results and create data structure
%Jonatan Marcher-Rorsted 17/05/21

HLcorrection = zeros(4,1);
Rawresults = acalos_raw;

validfreq = [500 1000 2000 4000];
x= [-10:2.5:120];
%delete data from 250 and 6000 Hz
Rawresults(find(Rawresults(:,1)==250 | Rawresults(:,1)==6000),:)=[];
freq = unique(Rawresults(:,1));

if sum(sum(validfreq==freq))~=length(validfreq)
    validData = sum(validfreq==freq); %% This indexes will be NAN data
else
    validData = [1 1 1 1];
end
data_tmp = zeros(length(freq),5);
% loop
for ii=1:length(validfreq)
    if validData(ii)
        iOK = Rawresults(:,1)==freq(ii);
        dataLoudFit = Rawresults(iOK,[2,3]);
        dataLoudFit(:,1) = dataLoudFit(:,1);
        %     dataLoudFit(:,1) = dataLoudFit(:,1);
        dataLoudFit = reshape(dataLoudFit', 1, []);
        
        select_fitting_function = 'BTUX';
        switch upper(select_fitting_function)
            case 'BY'
                % fitting method of Brand and Hohmann (2002)
                fitparameters = fit_loudness_function(dataLoudFit,'BY');
                %res{i}=[fit(1), fit(2), fit(3)];
            case {'BX','BTX','BTUX','BTUY','BTPX'}
                % fitting method of Oetting et al. (2014) with optional
                % threshold estimation and UCL estimation if not enough
                % reponse were measured in the upper loudness domain
                [fitparams,rawData] = fit_loudness_function(dataLoudFit,select_fitting_function);
                
                % calculate parameters as for the the BX fitting function
                %res{i}=fitparams;
                %resRaw{i}=[repmat(work.int_exppar1{i},size(rawData,1),1),rawData]';
                data_tmp(ii,1)=1; % RUN, Just to make it compatible
                
                data_tmp(ii,2)=freq(ii);
                data_tmp(ii,[3,4,5])=fitparams;
        end
        
        %NumRun = unique(Rawresults(:,1));
        %% Estimated parameters
        
        [y, failed] = loudness_function_bh2002(x, fitparams);
        iZero = (max(find(y==0)));
        iUCL = (min(find(y==50)));
        %%%%%%%%%%%%%%%%%%%%% HTL
        if  ~isempty(iZero)
            HTL(ii)=x(iZero+1);
        else
            if ~isempty(x(min(find(y<=1))+1))
                HTL(ii) =  x(min(find(y<=1))+1);
            else
                HTL(ii) = x(find(min(y)));
            end
        end
        %%%%%%%%%%%%%%%%%%%%% UCL
        if  ~isempty(iUCL)
            UCL(ii)=x(iUCL-1);
        else
            UCL(ii) = x(y==max(y));
        end
        
        Lcut =fitparams(1);
        m_lo = fitparams(2);
        %        HTL_A(ii) = Lcut - 22.5/m_lo;
        HTL_A(ii)= fitparams(1)-22.5/m_lo;
        MCL(ii) = x(max(find(y<=25)));
        
        %        %%%%%%%%%%%%%%%%%%%%%% OHC LOSS ESTIMATION %%%%%%%%%%%%%%%%
        %         if OHC_estimation
        %             if sum(freq==[500 1000])
        %                 [k_best,y_cu_model] = adapt_cu_model_to_data(x_HL,y,HTL(ii),'NB', freq);
        %                 OHL_LF_tmp= HTL(ii)*k_best;
        %                 OHL_LF = [OHL_LF OHL_LF_tmp];
        %             end
        %             if sum(freq==[2000 4000])
        %                 [k_best,y_cu_model] = adapt_cu_model_to_data(x_HL,y,HTL(ii),'NB', freq);
        %                 OHL_HF_tmp= HTL(ii)*k_best;
        %                 OHL_HF = [OHL_HF OHL_HF_tmp];
        %             end
        %         end
        %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        y_plot(:,ii) = y;
        x_plot(:,ii) = x;
    else
        
        y_plot(:,ii) = NaN*ones(1,length(x));
        x_plot(:,ii) = NaN*ones(1,length(x));
        
        UCL(ii) = NaN;
        
        
        Lcut =NaN;
        m_lo = NaN;%fitparams(2);
        
        HTL(ii)= NaN;
        MCL(ii) = NaN;
        data_tmp(ii,1)=1; % RUN, Just to make it compatible
        
        data_tmp(ii,2)=NaN;
        data_tmp(ii,[3,4,5])=NaN(1,3);
    end
end


Slope =data_tmp(:,4)';
Locut = data_tmp(:,3)';
m_high = data_tmp(:,5)';

%% Plot
% if exist('typePlot','var')
%     hh = [];
% else
%     hh=figure('visible','off');
% end
% hh=figure;
% hdata=[];
% col=[0 0 0 ; 0.5 0.5 0.5];
% hdata = plot(fliplr(x_plot),fliplr(y_plot));
% legend(fliplr({'0.25','0.5','1','2','4','6'}))


Results.Frequencies = validfreq;
Results.Lcut = Locut;
Results.m_low = Slope;
Results.m_high = m_high;
Results.HTL = HTL;
Results.MCL = MCL;
Results.UCL = UCL;
Results.RawData = Rawresults;

acalos_results = Results;
end
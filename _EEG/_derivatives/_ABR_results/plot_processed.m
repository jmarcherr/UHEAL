clear;close all;
cd(fileparts(matlab.desktop.editor.getActiveFilename))
% find all processed files
plotting=0
d = dir('*.mat');
clin_dir = '/work1/jonmarc/UHEAL_master/UHEAL/UHEAL_data/scraped';
% init
delay = 1.1e-3%1.1e-3% Check this delay
fs = 2^14;
blIDX = round(delay/(1/fs));
% conditions
ccc = 1:2;
%%
for s=1:length(d)%[26,11,28,61]%
    load([d(s).name])
    this_dir = cd;
    % load subject info
    cd(clin_dir);
    load([d(s).name(1:5) '.mat']);
    subinfo = dataalm.subinfo;
    cd(this_dir);
    % get age
    if isempty(subinfo.age)
        age(s)=nan;
    else
        age(s)=subinfo.age;
    end
    % get lab
    if strcmp(subinfo.lab,'PHY2')
        lab(s) = 1;
    elseif strcmp(subinfo.lab,'PHY1')
        lab(s) = 2;
    elseif strcmp(subinfo.lab,'CHBC')
        lab(s) = 3;
    else
        lab(s) = nan;
    end
    % get gender
    if isempty(subinfo.gender)
        gender(s)=nan;
    else
        gender(s)=subinfo.gender;
    end    
    
    cd(this_dir)
    %% Baseline correction
    %loop over rates (9,40)
    for kk=ccc
        data_w = data_abr.abr{kk};
        if length(data_w)<300
            data_w = zeros(327,1);
        end

        baseline = data_w(find(data_abr.time(data_abr.tidx)==0)+blIDX);%

        
        tn = data_abr.time(data_abr.tidx)-delay;
        tnIDX = find(tn>-1e-3);
        %figure(s)
        if plotting
            subplot(11,10,s)
        end
        if mean(data_w==0)
            sub_corrected(s,kk,:) = nan(size(sub_corrected,3),1);
        else
            
            sub_corrected(s,kk,:) =data_w(tnIDX)-baseline;
        end
        % plot all abr
        if plotting
            p1(kk)=plot(tn(tnIDX)*1000,data_w(tnIDX)-baseline);
            hold on
            plot(tn(tnIDX)*1000,zeros(size(tn(tnIDX))),'k--');
            box off
            set(gca,'Fontsize',10)            
            clear colormap
            title([d(s).name(3:5)])
            hold on     
            xlim([-.5e-3 8e-3]*1000)
            ylim([-.5 .7])
        end
    end
    clc
    disp(['proccesssing ' data_abr.subid(1:5)])
end
if plotting
    hleg = legend([p1(1) p1(2)],'9 Hz','40 Hz')
    hleg.Box = 'off'
    set(gcf,'position',[26   220   527   585]);
    hold off
end
%%
figure(1005)
subidx_preAD=(1:20);
subidx_postAD=(21:length(d));
p1=plot(tn(tnIDX)+0.00085,squeeze(nanmean(sub_corrected(subidx_preAD,:,:),1)));
hold on
p2=plot(tn(tnIDX),squeeze(nanmean(sub_corrected(subidx_postAD,:,:),1)));
plot(tn(tnIDX),zeros(size(tn(tnIDX))),'k--');
ylabel('ABR Amplitude [\mu V]');
xlabel('Time [s]')
box off
set(gca,'Fontsize',16)
hleg = legend([p1(1) p1(2) p2(1) p2(2)],'9 Hz','40 Hz','9 Hz after','40 Hz after')
hleg.Box = 'off'
%set(gcf,'position',[583 533 356 304]);
hold off
        xlim([-.5e-3 8e-3])%changed from 6e-3 to 8e-3
        ylim([-.2 .4])
 
        
        %% mean with delay
 close all
sub_tmp = sub_corrected(21:length(d),:,:);
sub_tmp=sub_tmp(:,:,14:end);
%sub_tmp = [sub_tmp zeros(42,2,14)];
 
sub_mean = nanmean([cat(1,sub_tmp,sub_corrected(1:20,:,1:244-13))])
 
plot(tn(tnIDX(14:end)),squeeze(sub_mean)')
        xlim([-.5e-3 8e-3])%changed from 6e-3 to 8e-3
        ylim([-.2 .4])
        hold on
        
        
 plot(tn(tnIDX),zeros(size(tn(tnIDX))),'k--');
ylabel('ABR Amplitude [\mu V]');
xlabel('Time [s]')       
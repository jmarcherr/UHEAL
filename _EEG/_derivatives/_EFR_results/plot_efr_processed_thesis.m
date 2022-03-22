%% plot EFR results
clear all;close all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
thisdir = cd;
cd('/work1/jonmarc/UHEAL_master/UHEAL')
addpath('/work1/jonmarc/UHEAL_master/UHEAL/_scripts/_tools')
UHEAL_startup
ft_defaults
cd(thisdir)
%subjects

thisdir = cd;
subs = dir('UH*');
load('/work1/jonmarc/UHEAL_master/UHEAL/UHEAL_data/scraped/uheal_data_table/uheal_data.mat');

%% get data
for s=1:length(subs)
    
    load(subs(s).name)
    clc
    disp(['sub ' subs(s).name(1:5) ' loaded...'])
    sub_num(s) = str2num(subs(s).name(3:5));
    % get FFR
    if isfield(data_efr,'FFR')
        FFR{s} = data_efr.FFR;
        FFR_SNR{s} = data_efr.FFR_SNR;
        F{s} = data_efr.F;
        F_crit{s} = data_efr.F_crit;
        subid{s} = data_efr.subid;
        f_fft{s} = data_efr.f_fft;
        fft_freq(s,:) = data_efr.fft_freq;
        TS{s} = data_efr.TS;
        %trials{s} = data_ffr.trials;
        time = data_efr.time;
        tidx = data_efr.tidx;
        fft_noise{s} = data_efr.fft_noise;
        noisef{s} = data_efr.noise.fft_sub_noise;
        stimear(s,:) = data_efr.stimear;
        subinfo{s} = data_efr.subinfo;
        if isempty(data_efr.subinfo.age)|isempty(data_efr.subinfo.gender)
            age(s) = nan;
            gender(s) = nan;
        else
        age(s) =data_efr.subinfo.age;
        gender(s) = data_efr.subinfo.gender;
        end

        if strcmp(data_efr.subinfo.lab,'PHY2')
            lab(s) = 1;
        elseif strcmp(data_efr.subinfo.lab,'PHY1')
            lab(s) = 2;
        elseif strcmp(data_efr.subinfo.lab,'CHBC')
            lab(s) = 3;
        else
            lab(s) = nan;
        end
        %CP(s) =  data_efr.subinfo.CP;
        CP(s) = uheal_data.CP_new(find(uheal_data.subid==sub_num(s)));
        HV(s) = data_efr.subinfo.HV;
        rjt_trials{s} = data_efr.rjt_trials;
        nr_reject(s) =data_efr.nr_reject;
        chan_labels{s} = data_efr.chan_labels;
        chans{s} = data_efr.channels;
        
    elseif ~isfield(data_efr.subinfo,'age')
        subinfo{s} = data_efr.subinfo;
        age(s) = nan;
    else
        
        subinfo{s} = data_efr.subinfo;
        age(s) = nan;
    end
    
    cd(thisdir)
    
end
 
%% NH idx & significant
for s=1:length(F)
    for cc=1:length(chans{s})
        if F{s}(cc)>=F_crit{s}(cc) %
            sig_idx(s,cc) = 1;
        else
            sig_idx(s,cc) = 0;
        end
    end
end
YNH_idx =find((CP==0 & ~isnan(age) & age<23 & sig_idx(:,10)'==1));
ONH_idx = find((CP==0 & ~isnan(age) & age>45) & sig_idx(:,10)'==1);
NH_idx = find(CP==0 & ~isnan(age) & sig_idx(:,10)'==1);

% simple sig (Cz)
sig_cz = find(sig_idx(:,10));
%% extract data (all chans)
close all

%cidx = [17,18]; %[5,11]%
TS_sub = []; fft_sub = [];
for s=1:length(TS)
    if ~isempty(chans{s})
        cidx = find(strcmp(chan_labels{s},'Cz'));
        chanoi = [1:16]%cidx;
        TS_sub(s,:,:) = TS{s}(chanoi,:); % timeseries
        fft_sub(s,:,:) = f_fft{s}(chanoi,:,:); % spectrum
        F_sub(s,:) = F{s}(chanoi); % F-stat
        FFR_sub(s,:) = FFR{s}(chanoi);
        SNR_sub(s,:) = FFR_SNR{s}(chanoi);
        noise_sub(s,:) = fft_noise{s}(chanoi);
        noise_f_sub(s,:,:) = noisef{s}(chanoi,:);
    else
        TS_sub(s,:,:) = nan(1,length(chanoi),2048);
        fft_sub(s,:,:) = nan(1,length(chanoi),1025);
        F_sub(s,:) = nan(1,16);
        FFR_sub(s,:) = nan(1,16);
        SNR_sub(s,:) = nan(1,16);
    end
    
end

%% save to UHEAL_data
load('/work1/jonmarc/UHEAL_master/UHEAL/UHEAL_data/scraped/uheal_data_table/uheal_data.mat')
uheal_data.EFR_SNR = nan(size(uheal_data.subid));
uheal_data.EFR_sig = nan(size(uheal_data.subid));

%%
for s=1:length(SNR_sub)
    % get this subid
    thisID = sub_num(s);
    this_idx = find(uheal_data.subid==thisID);
    
    uheal_data.EFR_SNR(this_idx) = db(F_sub(s,10));
    uheal_data.EFR_sig(this_idx) = sig_idx(s,10)';
end
thisdir = cd;
cd('/work1/jonmarc/UHEAL_master/UHEAL/UHEAL_data/scraped/uheal_data_table/')
uheal_table = struct2table(uheal_data)

writetable(uheal_table,'uheal_data.csv')  
save('uheal_data.mat','uheal_data')

cd(thisdir)
%% plotting F stats
close all
subplot(2,2,[1 2])

%Young
shadedErrorBar(1:16,nanmean(F_sub(YNH_idx,:)),nanstd(F_sub(YNH_idx,:))/sqrt(length(YNH_idx)),'transparent',0);
hold on
py=plot(1:16,nanmean(F_sub(YNH_idx,:)),'b')
% old
shadedErrorBar(1:16,nanmean(F_sub(ONH_idx,:)),nanstd(F_sub(ONH_idx,:))/sqrt(length(ONH_idx)),'transparent',0);
po=plot(1:16,nanmean(F_sub(ONH_idx,:)),'r')
fcrit=plot(1:16,F_crit{1}(1:16),'k--')

set(gca,'xtick',1:16,'xticklabel',chan_labels{1}(1:16))
hleg = legend([py, po,fcrit],'Young','Older','F_{Crit}');
hleg.Box = 'off'
%hleg.Position = [0.4738 0.6135 0.1518 0.1238];
xlabel('Channel')
ylabel('F-statistic')
%ylim([4 16.5])
%ylim([2 18])
subplot(2,2,3)
zlim = [10 35];
%zlim = [];
c=jm_topoplot(nanmean(F_sub(YNH_idx,:))',zlim,'YNH')
c.Label.String = 'F-stat';
%c.Label.Rotation = 0

subplot(2,2,4)
c=jm_topoplot(nanmean(F_sub(ONH_idx,:))',zlim,'ONH')
c.Label.String = 'F-stat';
%c.Label.Rotation = 0


set(gcf,'position',[441 318 560 420])
fig = gcf;
saveas(fig,'figs_thesis/F_yvo','epsc')
%% plotting FFR amp
close all
subplot(2,2,[1 2])

%Young
shadedErrorBar(1:16,nanmean(db(FFR_sub(YNH_idx,:))),nanstd(db(FFR_sub(YNH_idx,:)))/sqrt(length(YNH_idx)),'transparent',0);
hold on
py=plot(1:16,nanmean(db(FFR_sub(YNH_idx,:))),'b')
% old
shadedErrorBar(1:16,nanmean(db(FFR_sub(ONH_idx,:))),nanstd(db(FFR_sub(ONH_idx,:)))/sqrt(length(ONH_idx)),'transparent',0);
po=plot(1:16,nanmean(db(FFR_sub(ONH_idx,:))),'r')
%fcrit=plot(1:16,F_crit{1}(1:16),'k--')

set(gca,'xtick',1:16,'xticklabel',chan_labels{1}(1:16))
hleg = legend([py, po],'Young','Older');
hleg.Box = 'off'
hleg.Position = [0.4738 0.6135 0.1518 0.1238];
xlabel('Channel')
ylabel('EFR amplitude')
%ylim([4 16.5])
%ylim([2 18])
subplot(2,2,3)
zlim = [-90 -70];

%zlim = [];
c=jm_topoplot(nanmean(db(FFR_sub(YNH_idx,:)))',zlim,'YNH');
c.Label.String = 'dB mV';
subplot(2,2,4)
c=jm_topoplot(nanmean(db(FFR_sub(ONH_idx,:)))',zlim,'ONH');
c.Label.String = 'dB mV';


set(gcf,'position',[441 318 560 420])
fig = gcf;
saveas(fig,'figs_thesis/amp_yvo','epsc')
%% plotting SNR
close all
subplot(2,2,[1 2])

%Young
shadedErrorBar(1:16,nanmean(db(F_sub(YNH_idx,:))),nanstd(db(F_sub(YNH_idx,:)))/sqrt(length(YNH_idx)),'transparent',0);
hold on
py=plot(1:16,nanmean(db(F_sub(YNH_idx,:))),'b')
% old
shadedErrorBar(1:16,nanmean(db(F_sub(ONH_idx,:))),nanstd(db(F_sub(ONH_idx,:)))/sqrt(length(ONH_idx)),'transparent',0);
po=plot(1:16,nanmean(db(F_sub(ONH_idx,:))),'r')
%fcrit=plot(1:16,F_crit{1}(1:16),'k--')

set(gca,'xtick',1:16,'xticklabel',chan_labels{1}(1:16))
hleg = legend([py, po],'Young','Older');
hleg.Box = 'off'
hleg.Position = [0.4738 0.6135 0.1518 0.1238];
xlabel('Channel')
ylabel('SNR (dB)')
%ylim([4 16.5])
%ylim([2 18])
subplot(2,2,3)
zlim = [10 30];
%zlim = [];
c=jm_topoplot(nanmean(db(F_sub(YNH_idx,:)))',zlim,'YNH');
c.Label.String = 'SNR (dB)';
subplot(2,2,4)
c=jm_topoplot(nanmean(db(F_sub(ONH_idx,:)))',zlim,'ONH');
c.Label.String = 'SNR (dB)';


set(gcf,'position',[441 318 560 420])
fig = gcf;
saveas(fig,'figs_thesis/snr_yvo','epsc')

%% Bars F stat
YNH_non =find((CP==0 & ~isnan(age) & age<23 & sig_idx(:,10)'==0));
ONH_non = find((CP==0 & ~isnan(age) & age>45 & sig_idx(:,10)'==0));
close all
bar(1:2,[nanmean(F_sub(YNH_idx,10))' nanmean(F_sub(ONH_idx,10))'],'FaceColor',[0.6 0.6 0.6])
hold on
scatter(ones(size(YNH_idx)),F_sub(YNH_idx,10),'sizedata',85,'XJitter','rand','XJitterwidth',0.2,'marker','.','MarkerEdgeColor','k','MarkerFaceColor','k','markerfacealpha',0.5)
scatter(ones(size(ONH_idx))*2,F_sub(ONH_idx,10),'sizedata',85,'XJitter','rand','XJitterwidth',0.2,'marker','.','MarkerEdgeColor','r','MarkerFaceColor','r','markerfacealpha',0.5)
scatter(ones(size(YNH_non)),F_sub(YNH_non,10),'k+')
scatter(ones(size(ONH_non))*2,F_sub(ONH_non,10),'r+')
plot([0 3],[F_crit{1}(1:2)],'k--')

set(gca,'xtick',[1 2],'xticklabel',{'YNH','ONH'},'fontsize',14)
box off
ylabel('F-Statistic (Cz)')
set(gcf,'position',[441 502 270 223])
fig = gcf;
%saveas(fig,'figs2/f_bar','epsc')
%% 3 groups
YNH_idx =find((CP==0 & ~isnan(age) & age<=25 & sig_idx(:,10)'==1));
MNH_idx = find((CP==0 & ~isnan(age) & age>25 & age<60) & sig_idx(:,10)'==1);
ONH_idx = find((CP==0 & ~isnan(age) & age>60) & sig_idx(:,10)'==1);
NH_idx = find(CP==0 & ~isnan(age) & sig_idx(:,10)'==1);
cm = cbrewer('qual','Set1',10)
cmap = cm([1 2 10],:);
close all
figure('renderer','painter')
%subplot(2,3,[1 3])

%Young
shadedErrorBar(1:16,nanmean(F_sub(YNH_idx,:)),nanstd(F_sub(YNH_idx,:))/sqrt(length(YNH_idx)),'lineprops',{'color',cmap(3,:)},'transparent',0);
hold on
py=plot(1:16,nanmean(F_sub(YNH_idx,:)),'color',cmap(3,:))
% old
shadedErrorBar(1:16,nanmean(F_sub(MNH_idx,:)),nanstd(F_sub(MNH_idx,:))/sqrt(length(MNH_idx)),'lineprops',{'color',cmap(2,:)},'transparent',0);
pm=plot(1:16,nanmean(F_sub(MNH_idx,:)),'color',cmap(2,:))
shadedErrorBar(1:16,nanmean(F_sub(ONH_idx,:)),nanstd(F_sub(ONH_idx,:))/sqrt(length(ONH_idx)),'lineprops',{'color',cmap(1,:)},'transparent',0);
po=plot(1:16,nanmean(F_sub(ONH_idx,:)),'color',cmap(1,:))
%fcrit=plot(1:16,F_crit{1}(1:16),'k--')

set(gca,'xtick',1:16,'xticklabel',chan_labels{1}(1:16),'fontsize',14)
%xtickangle(45)
hleg = legend([py, pm, po],'Young','Middle-aged','Older');
hleg.Box = 'off'
hleg.Position = [0.4465 0.2512 0.1940 0.2347];
xlabel('Channel')
ylabel('SNR (dB)')
set(gcf,'position',[258 425 799 277])
fig = gcf;
saveas(fig,'figs_thesis/snr_YMO_chans','epsc')
%ylim([4 16.5])
%ylim([2 18])
figure('renderer','painter')
subplot(3,2,1)
zlim = [5 30];
c=jm_topoplot(nanmean(SNR_sub(YNH_idx,:))',zlim,'Young');
c.Label.String = 'SNR (dB)';
c.FontSize = 10;
subplot(3,2,3)
c=jm_topoplot(nanmean(SNR_sub(MNH_idx,:))',zlim,'Middle-aged');
c.Label.String = 'SNR (dB)';
c.FontSize = 10;
subplot(3,2,5)
c=jm_topoplot(nanmean(SNR_sub(ONH_idx,:))',zlim,'Older');
c.Label.String = 'SNR (dB)';
c.FontSize = 10;

set(gcf,'position',[441 318 560 420])
fig = gcf;
saveas(fig,'figs_thesis/snr_YMO_top','epsc')
%% bar amp
%% Bars F stat
YNH_non =find((CP==0 & ~isnan(age) & age<23 & sig_idx(:,10)'==0));
ONH_non = find((CP==0 & ~isnan(age) & age>45 & sig_idx(:,10)'==0));
close all
bar(1:2,[nanmean(db(FFR_sub(YNH_idx,10)))' nanmean(db(FFR_sub(ONH_idx,10)))'],'FaceColor',[0.6 0.6 0.6])
hold on
scatter(ones(size(YNH_idx)),db(FFR_sub(YNH_idx,10)),'sizedata',85,'XJitter','rand','XJitterwidth',0.2,'marker','.','MarkerEdgeColor','k','MarkerFaceColor','k','markerfacealpha',0.5)
scatter(ones(size(ONH_idx))*2,db(FFR_sub(ONH_idx,10)),'sizedata',85,'XJitter','rand','XJitterwidth',0.2,'marker','.','MarkerEdgeColor','r','MarkerFaceColor','r','markerfacealpha',0.5)
scatter(ones(size(YNH_non)),db(FFR_sub(YNH_non,10)),'k+')
scatter(ones(size(ONH_non))*2,db(FFR_sub(ONH_non,10)),'r+')
%plot([0 3],[F_crit{1}(1:2)],'k--')

set(gca,'xtick',[1 2],'xticklabel',{'YNH','ONH'},'fontsize',14)
box off
ylabel('FFR amplitude (Cz)')
set(gcf,'position',[441 502 270 223])
fig = gcf;
%saveas(fig,'figs2/amp_bar','epsc')

%% time series
%%
close all
figure('renderer','painter');
% filter for visualization
y_ts = squeeze(nanmean(TS_sub(YNH_idx,10,:)));
m_ts = squeeze(nanmean(TS_sub(MNH_idx,10,:)));
o_ts = squeeze(nanmean(TS_sub(ONH_idx,10,:)));
ts_all = squeeze(nanmean(TS_sub(find(sig_idx(:,10)),10,:)));

filt_coef = [100 200];
fs = 4096;
filt_def = designfilt('bandpassfir','FilterOrder',40, ...
    'CutoffFrequency1',filt_coef(1),'CutoffFrequency2',filt_coef(2), ...
    'SampleRate',fs);
y_ts =filtfilt(filt_def,y_ts);
m_ts = filtfilt(filt_def,m_ts);
o_ts = filtfilt(filt_def,o_ts);
ts_all = filtfilt(filt_def,ts_all);
t = time(find(time>=0 & time<=0.5));
plot(t,y_ts','color','k');
hold on
plot(t,m_ts','color',cmap(2,:));
plot(t,o_ts','color',cmap(1,:));
box off
hleg = legend('Young','Middle-aged','Older')
hleg.Box = 'off';
hleg.Position = [0.5518 0.6662 0.4607 0.3260];
xlabel('Time [s]');
ylabel('mV');set(gca,'fontsize',14)
set(gcf,'position',[441 475 369 227])
ylim([-0.15 0.15])
xlim([0 0.5])
title('')
fig = gcf;
saveas(fig,'figs_thesis/ts_ymo','epsc')

% all

figure('renderer','painter')
plot(t,ts_all,'k')
hold on
plot([0 0.3],[-0.08 -0.08],'-','color',[0.5 0.5 0.5],'linewidth',2)
xlabel('Time [s]');
ylabel('mV');set(gca,'fontsize',16)
set(gcf,'position',[441 475 369 227])
xlim([0 0.5])
ylim([-0.1 0.1])
box off
%title(['all sig. n= ' num2str(length(find(sig_idx(:,10))))])
fig = gcf;
saveas(fig,'figs_thesis/ts_all','epsc')
%% plotting spectra
figure
plot(fft_freq(1,:),squeeze(nanmean((db(fft_sub(YNH_idx,10,:))))),'k')
hold on
plot(fft_freq(1,:),squeeze(nanmean(db(fft_sub(ONH_idx,10,:)))),'r')
%ylim([-110 -50])
xlim([100 200])
box off
hleg = legend('YNH','ONH')
hleg.Box = 'off';
xlabel('Frequency [Hz]');
ylabel('db mV');set(gca,'fontsize',16)
set(gcf,'position',[441 498 560 227])
fig = gcf;
saveas(fig,'figs_thesis/freq_yvso','epsc')


%% age vs. FFR
close all
this_idx = find(sig_idx(:,10));
non_idx = find(sig_idx(:,10)==0);
scatter(age(NH_idx),db(FFR_sub(NH_idx,10)),'o','markerfacecolor',[0.6 0.6 0.6],'MarkerEdgecolor','k')
set(gca,'fontsize',12)
ll=lsline
set(ll,'linewidth',2,'color','k')
hold on
scatter(age(non_idx),db(FFR_sub(non_idx,10)),'r+')
ylabel('EFR [dB mV]')

xlabel('Age')
set(gcf,'Position',[228 420 280 209]);
[rho,pval]=corr(age(NH_idx)',db(squeeze(FFR_sub(NH_idx,10))))
fig = gcf;
saveas(fig,'figs_thesis/age_corr_ffr','epsc')

%% F-stat corr

close all
this_idx = find(sig_idx(:,10));
non_idx = find(sig_idx(:,10)==0);
scatter(age(this_idx),db(F_sub(this_idx,10)),'o','markerfacecolor',[0.6 0.6 0.6],'MarkerEdgecolor','k')
set(gca,'fontsize',12)
ll=lsline
set(ll,'linewidth',2,'color','k')
hold on
scatter(age(non_idx),db(F_sub(non_idx,10)),'r+')
ylabel('F-Statistic')

xlabel('Age')
set(gcf,'Position',[228 420 280 209]);
plot([10 80], db(F_crit{1}(:,1:2)),'k--')
[rho,pval]=corr(age(NH_idx)',db(F_sub(NH_idx,10)))
fig = gcf;
saveas(fig,'figs_thesis/age_corr_f','epsc')

%% SNR corr

close all
this_idx = find(sig_idx(:,10));
non_idx = find(sig_idx(:,10)==0);
scatter(age(:),SNR_sub(:,10),'o','markerfacecolor',[0.6 0.6 0.6],'MarkerEdgecolor','k')
set(gca,'fontsize',12)
ll=lsline
set(ll,'linewidth',2,'color','k')
hold on
scatter(age(non_idx),SNR_sub(non_idx,10),'r+')
ylabel('SNR (dB)')

xlabel('Age')
set(gcf,'Position',[228 420 280 209]);
plot([10 80], db(F_crit{1}(:,1:2)),'k--')
[rho,pval]=corr(age(NH_idx)',SNR_sub(NH_idx)')
fig = gcf;
saveas(fig,'figs_thesis/age_corr_snr','epsc')


%% noise floor
close all
scatter(age(NH_idx),db(noise_sub(NH_idx,10)),'o','markerfacecolor',[0.6 0.6 0.6],'MarkerEdgecolor','k')
xlabel('Age')
ylabel('Noise floor (dB mV)')
ll=lsline
set(gca,'fontsize',12)
set(ll,'linewidth',2,'color','k')
set(gcf,'Position',[228 420 280 209]);
[rho,pval]=corr(age(NH_idx)',db(noise_sub(NH_idx,10)))
fig = gcf;
saveas(fig,'figs/noise_floor','epsc')



%% exploratory subtracting noisefloor, then look at amplitudes
% noise estimation
close all
for s=1:size(fft_sub,1)
c_ffr(s,:) = squeeze(fft_sub(s,10,:)-noise_f_sub(s,10,:));
FFR_c(s,:) = c_ffr(s,find(fft_freq(1,:)==326));
end

%plotting spectra

plot(fft_freq(1,:),db(squeeze(nanmean((c_ffr(YNH_idx,:))))),'k')
hold on
plot(fft_freq(1,:),db(squeeze(nanmean(c_ffr(ONH_idx,:)))),'r')
%ylim([-110 -50])
xlim([100 450])
box off
hleg = legend('YNH','ONH')
hleg.Box = 'off';
xlabel('Frequency [Hz]');
ylabel('db mV');set(gca,'fontsize',16)
set(gcf,'position',[441 498 560 227])
fig = gcf;
saveas(fig,'figs/freq_yvso_corrected','epsc')

%% FFR corrected

close all
this_idx = find((sig_idx(:,10)' & ~isnan(age)));
non_idx = find(sig_idx(:,10)==0);
scatter(age(:),db(FFR_c(:)),'o','markerfacecolor',[0.6 0.6 0.6],'MarkerEdgecolor','k')
set(gca,'fontsize',12)
ll=lsline
set(ll,'linewidth',2,'color','k')
hold on
scatter(age(non_idx),db(FFR_c(non_idx)),'r+')
ylabel('FFR corrected-amp')

xlabel('Age')
set(gcf,'Position',[228 420 280 209]);
%plot([10 80], F_crit{1}(:,1:2),'k--')
[rho,pval]=corr(age(this_idx)',db(FFR_c(this_idx)))
fig = gcf;
%saveas(fig,'figs2/age_corr_ffr_corrected','epsc')
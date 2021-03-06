%% plot FFR results
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
%% get data
for s=1:length(subs)
    
    load(subs(s).name)
    clc
    disp(['sub ' subs(s).name(1:5) ' loaded...'])
    % get FFR
    if isfield(data_ffr,'itpc')
        itpc{s} = data_ffr.itpc;
        f_fft{s} = data_ffr.f_fft;
        fft_freq(s,:) = data_ffr.fft_freq;
        TS{s} = data_ffr.TS;
        %trials{s} = data_ffr.trials;
        time = data_ffr.time;
        tidx = data_ffr.tidx;
        stimear(s,:) = data_ffr.stimear;
        subinfo{s} = data_ffr.subinfo;
        if isempty(data_ffr.subinfo.age)|isempty(data_ffr.subinfo.gender)
            age(s) = nan;
            gender(s) = nan;
        else
        age(s) =data_ffr.subinfo.age;
        gender(s) = data_ffr.subinfo.gender;
        end

        if strcmp(data_ffr.subinfo.lab,'PHY2')
            lab(s) = 1;
        elseif strcmp(data_ffr.subinfo.lab,'PHY1')
            lab(s) = 2;
        elseif strcmp(data_ffr.subinfo.lab,'CHBC')
            lab(s) = 3;
        else
            lab(s) = nan;
        end
        CP(s) =  data_ffr.subinfo.CP;
        HV(s) = data_ffr.subinfo.HV;
        rjt_trials{s} = data_ffr.rjt_trials;
        nr_reject(s) =data_ffr.nr_reject;
        chan_labels{s} = data_ffr.chan_labels;
        chans{s} = data_ffr.channels;
        
    elseif ~isfield(data_ffr.subinfo,'age')
        subinfo{s} = data_ffr.subinfo;
        age(s) = nan;
    else
        
        subinfo{s} = data_ffr.subinfo;
        age(s) = nan;
    end
    
    cd(thisdir)
    
end
 

%% extract data (all chans)
close all

%cidx = [17,18]; %[5,11]%
TS_sub = []; fft_sub = [];
for s=1:length(TS)
    if ~isempty(chans{s})
        cidx = find(strcmp(chan_labels{s},'Cz'));
        chanoi = [1:16]%cidx;
        itpc_sub(s,:,:) = itpc{s}(chanoi,:); % timeseries
        fft_sub(s,:,:) = f_fft{s}(chanoi,:,:); % spectrum
        TS_sub(s,:,:) = TS{s}(chanoi,:,:);
    else
        itpc_sub(s,:,:) = nan(1,length(chanoi),1537);
        TS_sub(s,:,:) = nan(1,length(chanoi),3587);
        fft_sub(s,:,:) = nan(1,length(chanoi),1844);

    end
    
end
%%

YNH_idx = find(age<25 & ~CP);
ONH_idx = find(age>45 & ~CP);
chans = setdiff(1:16,[5 11]);
close all
figure('Renderer','painter')
eb1=shadedErrorBar(fft_freq(1,:),squeeze(nanmean(nanmean(itpc_sub(YNH_idx,chans,:)),2)),...
    squeeze(nanstd(nanmean(itpc_sub(YNH_idx,chans,:),2)))/sqrt(length(YNH_idx)))
hold on
eb2=shadedErrorBar(fft_freq(1,:),squeeze(nanmean(nanmean(itpc_sub(ONH_idx,chans,:)),2)),...
   squeeze(nanstd(nanmean(itpc_sub(ONH_idx,chans,:),2)))/sqrt(length(ONH_idx)),'lineprops', 'r')
xlim([1 20])
hleg = legend([eb1.mainLine eb2.mainLine],'YNH','ONH');
hleg.Box = 'off'
box off
ylabel('ITPC')
xlabel('frequency (hz)')
set(gca,'fontsize',12)
set(gcf,'position',[441 426 560 276])

fig = gcf;
saveas(fig,'figs/itpc_spec','epsc')
%%
close all
figure('renderer','painter')
top_idx ={YNH_idx,ONH_idx};
%subplot(1,4,1)
for fid = [2:2:8]
    subplot(3,4,fid/2)
    zlim = [0.1 0.4];
    jm_topoplot(squeeze(nanmean(itpc_sub(YNH_idx,1:16,find(fft_freq(1,:) == fid))))',zlim,[num2str(fid) ' Hz'],0);
    
    subplot(3,4,fid/2+4)
    jm_topoplot(squeeze(nanmean(itpc_sub(ONH_idx,1:16,find(fft_freq(1,:) == fid))))',zlim,[num2str(fid) ' Hz'],0);
end
subplot(3,4,12)
c=jm_topoplot(zeros(16,1),zlim,'scale',1)
c.Label.String = 'itpc';
fig = gcf;
saveas(fig,'figs/itpc_top','epsc')


%% time series
close all
figure('Renderer','painter')
eb1=shadedErrorBar(time,squeeze(nanmean(nanmean(TS_sub(YNH_idx,chans,:)),2)),...
    squeeze(nanstd(nanmean(TS_sub(YNH_idx,chans,:),2)))/sqrt(length(YNH_idx)))
hold on
eb2=shadedErrorBar(time,squeeze(nanmean(nanmean(TS_sub(ONH_idx,chans,:)),2)),...
   squeeze(nanstd(nanmean(TS_sub(ONH_idx,chans,:),2)))/sqrt(length(ONH_idx)),'lineprops', 'r')
%xlim([1 20])
hleg = legend([eb1.mainLine eb2.mainLine],'YNH','ONH');
hleg.Box = 'off'
hleg.Position = [0.8381 0.5020 0.1500 0.2978];
box off
ylabel('mV')
xlabel('time (s)')
set(gca,'fontsize',12)
set(gcf,'position',[441 426 560 276])
hA=gca;
hA.XRuler.MinorTick ='On';
%set(gcf,'position',[441 566 560 136])
set(gcf,'position',[441 587 560 115])
fig = gcf;
saveas(fig,'figs/itpc_ts','epsc')

%% age vs. ITPC frequencies
fid = [2:2:8];
close all
figure('renderer','painter')
for ii=1:length(fid)
    subplot(2,2,ii)
scatter(age,mean(itpc_sub(:,chans,find(fft_freq(1,:)==fid(ii))),2)','o','markerfacecolor',[0.6 0.6 0.6],'MarkerEdgecolor','k')
set(gca,'fontsize',12)
title([num2str(fid(ii)) ' Hz'])
ll=lsline
set(ll,'linewidth',2,'color','k')
ylim([0 .8])
[rho(ii),pval(ii)]=corr(age',mean(itpc_sub(:,chans,find(fft_freq(1,:)==fid(ii))),2))
end
set(gcf,'position',[441 244 421 458])
fig = gcf;
saveas(fig,'figs/age_itpc','epsc')


%% ratios
close all
scatter(age,...
    mean(itpc_sub(:,chans,find(fft_freq(1,:)==fid(end))),2)'./...
    mean(itpc_sub(:,chans,find(fft_freq(1,:)==fid(1))),2)',...
    'o','markerfacecolor',[0.6 0.6 0.6],'MarkerEdgecolor','k');
set(gca,'fontsize',12)
title(['8/2 Hz'])
ll=lsline
set(ll,'linewidth',2,'color','k')
ylabel('ITPC ratio')
xlabel('Age')
%ylim([0 .8])
ratio = mean(itpc_sub(:,chans,find(fft_freq(1,:)==fid(end))),2)'./...
    mean(itpc_sub(:,chans,find(fft_freq(1,:)==fid(1))),2)'
[rh,pva]=corr(age',ratio')
set(gcf,'position',[441 162 214 197])
fig = gcf;
saveas(fig,'figs/age_itpc_ratio','epsc')
%%
function c=jm_topoplot(var1,zlim,tit_string,coff)
load('/work1/jonmarc/UHEAL_master/UHEAL/_EEG/_func/topo_default.mat');
freq.powspctrm = var1;%nanmean(F_sub(YNH_idx,:))';
cfg = [];
cfg.comment = 'no';
cfg.marker = 'on';
cfg.maarkersymbol = '.';
cfg.layout = 'biosemi64.lay';
cfg.channel = freq.cfg.channel;
cfg.parameter = 'powspctrm';
cfg.style = 'straight';
cfg.zlim = zlim;
ft_topoplotER(cfg,freq);
title(tit_string)
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(brewermap(64,'YlOrRd')) % change the colormap
if coff
c=colorbar;
end
end
%% plot FFR results
clear all;close all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
thisdir = cd;
cd('/work1/jonmarc/UHEAL_master/UHEAL')
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
    if isfield(data_ffr,'FFR')
        FFR{s} = data_ffr.FFR;
        FFR_SNR{s} = data_ffr.FFR_SNR;
        F{s} = data_ffr.F;
        F_crit{s} = data_ffr.F_crit;
        subid{s} = data_ffr.subid;
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
    else
        TS_sub(s,:,:) = nan(1,length(chanoi),2048);
        fft_sub(s,:,:) = nan(1,length(chanoi),1025);
        F_sub(s,:) = nan(1,16);
        FFR_sub(s,:) = nan(1,16);
        SNR_sub(s,:) = nan(1,16);
    end
    
end


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
hleg.Position = [0.4738 0.6135 0.1518 0.1238];
xlabel('Channel')
ylabel('F-statistic')
ylim([4 16.5])
%ylim([2 18])
subplot(2,2,3)
zlim = [10 16];
c=jm_topoplot(nanmean(F_sub(YNH_idx,:))',zlim,'YNH')
c.Label.String = 'F-stat';
%c.Label.Rotation = 0

subplot(2,2,4)
c=jm_topoplot(nanmean(F_sub(ONH_idx,:))',zlim,'ONH')
c.Label.String = 'F-stat';
%c.Label.Rotation = 0


set(gcf,'position',[441 318 560 420])
fig = gcf;
saveas(fig,'figs/F_yvo','epsc')
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
ylabel('FFR amplitude')
%ylim([4 16.5])
%ylim([2 18])
subplot(2,2,3)
zlim = [-60 -45];
c=jm_topoplot(nanmean(db(FFR_sub(YNH_idx,:)))',zlim,'YNH');
c.Label.String = 'dB mV';
subplot(2,2,4)
c=jm_topoplot(nanmean(db(FFR_sub(ONH_idx,:)))',zlim,'ONH');
c.Label.String = 'dB mV';


set(gcf,'position',[441 318 560 420])
fig = gcf;
saveas(fig,'figs/amp_yvo','epsc')
%% plotting SNR
close all
subplot(2,2,[1 2])

%Young
shadedErrorBar(1:16,nanmean(SNR_sub(YNH_idx,:)),nanstd(SNR_sub(YNH_idx,:))/sqrt(length(YNH_idx)),'transparent',0);
hold on
py=plot(1:16,nanmean(SNR_sub(YNH_idx,:)),'b')
% old
shadedErrorBar(1:16,nanmean(SNR_sub(ONH_idx,:)),nanstd(SNR_sub(ONH_idx,:))/sqrt(length(ONH_idx)),'transparent',0);
po=plot(1:16,nanmean(SNR_sub(ONH_idx,:)),'r')
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
zlim = [18 24];
c=jm_topoplot(nanmean(SNR_sub(YNH_idx,:))',zlim,'YNH');
c.Label.String = 'SNR (dB)';
subplot(2,2,4)
c=jm_topoplot(nanmean(SNR_sub(ONH_idx,:))',zlim,'ONH');
c.Label.String = 'SNR (dB)';


set(gcf,'position',[441 318 560 420])
fig = gcf;
saveas(fig,'figs/snr_yvo','epsc')

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
saveas(fig,'figs/f_bar','epsc')

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
saveas(fig,'figs/amp_bar','epsc')

%% time series
%%
close all
% filter for visualization
y_ts = squeeze(nanmean(TS_sub(YNH_idx,10,:)));
o_ts = squeeze(nanmean(TS_sub(ONH_idx,10,:)));
ts_all = squeeze(nanmean(TS_sub(find(sig_idx(:,10)),10,:)));

filt_coef = [200 400];
fs = 4096;
filt_def = designfilt('bandpassfir','FilterOrder',40, ...
    'CutoffFrequency1',filt_coef(1),'CutoffFrequency2',filt_coef(2), ...
    'SampleRate',fs);
y_ts =filtfilt(filt_def,y_ts);
o_ts = filtfilt(filt_def,o_ts);
ts_all = filtfilt(filt_def,ts_all);
t = time(find(time>=0));
plot(t,y_ts','k')
hold on
plot(t,o_ts','r')
box off
hleg = legend('YNH','ONH')
hleg.Box = 'off';
xlabel('Time [s]');
ylabel('mV');set(gca,'fontsize',16)
set(gcf,'position',[441 498 560 227])
ylim([-0.15 0.15])
title('')
fig = gcf;
saveas(fig,'figs/ts_yvso','epsc')


% all

figure
plot(t,ts_all,'k')
xlabel('Time [s]');
ylabel('mV');set(gca,'fontsize',16)
set(gcf,'position',[441 498 560 227])

ylim([-0.15 0.15])
box off
title(['all sig. n= ' num2str(length(find(sig_idx(:,10))))])
fig = gcf;
saveas(fig,'figs/ts_all','epsc')
%% plotting spectra
figure
plot(fft_freq(1,:),squeeze(nanmean(nanmean((db(fft_sub(YNH_idx,:,:)))))),'k')
hold on
plot(fft_freq(1,:),squeeze(nanmean(nanmean(db(fft_sub(ONH_idx,:,:))))),'r')
ylim([-110 -50])
xlim([300 400])
box off
hleg = legend('YNH','ONH')
hleg.Box = 'off';
xlabel('Frequency [Hz]');
ylabel('db mV');set(gca,'fontsize',16)
set(gcf,'position',[441 498 560 227])
fig = gcf;
%saveas(fig,'figs/freq_yvso','epsc')


%% age vs. FFR
close all
this_idx = find(sig_idx(:,10));
non_idx = find(sig_idx(:,10)==0);
scatter(age(this_idx),db(FFR_sub(this_idx,10)),'o','markerfacecolor',[0.6 0.6 0.6],'MarkerEdgecolor','k')
set(gca,'fontsize',12)
ll=lsline
set(ll,'linewidth',2,'color','k')
hold on
scatter(age(non_idx),db(FFR_sub(non_idx,10)),'r+')
ylabel('FFR [dB mV]')

xlabel('Age')
set(gcf,'Position',[228 420 280 209]);

fig = gcf;
saveas(fig,'figs/age_corr_ffr','epsc')

%% SNR corr

close all
this_idx = find(sig_idx(:,10));
non_idx = find(sig_idx(:,10)==0);
scatter(age(this_idx),F_sub(this_idx,10),'o','markerfacecolor',[0.6 0.6 0.6],'MarkerEdgecolor','k')
set(gca,'fontsize',12)
ll=lsline
set(ll,'linewidth',2,'color','k')
hold on
scatter(age(non_idx),F_sub(non_idx,10),'r+')
ylabel('F-Statistic')

xlabel('Age')
set(gcf,'Position',[228 420 280 209]);
plot([10 80], F_crit{1}(:,1:2),'k--')

fig = gcf;
saveas(fig,'figs/age_corr_f','epsc')
%% plotting FFR amplitude and SNR
close all
cmap = flip(cbrewer('div','RdYlBu',77-17)); cmap = cmap./max(cmap(:));
jit = randn(size(FFR)).*0.2;
cidx = [1,10]; %[5,11]%
NH_idx = find(CP==0 & ~isnan(age) & sig_idx(:,10)'==1);
for s=1:length(FFR)
    subplot(1,2,1)
    if length(chans{s})>=18 %18 % all chans
        chanoi = cidx;
    FFR_sub(s,:) = FFR{s}(chanoi);
    FFR_SNR_sub(s,:) = FFR_SNR{s}(chanoi);
    %plot(1+jit(s),db(FFR_sub(s,:))','o','markerfacecolor',[cmap(age(s)-17,:)],'color',[cmap(age(s)-17,:)])
    %plot(1+jit(s),db(FFR(s))','o','markerfacecolor',[0.6 0.6 0.6],'color','k')
    xlim([0 2])
    hold on
    xlabel('')
    ylabel('FFR amplitude (dB mV)')
    set(gca,'xtick',[],'fontsize',16,'Xcolor','none')
    box off
 
    subplot(1,2,2)
    %plot(1+jit(s),FFR_SNR_sub(s,:)','o','markerfacecolor',[cmap(age(s)-17,:)],'color',[cmap(age(s)-17,:)])
    %plot(1+jit(s),FFR_SNR(s)','o','markerfacecolor',[0.6 0.6 0.6],'color','k')
    xlim([0 2])
    hold on
    xlabel('')
    ylabel('SNR (dB)')
    set(gca,'xtick',[],'fontsize',16,'Xcolor','none')
    box off
    else 
        FFR_sub(s) = nan;
        FFR_SNR_sub(s) = nan;
        
    end
end
cb=colorbar;
 
cb.FontSize = 12;
cb.Limits = [0 1]
cb.Ticks = [linspace(0,1,5)];
cb.TickLabels = {linspace(18,70,5)};
cb.Label.String = 'Age';
cb.Label.Rotation = 90;
cb.Label.FontSize = 16;
cb.Label.FontName = 'Arial';
cb.Location = 'eastoutside';
colormap(cmap)
 
set(gcf,'position',[305 412 432 299])
 
%NH_idx =(CP==0 & ~isnan(age) & ~isnan(FFR_sub)');
figure
subplot(2,2,1)
scatter(age(NH_idx),db(FFR_sub(NH_idx)))
[rho,pval]=corr(age(NH_idx)',db(FFR_sub(NH_idx))')
lsline 
subplot(2,2,2)
scatter(age(NH_idx),db(FFR_SNR_sub(NH_idx)))
[rho,pval]=corr(age(NH_idx)',FFR_SNR_sub(NH_idx)')
lsline
subplot(2,2,3)
%bar([1 2],[nanmean(db(FFR_sub(find(age<25)))) nanmean(db(FFR_sub(find(age>45))))])
ynh = [];onh=[];
ynh = db(FFR_sub(YNH_idx))';
onh = db(FFR_sub(ONH_idx))';
if length(ynh)<length(onh)
    ynh(length(ynh)+1:length(onh)) = nan;
else
    onh(length(onh)+1:length(ynh)) = nan;
end
bar([ynh;onh]')
%subplot(2,2,4)
%bar([1 2;1 2],[nanmean(FFR_SNR_sub(YNH_idx,:),1);nanmean(FFR_SNR_sub(YNH_idx,:),1)])
%% plot and highlight some
close all
 
jit = randn(size(FFR)).*0.2;
cc = {'r','b'};
out = [47 18];
 
for s=1:length(FFR)
  
    subplot(1,2,1)
    if any(s==out)
        plot(1+jit(s),db(FFR(s))','o','markerfacecolor','r','color','r')
        text(1+jit(s)-.2,db(FFR(s))',num2str(s))
    else
        plot(1+jit(s),db(FFR(s))','o','markerfacecolor','b','color','b')
    end
    xlim([0 2])
    hold on
    xlabel('')
    ylabel('FFR amplitude (dB mV)')
    set(gca,'xtick',[],'fontsize',16,'Xcolor','none')
    box off
 
    subplot(1,2,2)
    if any(s==out)
        plot(1+jit(s),FFR_SNR(s)','o','markerfacecolor','r','color','r')
        text(1+jit(s)-.2,FFR_SNR(s)',num2str(s))
    else
        plot(1+jit(s),FFR_SNR(s)','o','markerfacecolor','b','color','b')
    end
    xlim([0 2])
    hold on
    xlabel('')
    ylabel('SNR (dB)')
    set(gca,'xtick',[],'fontsize',16,'Xcolor','none')
    box off
end
 
%% old vs young
close all
%HI_idx =([14,24,31,42,45,46]);
yidx = find(age_FFR<30);
oidx = find(age_FFR>48);
oidx = setdiff(oidx,HI_idx);
[i,m]=maxk(FFR(oidx),4);
%oidx = setdiff(oidx,oidx(m));
%%
subplot(1,2,1)
%bar([1 2],db([mean(FFR(yidx)) mean(FFR(oidx))]))
hold on
scatter(ones(size(yidx)),db(FFR(yidx)))
hold on
scatter(ones(size(oidx))*2,db(FFR(oidx)))
 
xlim([0 3])
set(gca,'xtick',[1 2],'xticklabel',{'Young','Older'})
ylabel('FFR amplitude (dB)')
 
subplot(1,2,2)
bar([1 2],([mean(FFR(yidx)) mean(FFR(oidx))]))
hold on
scatter(ones(size(yidx)),(FFR(yidx)))
hold on
scatter(ones(size(oidx))*2,(FFR(oidx)))
xlim([0 3])
set(gca,'xtick',[1 2],'xticklabel',{'Young','Older'})
ylabel('FFR amplitude (mV)')
 
a ={FFR(yidx);FFR(oidx)}';
[t df pvals surog] = statcond(a,'paired','on', 'mode', 'perm', 'naccu', 10000);
pvals
 
set(gcf,'position',[680 563 329 242])
%% spectrum
 
close all
%for s=1:length(subs)
for kk=60
    figure(1)
    %subplot(6,7,s)
 
 
%ylim([0 6e-4])
xlim([100 1050])
hold on
semilogx(fft_freq(1,:),mean(db(f_fft)),'k')
xlabel('Frequency (Hz)')
ylabel('FFR amplitude')
set(gca,'fontsize',14,'xtick',[0 326 500 1000])
 
set(gcf,'position',[680   243   521   562])
%hleg = legend('326 Hz')
%hleg.Box = 'Off'
box off
%title(subjects{s})
%text(326,2e-4,['SNR ' num2str(round(FFR_SNR(s),2)) 'dB'],'fontsize',6)
% time domain plot
%plot(data_filt)
 
end
%end
%%
%boxplot
clear g
%close all
CP =abs(CP-1);
figure('Position',[440   577   210   221]);
g=gramm('y',mean(SNR(:,[3]),2),'x',CP+1,'color',CP+1)
g.stat_boxplot();
g.set_color_options('map','brewer_paired');
%g(1,1).set_title('Chirp-stim correlation');
g.set_names('x','','y','FFR SNR (dB)');
g.axe_property('XLim',[0 3]);
g.axe_property('Xtick',[1 2],'YLim',[0 62]);
g.axe_property('Xticklabels',{'Young','Old'},'fontsize',14);
no_legend(g(1,1))
g.update('color',CP+1);
%g.set_color_options('map','brewer_paired');
g.set_color_options('chroma',0,'lightness',75);
g.geom_jitter();
no_legend(g)
g.draw();
CP =abs(CP-1);
 
set(gcf,'position',[680   931   250   167])
 
%% plotting (old)
close all
h=figure(11)
subc_y = cbrewer('seq','Blues',length(subjects)+2);subc_y(1:2,:) = [];
subc_o = cbrewer('seq','Reds',length(subjects)+2);
 
%CP= [1,1,1,1,1,1,0,0,0,0,0,0];
% gramm box plot
clear g
close all
figure('Position',[440   499   680   299]);
% Generate label data
SNR_cat = [];n_label = [];f_label = [];yo_label=[];
for i=[1:4]
    SNR_cat = [SNR_cat SNR(:,i)'];
    n_label = [n_label abs(mod(ones(size(CP))*i,2)-1)];
end
f_label = ones(size(n_label));f_label(end/2+1:end)=2;f_label=flip(f_label);
yo_label = repmat(CP,1,4);
% groupnames
groupnames = ['Y','O']
gnames = [];
for i=1:length(yo_label)
    gnames{i} = groupnames(abs(yo_label(i)-2));
end
gnames = gnames'
 
%create gramm
clear g
%noise
g(1,2)=gramm('x',f_label(find(n_label)),'y',SNR_cat(find(n_label)),'color',gnames(find(n_label)));
    %'subset',strcmp(gnames(find(n_label)),'Y'));
%no noise
g(1,1)=gramm('x',f_label(find(~n_label)),'y',SNR_cat(find(~n_label)),'color',gnames(find(~n_label)));
%g(1,2).geom_jitter()
n_labelx = [find(~n_label);find(n_label)];
for i=1:2
    g(1,i).stat_boxplot();
    g(1,i).axe_property('XLim',[0 3],'YLim',[-20 80]);
    g(1,i).axe_property('Xtick',[1 2]);
    g(1,i).axe_property('Xticklabels',{'326','706'},'fontsize',14);
    g(1,i).set_color_options('map','brewer_paired');
    %g(1,i).draw()
    %no_legend(g(1,1))
    %g.draw
    g(1,i).update('x',f_label(n_labelx(i,:))-.25,'subset',strcmp(gnames(n_labelx(i,:)),'O'));
    g(1,i).set_color_options('chroma',0,'lightness',75);
    g(1,i).geom_point();
    g(1,i).draw()
    g(1,i).update('x',f_label(n_labelx(i,:))+.1,'subset',strcmp(gnames(n_labelx(i,:)),'Y'));
    g(1,i).set_color_options('chroma',0,'lightness',75);
    g(1,i).geom_point();
    g(1,i).set_names('x','FFR frequency','y','FFR SNR (dB)','color','Age group');
 
    %g(1,i).draw()
    
end
%end
%no_legend(g(1,1))
%g.set_names('x','cf','y','EFR(4Hz) SNR','color','Age group');
g(1,1).set_title('no noise')
g(1,2).set_title('noise')
g.draw()
 
 
%% old plot
close all
figure(10)
subcolors = cbrewer('seq','Blues',length(subjects)+2);
subred = cbrewer('seq','Reds',length(subjects));
subcolors(end-1,:) = subred(end-1,:);
subcolors(end,:) = subred(end,:)
subcolors(1:2,:) =[];
clear p
subjects={'Y1','Y2','Y3','Y4','Y5','O1','O2','O3'}
for s=1:length(subjects)
subplot(1,2,1)
p(s)=plot(1:2,SNR(s,[1 2]),'-o','color',subcolors(s,:))
hold on
plot(1:2,SNR(s,[3 4]),'-x','color',subcolors(s,:))
xlim([0 3.5])
%ylim([20 65])
set(gca,'xtick',[1,2],'xticklabels',{'Inf','5 dB'},'fontsize',14)
hleg = legend('706 Hz','326 Hz','location','best')
hleg.Box = 'off'
ylabel('dB SNR')
xlabel('Noise condition [SNR]')
grid on
box off
subplot(1,2,2)
plot(1:2,FFR(s,[1 2]),'-o','color',subcolors(s,:))
hold on
plot(1:2,FFR(s,[3 4]),'-x','color',subcolors(s,:))
xlim([0 3])
set(gcf,'position',[440   499   579   299])
ylabel(['FFR amplitude' newline '[dB rel. 1\muV]'])
xlabel('Noise condition [SNR]')
grid on
set(gca,'xtick',[1,2],'xticklabels',{'Inf','5 dB'},'fontsize',14)
ylim([-120 -60])
legend(p(:),subjects)
end
box off
 
%% mean fig
 
close all
figure(10)
subcolors = cbrewer('qual','Set1',length(subjects));
clear p
norm_subs = [1,2,4,5,6];
 
subplot(1,2,1)
p(s)=plot(1:2,mean(SNR(norm_subs,[1 2])),'-o','color',subcolors(1,:))
hold on
plot(1:2,mean(SNR(norm_subs,[3 4])),'-x','color',subcolors(1,:))
xlim([0 3])
%ylim([20 65])
set(gca,'xtick',[1,2],'xticklabels',{'Inf','5 dB'},'fontsize',14)
hleg = legend('706 Hz','326 Hz','location','best')
hleg.Box = 'off'
ylabel('dB SNR')
xlabel('Noise condition [SNR]')
grid on
subplot(1,2,2)
plot(1:2,mean(FFR(norm_subs,[1 2])),'-o','color',subcolors(1,:))
hold on
plot(1:2,mean(FFR(norm_subs,[3 4])),'-x','color',subcolors(1,:))
xlim([0 3])
set(gcf,'position',[440   499   579   299])
ylabel(['FFR amplitude' newline '[dB rel. 1\muV]'])
xlabel('Noise condition [SNR]')
grid on
set(gca,'xtick',[1,2],'xticklabels',{'Inf','5 dB'},'fontsize',14)
ylim([-120 -60])
%legend(p(:),subjects(1:5))

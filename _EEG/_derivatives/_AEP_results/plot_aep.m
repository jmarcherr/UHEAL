%% plot AEP
clear all;close all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
thisdir = cd;
cd '/work1/jonmarc/UHEAL_master/UHEAL'

UHEAL_startup
cd(thisdir)
%subjects
subs = dir('UH*');
 
%%
for s=1:length(subs)

        load(subs(s).name,'aep_avg','n1','timelock_time','tidx','subinfo','short_aep','long_aep','aep')%'timelock')
        
        aep_sub(s,:,:) = aep_avg;
        n1_sub(s,:) = n1;
        short_aep_sub(s,:) = short_aep;
        long_aep_sub(s,:) = long_aep;
        time = timelock_time;
        timeidx = tidx;
        age(s) = subinfo.age;
        gender(s) = subinfo.gender;
        CP_sub(s) = subinfo.CP;
        data_sub(s,:,:) = aep;
    clc
    disp(['sub ' subs(s).name(1:5) ' loaded...'])
end
%%

%%
rates = [0.5,1,1.5,2];
cmap = cbrewer('qual','Set1',4);
for s=1:length(subs)
%     figure(1)
%     subplot(6,7,s)
%     for kk=1:4
%         thisaep  = squeeze(aep_avg(s,kk,:));
%         p1(kk) = plot(time,thisaep,'color',cmap(kk,:))
%         hold on
%     end
%  
%     hold on
%     set(gcf,'position',[107    87   700   994])
%     
%     %
%     hold on
%     %ylabel('AEP Amplitude [\mu V]');
%     %xlabel('Time [s]')
%     box off
%     set(gca,'Fontsize',12)
%
%
%
%     xlim([-.01 .3])
%     ylim([-8 8])


%%
%     figure(2)
%     subplot(6,7,s)

%     p2(1) = plot(time,short_aep)
%     hold on
%     p2(2) = plot(time,long_aep)
%     %ylabel('AEP Amplitude [\mu V]');
%     %xlabel('Time [s]')
%     box off
%     set(gca,'Fontsize',12)
%
%
%     clear colormap
%     %title('AEP short vs. long')
%     cm=cbrewer('qual','Set1',5);
%     colormap(cm)
%     hold on
%
%
%
%     xlim([-.01 .3])
%     ylim([-8 8])
%     set(gcf,'position',[107    87   700   994])

end
% legend(p1(:),num2str(rates'))
% hleg = legend([p2(1) p2(2)],'Short ISI','Long ISI','location','best')
% hleg.Box = 'off'
%


%% mean
close all
figure('Renderer','painter')
cm = lines;
rates = [0.5 1 1.5 2]
YNH_idx = find(CP_sub==0 & age<25);
ONH_idx = find(CP_sub==0 & age>45);
%pmean = plot(time,squeeze(mean(sub_aep(4:19,:,:),1)))

for kk=1:4
    subplot(1,2,1)
    pmean(kk) = shadedErrorBar(time,squeeze(mean(aep_sub(YNH_idx,kk,:),1)),squeeze(std(aep_sub(YNH_idx,kk,:),1))/sqrt(length(find(YNH_idx))),'lineprops',{'-','color',cm(kk,:)},'transparent',1)
    hold on
    plot([time(timeidx(1)) time(timeidx(1))],[-6 5],'k--')
    plot([time(timeidx(end)) time(timeidx(end))],[-6 5],'k--')
    %plot(tn,zeros(size(tn)),'k--');
    %xlabel('Time(s)')
    ylabel('AEP Amplitude (\muV)');
    xlabel('Time [s]')
    box off
    set(gca,'Fontsize',16)
    xlim([-.01 .3])
    ylim([-6 4])
    %hleg = legend([pmean(:).mainLine],num2str(rates(:)),'location','best');
    %hleg.Box = 'off';
    title('YNH')
    
    subplot(1,2,2)
    pmean(kk) = shadedErrorBar(time,squeeze(mean(aep_sub(ONH_idx,kk,:),1)),squeeze(std(aep_sub(ONH_idx,kk,:),1))/sqrt(length(find(ONH_idx))),'lineprops',{'-','color',cm(kk,:)},'transparent',1)
    hold on
    plot([time(timeidx(1)) time(timeidx(1))],[-6 5],'k--')
    plot([time(timeidx(end)) time(timeidx(end))],[-6 5],'k--')
    %plot(tn,zeros(size(tn)),'k--');
    %xlabel('Time(s)')
    ylabel('AEP Amplitude (\muV)');
    xlabel('Time [s]')
    box off
    set(gca,'Fontsize',16)
    xlim([-.01 .3])
    ylim([-6 4])
    title('ONH')

end
    hleg = legend([pmean(:).mainLine],num2str(rates(:)),'location','best');
    hleg.Box = 'off';

%title('AEP Grand Average')
%set(gcf,'position',[371 726 363 304])
xlim([-.01 .3])
ylim([-6 4])
hleg = legend([pmean(:).mainLine],num2str(rates(:)),'location','best');
hleg.Box = 'off';

fig = gcf;
saveas(fig,'figs/yvo_aep','epsc')


%% all
close all
figure('Renderer','painter')
for kk=1:4
    %subplot(1,2,1)
    pmean(kk) = shadedErrorBar(time,squeeze(mean(aep_sub(:,kk,:),1)),squeeze(std(aep_sub(:,kk,:),1))/sqrt(size(aep_sub,1)),'lineprops',{'-','color',cm(kk,:)},'transparent',1)
    hold on
    plot([time(timeidx(1)) time(timeidx(1))],[-6 5],'k--')
    plot([time(timeidx(end)) time(timeidx(end))],[-6 5],'k--')
    %plot(tn,zeros(size(tn)),'k--');
    %xlabel('Time(s)')
    ylabel('AEP Amplitude (\muV)');
    xlabel('Time [s]')
    box off
    set(gca,'Fontsize',16)
    xlim([-.01 .3])
    ylim([-6 4])

    


end
    hleg = legend([pmean(:).mainLine],num2str(rates(:)),'location','best');
    hleg.Box = 'off';
    title('all')
    set(gcf,'position',[441 283 313 420])
 fig = gcf;
saveas(fig,'figs/all_aep','epsc')   
%% short and long
%%
figure('Renderer','painter')
pse(1) = shadedErrorBar(time,mean(short_aep_sub(:,:)),std(short_aep_sub(:,:),1)/sqrt(size(short_aep_sub,1)),'lineprops','-b','transparent',1)
hold on
pse(2) = shadedErrorBar(time,mean(long_aep_sub(:,:)),std(long_aep_sub(:,:),1)/sqrt(size(long_aep_sub,1)),'lineprops','-r','transparent',1)
ylabel('AEP Amplitude (\muV)');
xlabel('Time [s]')
box off
set(gca,'Fontsize',16)
plot([time(timeidx(1)) time(timeidx(1))],[-5 5],'k--')
plot([time(timeidx(end)) time(timeidx(end))],[-5 5],'k--')
 
clear colormap
title('AEP short vs. long')
cm=cbrewer('qual','Set1',5);
%cm(1:10,:)=repmat([1 1 1],10,1);
colormap(cm)
hold on
hleg = legend([pse(1).mainLine pse(2).mainLine],'Short ISI','Long ISI','location','best')
hleg.Box = 'off'
xlim([-.01 .3])
ylim([-4 4])
    title('all')
    set(gcf,'position',[441 283 313 420])
 saveas(fig,'figs/all_aep_SL','epsc')  
%set(gcf,'position',[371 726 363 304])
% hold off

%% group short long
close all
figure('Renderer','painter')
subplot(1,2,1)
pse(1) = shadedErrorBar(time,mean(short_aep_sub(YNH_idx,:)),std(short_aep_sub(YNH_idx,:),1)/sqrt(size(short_aep_sub(YNH_idx,:),1)),'lineprops','-b','transparent',1)
hold on
pse(2) = shadedErrorBar(time,mean(long_aep_sub(YNH_idx,:)),std(long_aep_sub(YNH_idx,:),1)/sqrt(size(long_aep_sub(YNH_idx,:),1)),'lineprops','-r','transparent',1)
ylabel('AEP Amplitude (\muV)');
xlabel('Time [s]')
box off
set(gca,'Fontsize',16)
plot([time(timeidx(1)) time(timeidx(1))],[-6 5],'k--')
plot([time(timeidx(end)) time(timeidx(end))],[-6 5],'k--')
xlim([-.01 .3])
ylim([-6 4])
title('YNH')
subplot(1,2,2)
pse(1) = shadedErrorBar(time,mean(short_aep_sub(ONH_idx,:)),std(short_aep_sub(ONH_idx,:),1)/sqrt(size(short_aep_sub(ONH_idx,:),1)),'lineprops','-b','transparent',1)
hold on
pse(2) = shadedErrorBar(time,mean(long_aep_sub(ONH_idx,:)),std(long_aep_sub(ONH_idx,:),1)/sqrt(size(long_aep_sub(ONH_idx,:),1)),'lineprops','-r','transparent',1)
ylabel('AEP Amplitude (\muV)');
xlabel('Time [s]')
box off
set(gca,'Fontsize',16)
plot([time(timeidx(1)) time(timeidx(1))],[-6 5],'k--')
plot([time(timeidx(end)) time(timeidx(end))],[-6 5],'k--')

clear colormap
title('AEP short vs. long')
cm=cbrewer('qual','Set1',5);
%cm(1:10,:)=repmat([1 1 1],10,1);
colormap(cm)
hold on
hleg = legend([pse(1).mainLine pse(2).mainLine],'Short ISI','Long ISI','location','best');
hleg.Box = 'off'
xlim([-.01 .3])
ylim([-6 4])
    title('ONH')
    %set(gcf,'position',[441 283 313 420])
    fig=gcf
 saveas(fig,'figs/YvO_aep_SL','epsc')  
%set(gcf,'position',[371 726 363 304])
% hold off


%% SLOPE
% figure(100)
for s=1:length(n1_sub)
    %subplot(6,7,s)
    figure(100)
    plot(rates,n1_sub(s,:),'o')
    h=lsline
    %xlabel('ISI (s)')
    %ylabel('n100 amplitude (\mu V)')
    set(gca,'fontsize',10,'xtick',[rates],'xticklabel',round(rates,1))
    box off
    xlim([rates(1) rates(end)])
    slope_sub(s,:) = (h.YData(2)-h.YData(1))/(h.XData(2)-h.XData(1));
    %text(rates(3),-3,['slope = ' num2str(slope)])
    %title(subjects{s})
    slopeX(s,:) = h.XData;
    slopeY(s,:) = h.YData;
end
% set(gcf,'position',[107    87   700   994])
     fig=gcf
 saveas(fig,'figs/Slope','epsc')  
 %%
% mean slope
figure('Renderer','painter')
% get individual slopesc
for s=1:size(aep_sub,1)%length(subs)
plot(slopeX(s,:),slopeY(s,:),'-','color',[.8 .8 .8])
hold on
end
 
hold on
errorbar(rates,mean(n1_sub),nanstd(n1_sub)/sqrt(size(~isnan(n1_sub),1)),'k.')
plot(rates,mean(n1_sub),'ko','markerfacecolor','k')
hmean=lsline
xlabel('ISI (s)')
ylabel('n100 amplitude (\muV)')
set(gca,'fontsize',14,'xtick',[rates],'xticklabel',{'.5','1','1.5','2'})
box off
xlim([0.25 2.25])
slope = (hmean.YData(2)-hmean.YData(1))/(hmean.XData(2)-hmean.XData(1));
%ylim([-6 2])
%text(rates(3),-3,['slope = ' num2str(slope)]),
title(['N1 Slope' ])
%set(gcf,'position',[371 726 363 304])
set(gca,'fontsize',12)
set(gcf,'position',[441 459 270 266])

%%
figure('Renderer','painter')
scatter(age,slope_sub,'o','markerfacecolor',[0.6 0.6 0.6],'MarkerEdgecolor','k')
ll=lsline
title(['N1 Slope' ])
%set(gcf,'position',[371 726 363 304])
set(gca,'fontsize',12)
set(ll,'linewidth',2,'color','k')
set(gcf,'position',[441 459 270 266])
     fig=gcf
 saveas(fig,'figs/Slope_age','epsc') 
 
 %% each isi n1 amp
 close all
 figure('renderer','painter')
 for ii=1:4
     subplot(1,4,ii)
 scatter(age,n1_sub(:,ii),'o','markerfacecolor',[0.6 0.6 0.6],'MarkerEdgecolor','k')
 ll{ii}=lsline
 hold on
 title(rates(ii))
 ylim([-20 4])
 xlim([18 80])
 set(gca,'fontsize',12)
set(ll{ii},'linewidth',2,'color','k')
 end

 set(gcf,'position',[441 506 917 196])
      fig=gcf
 saveas(fig,'figs/n100_age','epsc') 
 
 %% slope in two age groups
 close all
figure('Renderer','painter')
eb(1)=errorbar(rates,nanmean(n1_sub(ONH_idx,:)),nanstd(n1_sub(ONH_idx,:))/sqrt(length(ONH_idx)),'ko','color',[1 0.6 0.6],'markerfacecolor',[1 0.6 0.6],'MarkerEdgecolor','k')
hold on
plot(rates,nanmean(n1_sub(ONH_idx,:)))
ll=lsline
hold on
eb(2)=errorbar(rates,nanmean(n1_sub(YNH_idx,:)),nanstd(n1_sub(YNH_idx,:))/sqrt(length(YNH_idx)),'ko','markerfacecolor',[0.6 0.6 1],'MarkerEdgecolor','k')
plot(rates,nanmean(n1_sub(YNH_idx,:)))
%title(['N1 amplitude' ])
%set(gcf,'position',[371 726 363 304])
ylabel('N100 amplitude (\muV)')
xlabel('isi (s)')
set(gca,'fontsize',12)
set(ll,'linewidth',2,'color','k')
set(gcf,'position',[441 459 270 266])
hleg = legend([eb(2) eb(1)],'YNH','ONH')
xlim([0.25 2.25])
hleg.Box = 'off'
     fig=gcf
     box off
% saveas(fig,'figs/Slope_age','epsc') 
      fig=gcf
 saveas(fig,'figs/n100_age_group','epsc') 
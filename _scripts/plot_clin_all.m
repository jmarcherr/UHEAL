%% plot audiogram

close all
clear all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
cd ..
%addpath('O:\Public\Hearing-Systems-group\cahr\Temporary_ftp\UHEAL')
UHEAL_startup
cd(datadir)
cd([datadir '/scraped']);
freq_aud = [250 500 1000 2000 4000 8000 9000 10000 11200 12500 14000 16000];
d=dir('UH*.mat')

%%
% ages=[];c=0;
% for s=1:length(d)
%     load([d(s).name]);
%     if ~isempty(dataalm.aud)
%         c=c+1;
%         %get audiogram and age
%         [ages(c),aud_L,aud_R,aud_freq,gender] = get_aud(dataalm);
%     end
% end


%% Audiograms
ages = [18 77];
colormap
cmap = flip(cbrewer('div','RdYlBu',max(ages)-17));
cmap(find(cmap>1))=1;
ageidx = linspace(min(ages),max(ages),size(cmap,1));

%%
%figure(1)
for s=1:length(d)
    
    %load relevant file
    %load(['UH' num2str(subid(s))])
    
    load([d(s).name]);
            sub_id{s} = d(s).name;
    if ~isempty(dataalm.aud)
        %get audiogram and age
        [age,aud_L,aud_R,aud_freq,gender] = get_aud(dataalm);
        age = dataalm.subinfo.age;
        gender = dataalm.subinfo.gender;
        CP = dataalm.subinfo.CP;

        %get stim ear
        if strcmp(dataalm.id,'UH099')
            stimear = 1;
        else
            stimear = dataalm.stim.ffr.ear(1);
        end
        %plot stim ear
        if stimear ==1 %left ear
            
            p1 = semilogx(aud_freq,aud_L','-','color',[cmap(age-17,:) 0.8]);
            plot_aud_param(p1,aud_freq);
            hold on
            aud(s,:) = aud_L;
        else % right ear
            p1 = semilogx(aud_freq,aud_R','-','color',[cmap(age-17,:) 0.8]);
            plot_aud_param(p1,aud_freq);
            hold on
            aud(s,:) = aud_R;
        end
        
        age_sub(s) = age; %log age
        gender_sub(s) = gender;
        CP_sub(s) = CP;
    else
        age_sub(s) = nan;
        gender_sub(s) = nan;
        CP_sub(s) = nan;
    end
    clc
    disp(['Subject ' num2str(s) ' processed.'])
end
cb=colorbar;

cb.FontSize = 12;
cb.Limits = [0 1];
cb.Ticks = [linspace(0,1,5)];
cb.TickLabels = {linspace(18,70,5)};
cb.Label.String = 'Age';
cb.Label.Rotation = 90;
cb.Label.FontSize = 16;
cb.Label.FontName = 'Arial';
colormap(cmap)

set(gcf,'position',[305 412 432 299]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% group audiogram plots
close all
idx = find(~isnan(age_sub));
fig = plot_audiogram_groups(idx,age_sub,aud,aud_freq,cmap)
% save figure
title(['all, n=' num2str(length(idx))])
saveas(fig,'figs/audiogram_all','epsc')
%% only NH
close all
NH_idx = find(CP_sub==0);
fig = plot_audiogram_groups(NH_idx,age_sub,aud,aud_freq,cmap)
% save figure
title(['NH, n=' num2str(length(NH_idx))])
saveas(fig,'figs/audiograms_NH','epsc')

%% only HI
close all
HI_idx = find(CP_sub==1);
fig = plot_audiogram_groups(HI_idx,age_sub,aud,aud_freq,cmap)
% save figure
title(['HI, n=' num2str(length(HI_idx))])
saveas(fig,'figs/audiograms_HI','epsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% age plots
% all
close all
hist(age_sub,max(age_sub)-min(age_sub))
xlabel('Age (years)')
ylabel('nr. of participants')
hold on
plot([25 25],[0 8],'k--')
plot([60 60],[0 8],'k--')
set(gca,'fontsize',16)
xlim([16 70])
set(gcf,'position',[305 412 432 299])
%title(['all, n=' num2str(length(idx))])
fig= gcf
saveas(fig,'figs/age_hist_all','epsc')
%% grouped male female
figure(3)
Y   = gender_sub(find(age_sub<=25)); 
O1  = gender_sub(find(age_sub>25 & age_sub<=60));
O2  = gender_sub(find(age_sub>60));

b1=bar([1:3],[length(find(Y==1)) length(find(O1==1)) length(find(O2==1));...
    length(find(Y==2)) length(find(O1==2)) length(find(O2==2))]','stacked')

hold on
xlabel('Group')
xtickangle(45);
ylabel('Nr. of participants')
set(gca,'xtick',[1:3],'xticklabels',{'Y','O1','O2'},'fontsize',16)
hleg = legend(b1,'Female','Male')
b1(1).FaceColor = [0.2 0.2 0.5];
b1(2).FaceColor = [0.5 0.7 0.5];
hleg.Position = [0.6188 0.7419 0.2824 0.1689];
set(gcf,'position',[305 412 432 299])
fig= gcf
saveas(fig,'figs/age_hist_gender','epsc')

%% grouped male female and HI
HI_idx = (CP_sub==1)

figure(4)
Y   = gender_sub(find(age_sub<=25)); 
O1  = gender_sub(find(age_sub(~HI_idx)>25 & age_sub(~HI_idx)<=60));
O1HI = gender_sub(find(age_sub(HI_idx)>25 & age_sub(HI_idx)<=60));
O2  = gender_sub(find(age_sub(~HI_idx)>60));
O2HI = gender_sub(find(age_sub(HI_idx)>60 ));

b1=bar([1:5],[length(find(Y==1)) length(find(O1==1)) length(find(O2==1)) length(find(O1HI==1)) length(find(O2HI==1));...
    length(find(Y==2)) length(find(O1==2)) length(find(O2==2)) length(find(O1HI==2)) length(find(O2HI==2))]','stacked')

hold on
xlabel('Group')
ylabel('Nr. of participants')
set(gca,'xtick',[1:5],'xticklabels',{'Y','O1','O2','O1HI','O2HI'},'fontsize',16)
legend(b1,'Female','Male')
xtickangle(45)
b1(1).FaceColor = [0.2 0.2 0.5]
b1(2).FaceColor = [0.5 0.7 0.5]
%b1(3).FaceColor = [0.2 0.2 0.6];
%b1(4).Facecolor = [0.5 0.7 0.4];
set(gcf,'position',[305 412 432 299])
ylim([0 max(b1(1).YData)+max(b1(2).YData)+5]);
ymax = max(b1(1).YData)+max(b1(2).YData)+5;
fig= gcf
saveas(fig,'figs/age_hist_gender_HI','epsc')

%% grouped male female and HI from hvidovre
HI_idx = zeros(length(gender_sub),1);
HI_idx([14,24,31,42,45,46]) = 1; HI_idx=logical(HI_idx); % from who-what-when

% find who has been to hvidovre (from who-what-when 10-05-21)
DRCMR_idx = zeros(length(gender_sub),1);
DRCMR_idx([26,33,40,12,37,23,35,16,34,6,39,19,8,32,1,36,48,2,41,15,13,...
    11,14,25,24,10,20,42,46,53])=1; DRCMR_idx = logical(DRCMR_idx);

%
O_idx = intersect(find(~HI_idx),find(DRCMR_idx));
OHI_idx =intersect(find(HI_idx),find(DRCMR_idx));


figure(5)
Y   = gender_sub(find(age_sub(DRCMR_idx)<=25)); 
O1  = gender_sub(find(age_sub(O_idx)>25 & age_sub(O_idx)<=60));
O1HI = gender_sub(find(age_sub(OHI_idx)>25 & age_sub(OHI_idx)<=60));
O2  = gender_sub(find(age_sub(O_idx)>60));
O2HI = gender_sub(find(age_sub(OHI_idx)>60 ));




b1=bar([1:5],[length(find(Y==1)) length(find(O1==1)) length(find(O2==1)) length(find(O1HI==1)) length(find(O2HI==1));...
    length(find(Y==2)) length(find(O1==2)) length(find(O2==2)) length(find(O1HI==2)) length(find(O2HI==2))]','stacked')

hold on
xlabel('Group')
ylabel('Nr. of participants')
set(gca,'xtick',[1:5],'xticklabels',{'Y','O1','O2','O1HI','O2HI'},'fontsize',16)
legend(b1,'Female','Male')
xtickangle(45)
b1(1).FaceColor = [0.2 0.2 0.5]
b1(2).FaceColor = [0.5 0.7 0.5]
%b1(3).FaceColor = [0.2 0.2 0.6];
%b1(4).Facecolor = [0.5 0.7 0.4];
set(gcf,'position',[305 412 432 299])
ylim([0 ymax]);
%% MEMR


figure(2)
cmap = cbrewer('seq','YlGnBu',7+2);
cmap = cmap(2:end,:);
plotting = 0;
for s=98:length(d)
    
    % Load scaped data
    load([d(s).name])
    
    if ~isempty(dataalm.stim.ffr)
        stimear = dataalm.stim.ffr.ear(1);
    else
        stimear = 1;
    end
    % does MEMR excist?
    if ~isempty(dataalm.memr)
        % MEMR
        if stimear == 1 % left stim ear
            % does memr exist for stim ear?
            if isfield(dataalm.memr,'L')
                reflex = dataalm.memr.L.reflex_ipsi.response;
                levels = dataalm.memr.L.reflex_ipsi.labels;
                f_center = dataalm.memr.L.reflex_ipsi.f_center;
                freq = dataalm.memr.L.reflex_ipsi.freq;
            else
                reflex = dataalm.memr.R.reflex_ipsi.response;
                levels = dataalm.memr.R.reflex_ipsi.labels;
                f_center = dataalm.memr.R.reflex_ipsi.f_center;
                freq = dataalm.memr.R.reflex_ipsi.freq;
                
            end
        else            % right stim ear
            % does memr exist for stim ear?
            if isfield(dataalm.memr,'R')
                reflex = dataalm.memr.R.reflex_ipsi.response;
                levels = dataalm.memr.R.reflex_ipsi.labels;
                f_center = dataalm.memr.R.reflex_ipsi.f_center;
                freq = dataalm.memr.R.reflex_ipsi.freq;
            else
                reflex = dataalm.memr.L.reflex_ipsi.response;
                levels = dataalm.memr.L.reflex_ipsi.labels;
                f_center = dataalm.memr.L.reflex_ipsi.f_center;
                freq = dataalm.memr.L.reflex_ipsi.freq;
            end
        end
        if plotting
            subplot(10,12,s)
            for ll =1:length(levels)
                semilogx(f_center,reflex(:,ll),'color',cmap(ll,:))
                axis([freq(1) freq(end) -0.1 0.1])
                a=gca;
                a.XTick = [250,500,1000,2000,4000,8000];
                a.XTickLabel = [{'.25'},{'.5'},{'1'},{'2'},{'4'},{'8'}];
                %title(['\Delta absorbance = contracted - baseline. Target pressure: ', num2str(target_p), ' daPa'])
                xlabel('Frequency [kHz]')
                ylabel('\Delta Absorbance')
                grid on
                hold on
                if ~isempty(dataalm.per)
                    title(dataalm.per.PersonNumber)
                else
                    title('UH?')
                end
            end
        end
    else
        warning('No MEMR data for subject')
        if plotting
            subplot(10,12,s)
            a=gca;
            a.XTick = [250,500,1000,2000,4000,8000];
            a.XTickLabel = [{'.25'},{'.5'},{'1'},{'2'},{'4'},{'8'}];
            %title(['\Delta absorbance = contracted - baseline. Target pressure: ', num2str(target_p), ' daPa'])
            xlabel('Frequency [kHz]')
            ylabel('\Delta Absorbance')
            grid on
            hold on
            %title(dataalm.per.PersonNumber)
        end
        reflex = nan(16,7);
    end
    
    f_lim = [6:15];
    [m_fc,I(s)] = max((reflex(f_lim,end)'));
    growth_sub(s,:) = reflex(I(s)+f_lim(1)-1,:);
    growth_sub_alt(s,:) = sum(abs(reflex));
    reflex_sub(s,:,:) = reflex;
    
    MEM_slope(s,:) = polyfit(levels,growth_sub_alt(s,:),1);
end
for x = 1:length(levels)
    txt(x) = {[num2str(levels(x)),' dB SPL']};
end

hleg = legend(txt)
hleg.Box = 'off'
%% average plots
%colormap
cmap = cbrewer('seq','YlGnBu',7+2);
cmap = cmap(2:end,:);

% mean
MEMR_mean = squeeze(nanmean(reflex_sub,1));
% growth curve
figure(3)

p_y = plot(levels,nanmean(growth_sub_alt,1),'-k^','markeredgecolor','k','markerfacecolor','w')
hold on
for ll=1:length(levels)
    errorbar(levels(ll),nanmean(growth_sub_alt(:,ll),1),nanstd(growth_sub_alt(:,ll),1)/sqrt(length(d)),'color',cmap(ll,:))
    eby = plot(levels(ll),nanmean(growth_sub_alt(:,ll),1),'^','color',cmap(ll,:),'MarkerFaceColor',cmap(ll,:))
end

hold on
set(gca,'fontsize',16)
box off
set(gcf,'position',[228 839 432/2 209])
set(gca,'Xtick',levels,'Xticklabels',{'','80','','90','','100',''},'ytick',[0,.2,.4,.6],'yticklabels',{'0','.2','.4','.6'})
xlim([72 108])
% mean results
%%
figure(4)
subplot(1,2,1)
for i=1:length(levels)
    % Plot
    semilogx(f_center, nanmean(reflex_sub(:,:,i),1),'color',cmap(i,:))
    axis([freq(1) freq(end) -0.08 0.08])
    a=gca;
    a.XTick = [250,500,1000,2000,4000,8000];
    a.XTickLabel = [{'.25'},{'.5'},{'1'},{'2'},{'4'},{'8'}];
    hold on
    set(gca,'fontsize',16)
    box off
    for x = 1:length(levels)
        txt(x) = {[num2str(levels(x)),' dB SPL']};
    end
    set(gcf,'position',[228 839 432 209])
    hleg = legend(txt)
    hleg.Box = 'off'
    hleg.Position = [0.5547 0.1231 0.3704 0.8038]
end

%% MEMR slope

close all

figure(4)
MEM_slope(42,:) = nan;
jit = randn(size(MEM_slope,1)).*0.2;
for s=1:length(MEM_slope)
subplot(1,3,1)
%plot(1+jit(s),db(FFR(s))','o','markerfacecolor',[cmap(age(s)-17,:)],'color',[cmap(age(s)-17,:)])
plot(1+jit(s),MEM_slope(s,1)','o','markerfacecolor',[0.6 0.6 0.6],'color','k')
xlim([0 2])
hold on
xlabel('')
ylabel('Linear coefficient')
set(gca,'xtick',[],'fontsize',16,'Xcolor','none')
box off
ylim([-0.02 0.05])

end
subplot(1,3,[2 3])
scatter(age_sub,MEM_slope(:,1))
set(gca,'fontsize',16,'yticklabel','')
xlabel('Age')
ylim([-0.02 0.05])
%set(gcf,'position',[228 839 432 209])
set(gcf,'position',[228 416 441 609])
lsline
%% TEOAE

close all
for s=setdiff(1:length(d),[22,33])
    
    % Load scaped data
    load([d(s).name])
    figure(1)
    subplot(6,8,s)
    teoae_amp = dataalm.teoae{1};
    teoae_resp = [dataalm.teoae{3},dataalm.teoae{2}];
    stimear = dataalm.stim.ffr.ear(1);
    if isnan(teoae_amp(:,1,stimear))
        if stimear ==1
            stimear = 2;
        else
            stimear = 1;
        end
    end
    
    bar(squeeze(teoae_amp(:,1,stimear)))
    % save
    teoae_sub_amp(s,:,:,:) = squeeze(teoae_amp(:,:,stimear));
    teoae_sub_resp(s,:,:) = teoae_resp;
    %set(gca,'xtick',[1:5],'xticklabel',{'.5-1.5','1.5-2.5','2.5-3.5','3.5-4.5','4.5-5.5'})
    xtickangle(45)
    hold on
    plot(squeeze(teoae_amp(:,2,stimear)));
    plot([0 6],[-10 -10],'k--')
    plot(find(teoae_amp(:,3,stimear)),ones(size(find(teoae_amp(:,3,stimear))))*20,'kx')
    ylim([-25 25])
    if ~isempty(dataalm.per)
        title(dataalm.per.PersonNumber)
    else
        title('UH?')
    end
    figure(2)
    subplot(6,8,s)
    plot(teoae_resp(:,1),teoae_resp(:,[stimear+1]));
    ylim([-0.8 0.8])
    %xlim([0 13])
    box off
    if ~isempty(dataalm.per)
        title(dataalm.per.PersonNumber)
    else
        title('UH?')
    end

end%


%% mean
close all
bar(squeeze(mean(teoae_sub_amp(:,:,1),1)))
hold on
sig = zeros(size(teoae_sub_amp,1),5);
for s=1:size(teoae_sub_amp,1)
    jit_s = randn*0.1;
    for kk=1:5
        if teoae_sub_amp(s,kk,3)
            p_sig=plot(kk+jit_s,squeeze(teoae_sub_amp(s,kk,1)),'k+')
            sig(s,kk)=1;
        else
            p_nonsig=plot(kk+jit_s,squeeze(teoae_sub_amp(s,kk,1)),'ko')
        end
    end
end
        
nf=errorbar(1:5,mean(teoae_sub_amp(:,:,2),1),std(teoae_sub_amp(:,:,2),1)/sqrt((length(d)-2)),'r')
ylim([-35 30])

% %significant
for kk=1:5
    this_sig = length(find(sig(:,kk)))/size(sig,1);
    text(kk-0.1,25,[num2str(round(this_sig,2)*100) '%'])
end
set(gca,'fontsize',16)

legend([p_sig,p_nonsig,nf],{'Significant','N.S','Noise Floor'})
    

%% mean wave form

plot(squeeze(teoae_sub_resp(1,:,1)),teoae_sub_resp(:,:,3),'color',[0.5 0.5 0.5])
hold on
plot(squeeze(teoae_sub_resp(1,:,1)),mean(teoae_sub_resp(:,:,3),1),'k')
set(gca,'fontsize',16)

%plot(squeeze(teoae_sub_resp(1,:,1)),teoae_sub_resp(:,:,2))


%% functions
function fig = plot_audiogram_groups(idx,age_sub,aud,aud_freq,cmap)
for s=1:length(idx)
            p1 = semilogx(aud_freq,aud(idx(s),:)','-','color',[cmap(age_sub(idx(s))-17,:) 0.8]);
            plot_aud_param(p1,aud_freq);
            hold on
            
            cb=colorbar;
            
            cb.FontSize = 12;
            cb.Limits = [0 1];
            cb.Ticks = [linspace(0,1,5)];
            cb.TickLabels = {linspace(18,70,5)};
            cb.Label.String = 'Age';
            cb.Label.Rotation = 90;
            cb.Label.FontSize = 16;
            cb.Label.FontName = 'Arial';
            colormap(cmap)
end
set(gcf,'position',[305 412 432 299]);
fig = gcf
fig=gcf;
end



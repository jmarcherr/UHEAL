%% plot clinical measures

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


%% get clin data

for s=1:length(d)
   results{s} =  get_data(d,s);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% group audiogram plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%init
ages = [18 77];
colormap
cmap = cbrewer('seq','YlGnBu',max(ages)-17);
cmap(find(cmap>1))=1;cmap(find(cmap<0))=0;
ageidx = linspace(min(ages),max(ages),size(cmap,1));

close all
% get subject parameters
for s=1:length(d)
    age_sub(s) = results{s}.age_sub;
    aud(s,:,:) = results{s}.aud;
    aud_freq = results{s}.aud_freq;
    PTA_HF(s) = nanmean(results{s}.aud(7:end)); % PTA 9-16k
    PTA_LF(s) = nanmean(results{s}.aud(1:6));   % PTA 250-8k
    CP_sub(s,:) = results{s}.CP_sub;
    gender_sub(s,:) = results{s}.gender_sub;
    HV_sub(s) = results{s}.HV_sub;
    subid(s) = str2double(results{s}.sub_id(3:5));
end
uheal_data = struct;
uheal_data.subid = subid';
uheal_data.Age = age_sub';
uheal_data.CP = CP_sub;
uheal_data.gender = gender_sub;
uheal_data.DRCMR = HV_sub';
uheal_data.PTA_lf = PTA_LF';
uheal_data.PTA_hf = PTA_HF';

%%
idx = find(~isnan(age_sub));
fig = plot_audiogram_groups(idx,age_sub,aud,aud_freq,cmap)
% save figure
title(['all, n=' num2str(length(idx))])
saveas(fig,'figs/audiogram_all','epsc')
%% only NH
close all
NH_idx = setdiff(find(CP_sub==0));
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
%% age plots (only NH)
% all
close all
figure('renderer','painter')
[N,X]=hist(age_sub(NH_idx),max(age_sub(NH_idx))-min(age_sub(NH_idx)))
Bh = bar(X,N,'facecolor',[0.5 0.5 0.5],'EdgeColor', [0.5 0.5 0.5],'BarWidth',0.5);
xlabel('Age')
ylabel('n')
ax = ancestor(gca, 'axes');
xrule = ax.XAxis;

hold on

plot([25 25],[0 8],'k--')
plot([60 60],[0 8],'k--')
set(gca,'fontsize',16)
xlim([16 70])
xrule.FontSize = 14;
%set(gcf,'position',[305 412 432 299])
set(gcf,'position',[305 403 191 299])
box off
%title(['all, n=' num2str(length(idx))])
fig= gcf
saveas(fig,'figs/age_hist_all','epsc')

%%
%% age plots sex divded (only NH)
% all
close all
figure('renderer','painter')
maleidx = find(gender_sub==2 & CP_sub==0);
femaleidx =find(gender_sub==1 & CP_sub==0);

[N,X]=hist(age_sub(maleidx),6);
[Nf,Xf] = hist(age_sub(femaleidx),6)%max(age_sub(NH_idx))-min(age_sub(NH_idx)));
%Bh = bar([Xf],[Nf;N],'facecolor',[0.8 0.5 0.5],'EdgeColor', [0.8 0.5 0.5],'BarWidth',0.5);
%hold on
Bh = bar([Xf],[N;Nf],'stacked','BarWidth',.8);
Bh(1).FaceColor = [0.5 0.5 0.5];Bh(1).EdgeColor = [0.5 0.5 0.5];
Bh(2).FaceColor = [0.8 0.5 0.5];Bh(2).EdgeColor = [0.8 0.5 0.5];
xlabel('Age')
ylabel('n')
ax = ancestor(gca, 'axes');
xrule = ax.XAxis;

hold on

%plot([25 25],[0 8],'k--')
%plot([60 60],[0 8],'k--')
set(gca,'fontsize',16)
xlim([16 70])
xrule.FontSize = 14;
%set(gcf,'position',[305 412 432 299])
set(gcf,'position',[305 403 191 299])
hleg = legend('Male','Female');
hleg.FontSize = 10;
hleg.Box = 'off';
hleg.Position = [0.4634    0.7620    0.5236    0.1288];
box off
%title(['all, n=' num2str(length(idx))])
fig= gcf
saveas(fig,'figs/age_hist_male_female','epsc')


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
ylim([0 70])
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
set(gcf,'position',[305 412 432 299])
ylim([0 max(b1(1).YData)+max(b1(2).YData)+5]);
ymax = max(b1(1).YData)+max(b1(2).YData)+5;
fig= gcf
saveas(fig,'figs/age_hist_gender_HI','epsc')


%% grouped male female and HI from hvidovre
for i=1

% find who has been to hvidovre (from who-what-when 10-05-21)
DRCMR_idx = zeros(length(gender_sub),1);
DRCMR_idx = find(HV_sub),

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
fig= gcf
saveas(fig,'figs/age_hist_gender_DRCMR','epsc')
end
%% MEMR
close all
cmap = cbrewer('seq','YlGnBu',7+2);
cmap = cmap(2:end,:);

for s=1:length(d)
    reflex_sub(s,:,:) = results{s}.memr.reflex_sub;
    growth_sub_alt(s,:) = results{s}.memr.growth_sub;
    MEM_slope(s,:) =results{s}.memr.MEM_slope;
    disp(num2str(s))
end
    levels = results{1}.memr.levels;
    f_center = results{1}.memr.f_center
    freq = results{1}.memr.freq;
    
 uheal_data.memr = MEM_slope(:,1);   
    
    
%% average plots

close all
figure('Renderer','painter')
subplot(1,3,2)
cmap = cbrewer('seq','YlGnBu',7+2);
cmap = cmap(2:end,:);

% mean
MEMR_mean = squeeze(nanmean(reflex_sub,1));
% growth curve

p_y = plot(levels,nanmean(growth_sub_alt,1),'-k^','markeredgecolor','k','markerfacecolor','w');
hold on
for ll=1:length(levels)
    errorbar(levels(ll),nanmean(growth_sub_alt(:,ll),1),nanstd(growth_sub_alt(:,ll),1)/sqrt(length(d)),'color',cmap(ll,:))
    eby = plot(levels(ll),nanmean(growth_sub_alt(:,ll),1),'^','color',cmap(ll,:),'MarkerFaceColor',cmap(ll,:));
end

hold on
set(gca,'fontsize',12)
box off
set(gcf,'position',[228 839/2 432/2 209])
set(gca,'Xtick',levels,'Xticklabels',{'','80','','90','','100',''},'ytick',[0,.2,.4,.6],'yticklabels',{'0','.2','.4','.6'})
xlim([72 108])
xlabel('Elicitor level [dB]')
ylabel('\Sigma |\Delta Absorbance|')
% mean results
%figure(2)
subplot(1,3,1)
for i=1:length(levels)
    % Plot
    semilogx(f_center, nanmean(reflex_sub(:,:,i),1),'color',cmap(i,:));
    axis([freq(1) freq(end) -0.08 0.08])
    a=gca;
    a.XTick = [250,500,1000,2000,4000,8000];
    a.XTickLabel = [{'.25'},{'.5'},{'1'},{'2'},{'4'},{'8'}];
    hold on
    set(gca,'fontsize',12);
    box off
    for x = 1:length(levels)
        txt(x) = {[num2str(levels(x)),' dB SPL']};
    end
    set(gcf,'position',[228 839/2 432 209]);
    hleg = legend(txt);
    hleg.Box = 'off';
    hleg.Position = [0.5547 0.1231 0.3704 0.8038];
    xlabel('Frequency [Hz]');
    ylabel('\Delta Absorbance');
end
set(gcf,'Position',[228 420 618 209]);
hleg.Position= [0.6509 0.1470 0.2589 0.8038];

fig= gcf
saveas(fig,'figs/memr_summary','epsc')
%% MEMR slope

close all

figure('Renderer','painter')
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
set(gca,'xtick',[],'fontsize',12,'Xcolor','none')
box off
ylim([-0.02 0.05])

end
close all
figure('Renderer','painter')
%subplot(1,3,[2 3])
scatter(age_sub,MEM_slope(:,1),'o','markerfacecolor',[0.6 0.6 0.6],'MarkerEdgecolor','k')
set(gca,'fontsize',12)
ylabel('Linear coefficient')

xlabel('Age')
ylim([-0.02 0.05])
%set(gcf,'position',[228 839 432 209])
%set(gcf,'position',[228 416 441 609])
set(gcf,'Position',[228 420 618/3 209]);
ll=lsline
set(ll,'linewidth',2,'color','k')

fig= gcf
saveas(fig,'figs/memr_age','epsc')
%% TEOAE

close all
for s=1:length(d)
    
    % Load scaped data
    %load([d(s).name])
    teoae_sub_amp(s,:,:) = results{s}.teoae.amp;
    teoae_sub_resp(s,:,:,:) = results{s}.teoae.resp;
end
%%
%% mean
close all
% valid measurements
figure('Renderer','painter')
idx = find(~isnan(teoae_sub_amp(:,1,1,1)));
bar(squeeze(nanmean(teoae_sub_amp(idx,:,1),1)),'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5])
hold on
sig = zeros(size(teoae_sub_amp(idx,:,:,:),1),5);
for s=1:size(teoae_sub_amp(idx,:,:,:),1)
    jit_s = randn*0.1;
    for kk=1:5
        if teoae_sub_amp(idx(s),kk,3)
            p_sig=plot(kk+jit_s,squeeze(teoae_sub_amp(idx(s),kk,1)),'o',...
                'color',cmap(age_sub(s)-17,:),'markerfacecolor',cmap(age_sub(s)-17,:),...
                'markersize',4)
            sig(s,kk)=1;
        else
            p_nonsig=plot(kk+jit_s,squeeze(teoae_sub_amp(idx(s),kk,1)),'+',...
                'color',cmap(age_sub(s)-17,:),'markersize',4)
        end
    end
end
        
nf=errorbar(1:5,nanmean(teoae_sub_amp(:,:,2),1),nanstd(teoae_sub_amp(:,:,2),1)/sqrt((length(idx))),'r')
%ylim([-35 30])

% %significant
for kk=1:5
    this_sig = length(find(sig(:,kk)))/size(sig,1);
    text(kk-0.1,27,[num2str(round(this_sig,2)*100) '%'])
end
set(gca,'fontsize',16)

hleg = legend([p_sig,p_nonsig,nf],{'Significant','N.S','Noise Floor'})
hleg.Position = [0.2044    0.2280    0.3013    0.1829];
hleg.Box = 'off';
hleg.FontSize = 12;

ylabel('Amplitude (dB SPL)')
xlabel('Frequency bands (+/- .5 kHz)')
set(gcf,'Position',[441 386 468 339])
box off
fig= gcf
saveas(fig,'figs/teoae_overview','epsc')

%% alternative
close all
ages = [18 77];
colormap
cmap = cbrewer('seq','YlGnBu',max(ages)-17);
cmap(find(cmap>1))=1;cmap(find(cmap<0))=0;
ageidx = linspace(min(ages),max(ages),size(cmap,1));
% valid measurements
figure('Renderer','painter')
idx = find(~isnan(teoae_sub_amp(:,1,1,1)) & CP_sub==0);
bar(squeeze(nanmean(teoae_sub_amp(idx,:,1),1)),'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5])
hold on
sig = zeros(size(teoae_sub_amp(idx,:,:,:),1),5);
for s=1:size(teoae_sub_amp(idx,:,:,:),1)
    jit_s = randn*0.1;
    %for kk=1:5
        %if teoae_sub_amp(idx(s),kk,3)
            p_sig=plot(1:5,squeeze(teoae_sub_amp(idx(s),:,1)),'color',[cmap(age_sub(s)-17,:)])
            sig(s,kk)=1;
        %else
            %p_nonsig=plot(kk+jit_s,squeeze(teoae_sub_amp(idx(s),kk,1)),'ko')
        %end
    %end
end
        
nf=errorbar(1:5,nanmean(teoae_sub_amp(:,:,2),1),nanstd(teoae_sub_amp(:,:,2),1)/sqrt((length(idx))),'r')
ylim([-35 30])

% %significant
%for kk=1:5
%    this_sig = length(find(sig(:,kk)))/size(sig,1);
%    text(kk-0.1,27,[num2str(round(this_sig,2)*100) '%'])
%end
set(gca,'fontsize',16)

%hleg = legend([p_sig,p_nonsig,nf],{'Significant','N.S','Noise Floor'})
%hleg.Position = [0.2044    0.2280    0.3013    0.1829];
%hleg.Box = 'off';
%hleg.FontSize = 12;

ylabel('Amplitude (dB SPL)')
xlabel('Frequency bands (+/- .5 kHz)')
set(gcf,'Position',[441 386 468 339])
%% mean wave form
close all
figure('Renderer','painter')
%plot(squeeze(teoae_sub_resp(1,:,1)),teoae_sub_resp(:,:,3),'color',[0.5 0.5 0.5])
hold on
plot(squeeze(teoae_sub_resp(1,:,1)),nanmean(teoae_sub_resp(:,:,3),1),'k')
shadedErrorBar(squeeze(teoae_sub_resp(1,:,1)),nanmean(teoae_sub_resp(:,:,3),1),nanstd(teoae_sub_resp(:,:,3),1)/sqrt(size(teoae_sub_resp,1)))
set(gca,'fontsize',16)
set(gcf,'Position',[441 552 468 173])
xlabel('Time [ms]')
ylabel('mPA')
xlim([4 12.5])
box off
%plot(squeeze(teoae_sub_resp(1,:,1)),teoae_sub_resp(:,:,2))
fig= gcf
saveas(fig,'figs/teoae_waveform','epsc')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reverse digit span, nesi, ssq
close all
for s=1:length(d)
    % rds
    if ~isfield(results{s}.rds,'Fresp')
        fds(s) = nan;
    else
        fds(s) = sum([results{s}.rds.Fresp.Total_score{:}]);
    end
    if ~isfield(results{s}.rds,'Bresp')
        rds(s) = nan;
    else
        rds(s) = sum([results{s}.rds.Bresp.Total_score{:}])
    end
    %nesi
    if isempty(results{s}.nesi)
        nesi(s) = nan;
    else
    nesi(s) = results{s}.nesi;
    end
    if isempty(results{s}.tts);
        tts(s) = nan;
    else
        tts(s)= results{s}.tts;
    end
    % ssq
    if isempty(results{s}.ssq)
        ssq(s,:) = nan(1,12);
    else
    ssq(s,:)= normalize(results{s}.ssq(:,2));
    end
end
%% save to uheal struct
uheal_data.rds = rds';
uheal_data.fds = fds';
uheal_data.nesi =nesi';
uheal_data.tts = tts';
uheal_data.ssq12 = nanmean(ssq,2);

%%
figure('Renderer','painter')
scatter(age_sub,rds,'o','markerfacecolor',[0.6 0.6 0.6],'MarkerEdgecolor','k')
ll=lsline
title('Digit-span')
xlabel('Age')
ylabel('Digit-span score')
set(ll,'linewidth',2,'color','k')
set(gca,'fontsize',12)
set(gcf,'position',[441 459 270 266])

fig= gcf
saveas(fig,'figs/rds_age','epsc')

figure('Renderer','painter')
scatter(age_sub,nesi,'o','markerfacecolor',[0.6 0.6 0.6],'MarkerEdgecolor','k')
ll=lsline
title('NESI')
xlabel('Age')
ylabel('NESI score')
set(ll,'linewidth',2,'color','k')
set(gca,'fontsize',12)
set(gcf,'position',[441 459 270 266])

fig= gcf
saveas(fig,'figs/NESI_age','epsc')

figure('Renderer','painter')
for ii=1:12
    subplot(3,4,ii)
scatter(age_sub,ssq(:,ii)','o','markerfacecolor',[0.6 0.6 0.6],'MarkerEdgecolor','k')
ll=lsline
title(['SSQ-12 - ' num2str(ii)])
xlabel('Age')
ylabel('SSQ score')
set(ll,'linewidth',2,'color','k')
set(gca,'fontsize',12)
set(gcf,'position',[184 106 1049 596])
end
fig= gcf
saveas(fig,'figs/ssq_age','epsc')

% tts
figure('Renderer','painter')
scatter(age_sub,tts,'o','markerfacecolor',[0.6 0.6 0.6],'MarkerEdgecolor','k')
ll=lsline
title('TTS')
xlabel('Age')
ylabel('TTS score')
set(ll,'linewidth',2,'color','k')
set(gca,'fontsize',12)
set(gcf,'position',[441 459 270 266])

fig= gcf
saveas(fig,'figs/tts_age','epsc')
%%


%% acalos
fs = [500,1000,2000,4000];
x= [-10:2.5:120];
addpath('/work1/jonmarc/UHEAL_master/UHEAL/_scripts/_tools/ACALOS')
for s=1:length(d)
    if ~isempty(results{s}.acalos);
        tmp_data = results{s}.acalos.raw.RawData;
        freq_idx = tmp_data(:,1);
        
        for ff=1:4
        this_f = find(freq_idx==fs(ff));
        acalos_f{ff} = tmp_data(this_f,:,:);
        dataLoudFit = reshape(acalos_f{ff}(:,2:3)', 1, []);
        
        select_fitting_function = 'BTUX';
        
        [fitparams,rawData] = fit_loudness_function(dataLoudFit,select_fitting_function);
        data_tmp(ff,2)=fs(ff);
        data_tmp(ff,[3,4,5])=fitparams;
        %% Estimated parameters
        
        [y, failed] = loudness_function_bh2002(x, fitparams);
        iZero = (max(find(y==0)));
        iUCL = (min(find(y==50)));
        %%%%%%%%%%%%%%%%%%%%% HTL
        if  ~isempty(iZero)
            HTL(ff)=x(iZero+1);
        else
            if ~isempty(x(min(find(y<=1))+1))
                HTL(ff) =  x(min(find(y<=1))+1);
            else
                HTL(ff) = x(find(min(y)));
            end
        end
        %%%%%%%%%%%%%%%%%%%%% UCL
        if  ~isempty(iUCL)
            UCL(ff)=x(iUCL-1);
        else
            UCL(ff) = x(y==max(y));
        end
        
        Lcut =fitparams(1);
        m_lo = fitparams(2);
        %        HTL_A(ii) = Lcut - 22.5/m_lo;
        HTL_A(ff)= fitparams(1)-22.5/m_lo;
        MCL(ff) = x(max(find(y<=25)));
        
        y_plot(:,ff) = y;
        x_plot(:,ff) = x;
        
    end
    % save results
    acalos_yplot(s,:,:) = y_plot;
    HTL_A_sub(s,:) = HTL_A;
    MCL_sub(s,:) = MCL;
    AC_slope(s,:) = data_tmp(:,4)';
    L_cut(s,:) = data_tmp(:,3)';
    m_high(s,:) = data_tmp(:,5)';
    %acalos_ = results{s}.acalos.proc;
    else
        acalos_yplot(s,:,:) = nan(size(acalos_yplot(1,:,:)));
        HTL_A_sub(s,:) = nan(size(HTL_A_sub(1,:)));
        MCL_sub(s,:) = nan(size(MCL_sub(1,:)));
        AC_slope(s,:) = nan(size(AC_slope(1,:)));
        L_cut(s,:) = nan(size(L_cut(1,:)));
        m_high(s,:) = nan(size(m_high(1,:)));
    end
    clc
    disp(['subject ' num2str(s) ' processed'])
end

uheal_data.acalos_slope = AC_slope;
%% plotting
close all
YNH = find(age_sub<25);
ONH = find(age_sub>45);
figure('Renderer','painter')
plot(x,squeeze(nanmean(acalos_yplot,1)),'linewidth',2)
xlim([0 120])
xlabel('dB SPL')
ylabel('CU')
box off
hleg = legend(num2str(fs(:)),'location','best');
hleg.Box = 'off'
set(gca,'fontsize',12)
set(gcf,'position',[150 592 251 410])
fig= gcf
saveas(fig,'figs/ACALOS_all','epsc')

figure('Renderer','painter')
for ii=1:4
subplot(1,4,ii)
plot(x,squeeze(nanmean(acalos_yplot(YNH,:,ii))),'k')
hold on
plot(x,squeeze(nanmean(acalos_yplot(ONH,:,ii))),'r')
xlim([0 105])
xlabel('dB SPL')
ylabel('CU')
box off
title([num2str(fs(ii)) ' Hz'])
set(gcf,'position',[681 838 560 164])
end
hleg = legend('Young','Older');
hleg.Box = 'off'
fig= gcf
saveas(fig,'figs/ACALOS_yvo','epsc')

figure('Renderer','painter')
for ii=1:4

subplot(1,4,ii)
scatter(age_sub,AC_slope(:,ii),'o','markerfacecolor',[0.6 0.6 0.6],'MarkerEdgecolor','k')
ll=lsline
title([num2str(fs(ii)) ' Hz'])
xlabel('age')
ylabel('Slope (CU/dB)')
ylim([0.2 0.8])
box off
set(ll,'linewidth',2,'color','k')
set(gcf,'position',[681 838 560 164])
hold on
scatter(age_sub(find(CP_sub)),AC_slope(find(CP_sub),ii),'rx')
end
fig= gcf
saveas(fig,'figs/ACALOS_slope','epsc')

%% save uheal data
uheal_table = struct2table(uheal_data)
writetable(uheal_table,'uheal_data_clin.csv')  
save('uheal_data_clin.mat','uheal_data')

%% functions
function fig = plot_audiogram_groups(idx,age_sub,aud,aud_freq,cmap)
    figure('renderer','painter')
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



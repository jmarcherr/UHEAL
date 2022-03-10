clear;close all;
cd(fileparts(matlab.desktop.editor.getActiveFilename))
% find all processed files
addpath('/work1/jonmarc/UHEAL_master/UHEAL/_scripts/_tools');
addpath('/work1/jonmarc/UHEAL_master/UHEAL/_scripts/_tools/cbrewer/cbrewer');
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
    sub_id{s} = d(s).name(1:5)
    % load subject info
    %cd(clin_dir);
    %load([d(s).name(1:5) '.mat']);
    subinfo = data_abr.subinfo;
%     nr_reject(s,:) = data_abr.nr_reject;
%     valid_trials(s,:) = [size(data_abr.trials{1},2) size(data_abr.trials{2},2)];
    %cd(this_dir);
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
    if isempty(subinfo.CP)
        CP(s) = nan;
    else
        CP(s) = subinfo.CP;
    end
    
    % get nr. og trials
    % 9Hz
    for ii=1:2
    if size(data_abr.trials{ii},2)<1000
    n_trials(s,ii) = nan;%size(data_abr.trials{1},2);       
    else
    n_trials(s,ii) = size(data_abr.trials{ii},2);
    end
    end
        
    % get rejected
    nr_reject(s,:) = data_abr.nr_reject;
    
    cd(this_dir)
    %%
    clc
    disp(['9 Hz: avg nr trials, ' num2str(floor(nanmean(n_trials(:,1)))) ,...
        ' +/- ' num2str(nanstd(n_trials(:,1)))])
    disp(['40 Hz: avg nr trials, ' num2str(floor(nanmean(n_trials(:,2)))),...
        '+/- ' num2str(nanstd(n_trials(:,2)))])
    disp(['9 Hz: avg % reject, ' num2str((nanmean(nr_reject(:,1)))) ,...
        ' +/- ' num2str(nanstd(nr_reject(:,1)))])
    disp(['40 Hz: avg % reject, ' num2str((nanmean(nr_reject(:,2)))),...
        '+/- ' num2str(nanstd(nr_reject(:,2)))])
    
    % participants without abr
    disp([num2str(length(find(isnan(n_trials(:,1))))) ' subjects missing 9 Hz'])
    disp([num2str(length(find(isnan(n_trials(:,2))))) ' subjects missing 40 Hz'])
    %% Baseline correction
    %loop over rates (9,40)
    for kk=ccc
        data_w = data_abr.abr{kk};
        if length(data_w)<300
            data_w = zeros(327,1);
        end

        baseline = data_w(find(data_abr.time(data_abr.tidx)==0)+blIDX);%

        
        tn = data_abr.time(data_abr.tidx)-delay;
        tnIDX = find(tn>-5e-3);
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
            hold on addpath('/work1/jonmarc/UHEAL_master/UHEAL/_scripts/_tools/cbrewer');    
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

%% Correct for AD break
% plot both
sub_preAD = setdiff(1:20,2) % UH002 re-recorded
plot(squeeze(nanmean(sub_corrected(sub_preAD,1,1:244-13))))
hold on
plot(squeeze(nanmean(sub_corrected([2 21:end],1,14:end))))
sub_abr = [];
sub_abr(sub_preAD,:,:) = sub_corrected(sub_preAD,:,1:end-13);
sub_abr([2 21:size(sub_corrected,1)],:,:) = sub_corrected([2 21:end],:,14:end);

tn_corr = tnIDX(14:end);
% delayed time:
t_abr = tn(tn_corr);
close all
plot(t_abr,squeeze(nanmean(sub_abr,1))')

% baseline correction
baseline = squeeze(nanmean(sub_abr(:,:,find(t_abr>-1e-3 & t_abr<0)),3))'
for s=1:length(d)
sub_abr_b(s,:,:) = squeeze(sub_abr(s,:,:)-baseline(:,s)');
end

%% missing abr:
missing_abr     = find(isnan(sub_abr(:,1,1)))';
missing_age     = find(isnan(age));
missing_lab     = find(isnan(lab));
missing_gender  = find(isnan(gender));
missing_CP      = find(isnan(CP));

%% visual rejection: Find artifact subjects
reject_vis = 0;
if reject_vis
    for s=1:length(d)
        close all
        fig=figure
        plot(t_abr,squeeze(sub_abr_b(s,2,:))','color',[0.5 0.5 0.1])
        hold on
        plot(t_abr,squeeze(sub_abr_b(s,1,:))','b')
        xlabel('Time [s]')
        ylabel('Amplitude [\muV]')
        xlim([-2e-3 8e-3])
        grid on
        title(['subid: ' sub_id{s}])
        % key press for reject, click for keep
        keydown = waitforbuttonpress;
        if (keydown == 0)
            disp('Mouse button was pressed');
            disp('reject')
            rjt_sub(s) = 1;
        else
            disp('Key was pressed');
            rjt_sub(s) = 0;
        end
        %pause()
        savefile = [sub_id{s}];
        saveas(fig,['subject_ABRs' filesep savefile],'epsc')
        
        
    end
    %%
    cd('subject_ABRs')
    save('rjt_sub.mat','rjt_sub')
    cd ..
else
    
    %% load rjt_sub
    load('subject_ABRs/rjt_sub.mat')
end


%% peak detect
close all
peakd=1;
if peakd
        load('peaks/WVpos.mat');load('peaks/WIpos.mat');load('peaks/SP.mat');
    for s=1:length(d)%46:length(d)
        
        plot(t_abr,squeeze(sub_abr_b(s,2,:))','color',[0.5 0.5 0.1])
        hold on
        plot(t_abr,squeeze(sub_abr_b(s,1,:))','b')
        xlim([-2e-3 8e-3])
        title(['subid: ' sub_id{s}])
        grid on
        hold off
        if isnan(sub_abr(s,1,1)) | rjt_sub(s)==1
            SP(s,:) = [nan nan];
            WIneg(s,:) = [nan nan];
            WIpos(s,:) = [nan nan];
            WVneg(s,:) = [nan nan];
            WVpos(s,:) = [nan nan];
        else
            %[SP(s,:)]=getABRwaves(gcf,'SP');
            %[WIpos(s,:)]=getABRwaves(gcf,'Wave I pos');
            [WIneg(kk,:)]=getABRwaves(gcf,'Wave I neg');
            %[WVpos(s,:)]=getABRwaves(gcf,'Wave V pos');
            [WVneg(kk,:)]=getABRwaves(gcf,'Wave V neg');
        end
        
    end
    %%
    cd('peaks')
    mkdir(date)
    cd(date)
    save('WVpos.mat','WVpos');save('WIpos.mat','WIpos');save('SP.mat','SP');save('WIneg.mat','WIneg');save('WVneg.mat','WVneg');
    cd ..
    cd ..
    %%
else
    load('peaks/WVpos.mat');load('peaks/WIpos.mat');load('peaks/SP.mat');
end

%%
close all
savefile = 1;
% latencies and amplitudes
cm = cbrewer('qual','Set1',10)
cmap = cm([1 2 10],:);

%load('peaks/WVpos.mat');load('peaks/WIpos.mat');load('peaks/SP.mat');
NH_idx = find(CP==0);NH_idx(end)=[];
NH_idx = setdiff(NH_idx,[find(rjt_sub),missing_abr]);
figure('Renderer','painter')
subplot(1,3,1)
splot(1)=scatter(SP(NH_idx,1),SP(NH_idx,2),'o',...
    'MarkerEdgeColor','k','markerfacecolor',[cmap(1,:)],'MarkerFaceAlpha',.5)
hold on
splot(2)=scatter(WIpos(NH_idx,1),WIpos(NH_idx,2),'o'...
    ,'MarkerEdgeColor','k','markerfacecolor',cmap(2,:),'MarkerFaceAlpha',.5)
hold on
splot(3)=scatter(WVpos(NH_idx,1),WVpos(NH_idx,2),'o'...
    ,'MarkerEdgeColor','k','markerfacecolor',cmap(3,:),'MarkerFaceAlpha',.5)

ylabel('Amplitude [\muV]');
xlabel('Time [s]')     
xlim([-1e-3 10e-3])
hleg = legend([splot(1) splot(2) splot(3)],'SP','AP','Wave V');
hleg.Position = [0.2831 0.7789 0.0877 0.1471];
%hleg.Box = 'off';

set(gca,'fontsize',12,'xtick',[0:1:10]*1e-3)
%set(gca,'TickLabelInterpreter','latex');

%latencies
subplot(1,3,2)
scatter(age(NH_idx),SP(NH_idx,1),'o',...
        'MarkerEdgeColor','k','markerfacecolor',cmap(1,:),'MarkerFaceAlpha',.5)
hold on
scatter(age(NH_idx),WIpos(NH_idx,1),'o',...
        'MarkerEdgeColor','k','markerfacecolor',cmap(2,:),'MarkerFaceAlpha',.5)
scatter(age(NH_idx),WVpos(NH_idx,1),'o',...
        'MarkerEdgeColor','k','markerfacecolor',cmap(3,:),'MarkerFaceAlpha',.5)
    
ll=lsline
idx = [3,2,1]
for ii=1:3
set(ll(ii),'color',cmap(idx(ii),:),'Linewidth',2)
end

ylabel('Latency [s]');
xlabel('Age')     
xlim([0 100])
set(gca,'fontsize',12)
%set(gca,'TickLabelInterpreter','latex');


% amplitudes
subplot(1,3,3)
scatter(age(NH_idx),SP(NH_idx,2),'o',...
    'MarkerEdgeColor','k','markerfacecolor',cmap(1,:),'MarkerFaceAlpha',.5)
hold on
scatter(age(NH_idx),WIpos(NH_idx,2),'o',...
    'MarkerEdgeColor','k','markerfacecolor',cmap(2,:),'MarkerFaceAlpha',.5)
scatter(age(NH_idx),WVpos(NH_idx,2),'o',...
    'MarkerEdgeColor','k','markerfacecolor',cmap(3,:),'MarkerFaceAlpha',.5)

ll=lsline
idx = [3,2,1]
for ii=1:3
set(ll(ii),'color',cmap(idx(ii),:),'Linewidth',2)
end

ylabel('Amplitude [\muV]');
xlabel('Age')     
xlim([0 100])
set(gca,'fontsize',12)
%set(gca,'TickLabelInterpreter','latex');


% save

set(gcf,'position',[441 324 867 401])
if savefile
fig=gcf;
saveas(fig,"figs/abr_scatter3",'epsc')
end

%% save to UHEAL_Data file
load('/work1/jonmarc/UHEAL_master/UHEAL/UHEAL_data/scraped/uheal_data_table/uheal_data.mat')
uheal_data.SP_amp = nan(size(uheal_data.subid));
uheal_data.SP_lat = nan(size(uheal_data.subid));
uheal_data.AP_amp = nan(size(uheal_data.subid));
uheal_data.AP_lat =  nan(size(uheal_data.subid));
uheal_data.WV_amp = nan(size(uheal_data.subid));
uheal_data.WV_lat =  nan(size(uheal_data.subid));
%%
for s=1:length(SP)
    % get this subid
    thisID = str2double(sub_id{s}(3:5))
    this_idx = find(uheal_data.subid==thisID);
    uheal_data.SP_amp(this_idx) = SP(s,2);
    uheal_data.SP_lat(this_idx) = SP(s,1);
    uheal_data.AP_amp(this_idx) = WIpos(s,2);
    uheal_data.AP_lat(this_idx) = WIpos(s,1);
    uheal_data.WV_amp(this_idx) = WVpos(s,2);
    uheal_data.WV_lat(this_idx) = WVpos(s,1);
    
end

uheal_table = struct2table(uheal_data)

writetable(uheal_table,'uheal_data.csv')  
save('uheal_data.mat','uheal_data')
%%
%ratios WI/WV
close all
IVratio = WVpos(NH_idx,2)./WIpos(NH_idx,2);
idx = find(IVratio<5 & IVratio>-3);
figure
scatter(age(idx),IVratio(idx))
lsline
% ratio SP/AP
figure
spapratio = SP(:,2)./WIpos(:,2)
idx = find(spapratio<1 & spapratio>0);
scatter(age(idx),SP(idx,2)./WIpos(idx,2))
lsline
%lsline
%ylim([-1 1])
         
%% plot average with errorbars
save = 1;
close all
% all subjects
figure('Renderer','painter')
[fig1a,fig1b]=plot_abr(t_abr,sub_abr_b,age)
%suptitle(['all, n=' num2str(length(find(~isnan(sub_abr(:,1,1)))))])

% only NH
figure('Renderer','painter')
sub_abr_NH = sub_abr_b(find(CP==0),:,:);
[fig2a,fig2b]=plot_abr(t_abr,sub_abr_NH,age(CP==0))
%suptitle(['NH, n=' num2str(length(find(~isnan(sub_abr_NH(:,1,1)))))])


% only NH young
figure('Renderer','painter')
sub_abr_YNH = sub_abr_b(find(CP==0 & age<25),:,:);
[fig3a,fig3b]=plot_abr(t_abr,sub_abr_YNH,age(find(CP==0 & age<25)))
%suptitle(['YNH, n=' num2str(length(find(~isnan(sub_abr_YNH(:,1,1)))))])

% only NH old
figure('Renderer','painter')
%subplot(1,3,[1 2])
sub_abr_ONH = sub_abr_b(find(CP==0 & age>45),:,:);
[fig4a,fig4b]=plot_abr(t_abr,sub_abr_ONH,age(find(CP==0 & age>45)))
%suptitle(['ONH, n=' num2str(length(find(~isnan(sub_abr_ONH(:,1,1)))))])



% only HI
% figure
% sub_abr_HI = sub_abr(find(CP==1),:,:);
% plot_abr(t_abr,sub_abr_HI)
% title(['HI, n=' num2str(length(find(~isnan(sub_abr_HI(:,1,1)))))])



if save
    %waveforms
saveas(fig1a,"figs/abr_ONH",'epsc')
saveas(fig2a,"figs/abr_YNH",'epsc')
saveas(fig3a,"figs/abr_NH",'epsc')
saveas(fig4a,"figs/abr_all",'epsc')
    %hist
saveas(fig1b,"figs/abr_ONH_h",'epsc')
saveas(fig2b,"figs/abr_YNH_h",'epsc')
saveas(fig3b,"figs/abr_NH_h",'epsc')
saveas(fig4b,"figs/abr_all_h",'epsc')    
end
%% Functions



function [fig1,fig2]=plot_abr(t_abr,sub_abr,age)
cm = cbrewer('qual','Set1',10)
cmap = cm([1 2 10],:);
rate_colors = {'k','b'};
%subplot(1,3,[1 2])
% loop over rates
for kk=1:2
    
    abr_var=squeeze(nanstd(sub_abr(:,kk,:)))/sqrt(length(find(~isnan(sub_abr(:,1,1)))));
    abr_mean = squeeze(nanmean(sub_abr(:,kk,:)))';
    p(kk)=plot(t_abr,abr_mean,'color',rate_colors{kk});
    shadedErrorBar(t_abr,abr_mean,abr_var,...
    'lineprops',['-' rate_colors{kk}],'transparent',1);
    hold on
end
        xlim([-.5e-3 8e-3])%changed from 6e-3 to 8e-3
        ylim([-.2 .4])
        hold on
        box off
        grid on
        
        
 plot(t_abr,zeros(size(t_abr)),'k--');
   

hleg = legend([p(1) p(2)],'9 Hz','40Hz');
hleg.Box = 'off';

set(gca,'fontsize',12,'xtick',[0:1:8]*1e-3)
set(gca,'TickLabelInterpreter','latex');
ylabel('ABR Amplitude [\muV]');
xlabel('Time [s]')  
set(gcf,'position',[293 510 315 215])
fig1=gcf;
%axis tight
figure
%subplot(1,3,3)
hist(age,100)
xlim([0 99])
ylim([0 8])
xlabel('Age')
ylabel('n')
set(gca,'YAxisLocation','left','fontsize',12,'xtick',[25,50,75])
set(gcf,'position',[612 510 174 215])
 box off
 set(gca,'TickLabelInterpreter','latex');

fig2=gcf;
end

% get peaks
function [output] = getABRwaves(fig,titstring)
%fig = figure;
%plot(x,y,'Color',[0,0.7,0.9])
title(titstring)
% Enable data cursor mode
datacursormode on
dcm_obj = datacursormode(fig);
% Set update function
set(dcm_obj,'UpdateFcn',@myupdatefcn)
% Wait while the user to click
disp('Click line to display a data tip, then press "Return"')
pause
% Export cursor to workspace
info_struct = getCursorInfo(dcm_obj);
if isfield(info_struct, 'Position')
    disp('Clicked positioin is')
    disp(info_struct.Position)
end
output =info_struct.Position

    function [output_txt, output_value] = myupdatefcn(~,event_obj)
        % ~            Currently not used (empty)
        % event_obj    Object containing event data structure
        % output_txt   Data cursor text
        pos = get(event_obj, 'Position');
        output_txt = {['x: ' num2str(pos(1))], ['y: ' num2str(pos(2))]};
        output_value = [pos(1) pos(2)];
    end


%xwave = output_value(1);
%ywave = output_value(2);
end
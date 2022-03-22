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
ccc = 1;
%%
for s=1:length(d)%[26,11,28,61]%
    load([d(s).name])
    this_dir = cd;
    sub_id{s} = d(s).name(1:5)
    % load subject info
    %cd(clin_dir);
    %load([d(s).name(1:5) '.mat']);
    subinfo = data_abr.subinfo;
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
    %n_trials(s,1) = size(data_abr.trials{1},2);
    %n_trials(s,2) = size(data_abr.trials{2},2);
    % get rejected
    nr_reject(s) = data_abr.nr_reject(1);
    fs = data_abr.fs;
    cd(this_dir)
    %% Baseline correction
    %loop over rates (9,40)
    for kk=ccc
        data_w = data_abr.abr{kk};
        if length(data_w(1,:))<300
            data_w = zeros(16,327);
        end

        baseline = data_w(:,find(data_abr.time(data_abr.tidx)==0)+blIDX);%

        
        tn = data_abr.time(data_abr.tidx)-delay;
        tnIDX = find(tn>-10e-3);
        %figure(s)
        if plotting
            subplot(11,10,s)
        end
        if mean(data_w==0)
            sub_corrected(s,:,:) = nan(16,size(sub_corrected,3));
        else
            
            sub_corrected(s,:,:) =data_w(:,tnIDX)-baseline;
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
close all
sub_preAD = 1:20
plot(squeeze(nanmean(sub_corrected(sub_preAD,1,1:end))))
hold on
plot(squeeze(nanmean(sub_corrected(21:end,1,4:end))))
sub_abr = [];
sub_abr(1:20,:,:) = sub_corrected(sub_preAD,:,1:end-4);
sub_abr(21:size(sub_corrected,1),:,:) = sub_corrected(21:end,:,5:end);
%%
%sub_abr = sub_corrected;
tn_corr = tnIDX(5:end);
% delayed time:
t_abr = tn(tn_corr);
%%
%% missing abr:
missing_abr     = find(isnan(sub_abr(:,1,1)))';
missing_age     = find(isnan(age));
missing_lab     = find(isnan(lab));
missing_gender  = find(isnan(gender));
missing_CP      = find(isnan(CP));
%% mean over all
close all
plot(t_abr,squeeze(nanmean(sub_abr,1))')
hold on
%shadedErrorBar(t_abr,squeeze(nanmean(sub_abr,1)),squeeze(nanstd(sub_abr,1))/sqrt(length(find(~isnan(sub_abr(:,1,1))))))
xlim([-.002 .060])
%ylim([-1 0.5])
box off
ylabel('amplitude \muV')
xlabel('time s')
% baseline correction
baseline = squeeze(nanmean(sub_abr(:,:,find(t_abr>-1e-3 & t_abr<0)),3))'
for s=1:length(d)
sub_abr_b(s,:,:) = squeeze(sub_abr(s,:,:));%;-baseline(:,s)');

end
savefile = 'midlate_all'
fig = gcf;
%set(gcf,'position',[681 859 313 165])
%saveas(fig,savefile,'epsc')

%% mean over NH
NH_idx = find(CP==0);NH_idx(end)=[];
NH_idx = setdiff(NH_idx,[missing_abr]);


%%
close all
% individual traces
plot(t_abr,-squeeze(nanmean(sub_abr(NH_idx,:,:),2))','color',[0.5 0.5 0.5 0.5])
hold on
plot(t_abr,-squeeze(nanmean(mean(sub_abr(NH_idx,:,:)),2)))


%% with filtering
close all
filt_coef = [10 100];
filt_def = designfilt('bandpassfir','FilterOrder',40, ...
    'CutoffFrequency1',filt_coef(1),'CutoffFrequency2',filt_coef(2), ...
    'SampleRate',fs);

mlr_filt = filtfilt(filt_def,squeeze(mean(sub_abr(NH_idx,:,:),2))');

% individual traces
%for ff=4
    filt_coef = [10 1000];
    filt_def = designfilt('bandpassfir','FilterOrder',40, ...
        'CutoffFrequency1',filt_coef(1),'CutoffFrequency2',filt_coef(2), ...
        'SampleRate',fs);
    
    mlr_filt = filtfilt(filt_def,squeeze(sub_abr(NH_idx,1,:))');
%subplot(1,4,ff)
plot(t_abr*1000,-squeeze(mlr_filt)','color',[0.5 0.5 0.5 0.2])
hold on
plot(t_abr*1000,-mean(mlr_filt'),'linewidth',2)
xlim([-.002 .060]*1000)
ylim([-2 2])
set(gca,'fontsize',14,'xtick',[0 0.025 0.05 0.075]*1000)
box off
ylabel('amplitude \muV')
xlabel('time ms')

%end
%set(gcf,'position',[681 642 313 360])
savefile = 'figs/midlate_ind'
fig = gcf;
%set(gcf,'position',[681 859 313 165])
%saveas(fig,savefile,'epsc')
%%
%% visual rejection: Find artifact subjects
for s=1:length(d)
    close all
    fig=figure
    plot(t_abr,squeeze(sub_abr(s,1,:))','color',[0.5 0.5 0.1])
    %hold on
    %plot(t_abr,squeeze(sub_abr_b(s,1,:))','b')
    xlabel('Time [s]')
    ylabel('Amplitude [\muV]')
    xlim([-10e-3 110e-3])
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
    saveas(fig,['subject_MLRs' filesep savefile],'epsc')
    

end
%%
cd('subject_MLRs')
save rjt_sub
cd ..
%% peak detect (manual)
close all
for s=96:length(d)%
    
    plot(t_abr,squeeze(sub_abr(s,2,:))','color',[0.5 0.5 0.1])
    hold on
    %plot(t_abr,squeeze(sub_abr_b(s,1,:))','b')
    plot(t_abr,mean(mlr_filt'),'linewidth',2) % mean wave-form
    xlabel('Time [s]')
    ylabel('Amplitude [\muV]')
    xlim([-10e-3 110e-3])
    title(['subid: ' sub_id{s}])
    grid on
    hold off
    if isnan(sub_abr(s,1,1))
        N0(s,:) = [nan nan];
        P0(s,:) = [nan nan];
        Na(s,:) = [nan nan];
        Pa(s,:) = [nan nan];
    else
         [N0(s,:)]=getABRwaves(gcf,'N0');
         [P0(s,:)]=getABRwaves(gcf,'P0');
         [Na(s,:)]=getABRwaves(gcf,'Na');
         [Pa(s,:)]=getABRwaves(gcf,'Pa');

    end

end

%% 
save =0;
if save
    cd('peaks')
    save N0
    save P0
    save Na
    save Pa
    cd ..
end

%% plotting
load('peaks/N0');load('peaks/P0');load('peaks/Na');load('peaks/Pa')
close all
figure('renderer','painter')
scatter(N0(:,1),N0(:,2))
hold on
scatter(P0(:,1),P0(:,2))
scatter(Na(:,1),Na(:,2))
scatter(Pa(:,1),Pa(:,2))
plot(t_abr,mean(mlr_filt'),'linewidth',2) % mean wave-form
legend('N0','P0','Na','Pa','Mean MLR')
xlabel('time (s)')
ylabel('amplitude \muV')
set(gcf,'position',[441 279 473 423]);
fig = gcf;
    savefile = 'peak_picking';
    saveas(fig,['figs' filesep savefile],'epsc')

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
%set(gca,'TickLabelInterpreter','latex');
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
 %set(gca,'TickLabelInterpreter','latex');

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
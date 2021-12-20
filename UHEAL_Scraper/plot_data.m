%% audiograms for UHEAL WP2
clear all
close all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
cd ..
addpath('_analysis')
datadir = '_clin_data';
cd(datadir)
d = dir('*.mat');subjects = char(d.name);

%% Audiometry
offset = 50;
cmap = flip(cbrewer('div','RdYlBu',70-17));
figure(1)
for s=setdiff(1:length(subjects)-1,8)
    
    % Load scaped data
    load([subjects(s,:)])
    
    % Audiometry
    [age,aud_L,aud_R,aud_freq] = get_aud(dataalm);
    
    %p1 = plot(1:length(aud_freq),nanmean([aud_L,aud_R]'),'-','color',cmap(age-17,:));
    p1 = semilogx(aud_freq,nanmean([aud_L,aud_R]'),'-','color',cmap(age-17,:));
    hold on
    
    
    
end
plot_aud_param

%% MEMR

figure(2)
cmap = cbrewer('seq','YlGnBu',7+2);
cmap = cmap(2:end,:);
for s=1:length(subjects)-1
    
    % Load scaped data
    load([subjects(s,:)])
    
    stimear = dataalm.stim.ffr.ear(1);
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
        subplot(5,5,s)
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
            title(dataalm.per.PersonNumber)
        end
        
    else
        warning('No MEMR data for subject')
        subplot(5,5,s)
        a=gca;
        a.XTick = [250,500,1000,2000,4000,8000];
        a.XTickLabel = [{'.25'},{'.5'},{'1'},{'2'},{'4'},{'8'}];
        %title(['\Delta absorbance = contracted - baseline. Target pressure: ', num2str(target_p), ' daPa'])
        xlabel('Frequency [kHz]')
        ylabel('\Delta Absorbance')
        grid on
        hold on
        title(dataalm.per.PersonNumber)
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
    errorbar(levels(ll),nanmean(growth_sub_alt(:,ll),1),nanstd(growth_sub_alt(:,ll),1)/sqrt(length(subjects)),'color',cmap(ll,:))
    eby = plot(levels(ll),nanmean(growth_sub_alt(:,ll),1),'^','color',cmap(ll,:),'MarkerFaceColor',cmap(ll,:))
end

hold on
set(gca,'fontsize',16)
box off
set(gcf,'position',[680   949   316/2   149])
set(gca,'Xtick',levels,'Xticklabels',{'','80','','90','','100',''},'ytick',[0,.2,.4,.6],'yticklabels',{'0','.2','.4','.6'})

% mean results

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
    set(gcf,'position',[680   949   316   149])
    hleg = legend(txt)
    hleg.Box = 'off'
end


%% TEOA
close all
for s=1:length(subjects)-1
    
    % Load scaped data
    load([subjects(s,:)])
    figure(1)
    subplot(5,5,s)
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
    %set(gca,'xtick',[1:5],'xticklabel',{'.5-1.5','1.5-2.5','2.5-3.5','3.5-4.5','4.5-5.5'})
    xtickangle(45)
    hold on
    plot(squeeze(teoae_amp(:,2,stimear)));
    plot([0 6],[-10 -10],'k--')
    plot(find(teoae_amp(:,3,stimear)),ones(size(find(teoae_amp(:,3,stimear))))*20,'kx')
    ylim([-25 25])
    title(dataalm.per.PersonNumber)
    figure(2)
    subplot(5,5,s)
    plot(teoae_resp(:,1),teoae_resp(:,[stimear+1]));
    ylim([-0.8 0.8])
    %xlim([0 13])
    box off
    title(dataalm.per.PersonNumber)

end

%% acalos


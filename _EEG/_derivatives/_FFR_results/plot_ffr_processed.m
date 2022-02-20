%% plot FFR results
clear all;close all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
thisdir = cd;
cd ..
cd ..
cd ..
ft_defaults
UHEAL_startup
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
        age(s) =data_ffr.subinfo.age;
        gender(s) = data_ffr.subinfo.gender;
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
        
    else
        
        subinfo{s} = data_ffr.subinfo;
        age(s) = nan;
    end
    
    cd(thisdir)
    
end
 
%% NH idx
NH_idx =(CP==0 & ~isnan(age));
%% plotting time series
close all
cidx = [17,18]; %[5,11]%
for s=1:length(TS)
    if length(chans{s})>=18 & F{s}(cidx(stimear(s)))>= F_crit{s}(1)%18 % all chans
        chanoi = cidx(stimear(s));
    TS_sub(s,:,:) = TS{s}(chanoi,:);
    fft_sub(s,:,:) = f_fft{s}(chanoi,:,:);
    else
        TS_sub(s,:,:) = nan;
        fft_sub(s,:,:) = nan;
    end
end

plot(nanmean(squeeze(TS_sub(age<25,:,:))),'k')
hold on
plot(nanmean(squeeze(TS_sub(age>45,:,:))),'r')


%% plotting spectra
figure
plot(fft_freq(1,:),nanmean(squeeze(db(fft_sub(age<25,:,:)))),'k')
hold on
plot(fft_freq(1,:),nanmean(squeeze(db(fft_sub(age>45,:,:)))),'r')
%% plotting FFR amplitude and SNR
close all
cmap = flip(cbrewer('div','RdYlBu',77-17)); cmap = cmap./max(cmap(:));
jit = randn(size(FFR)).*0.2;
cidx = [17,18]; %[5,11]%

for s=1:length(FFR)
    subplot(1,2,1)
    if length(chans{s})>=18 & F{s}(cidx(stimear(s)))>= F_crit{s}(1)%18 % all chans
        chanoi = cidx(stimear(s));
    FFR_sub(s,:) = FFR{s}(chanoi);
    FFR_SNR_sub(s,:) = FFR_SNR{s}(chanoi);
    plot(1+jit(s),db(FFR_sub(s,:))','o','markerfacecolor',[cmap(age(s)-17,:)],'color',[cmap(age(s)-17,:)])
    %plot(1+jit(s),db(FFR(s))','o','markerfacecolor',[0.6 0.6 0.6],'color','k')
    xlim([0 2])
    hold on
    xlabel('')
    ylabel('FFR amplitude (dB mV)')
    set(gca,'xtick',[],'fontsize',16,'Xcolor','none')
    box off
 
    subplot(1,2,2)
    plot(1+jit(s),FFR_SNR_sub(s,:)','o','markerfacecolor',[cmap(age(s)-17,:)],'color',[cmap(age(s)-17,:)])
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
 
NH_idx =(CP==0 & ~isnan(age) & ~isnan(FFR_sub)');
figure
subplot(2,2,1)
scatter(age(NH_idx),db(FFR_sub(NH_idx)))
[rho,pval]=corr(age(NH_idx)',db(FFR_sub(NH_idx)))
lsline 
subplot(2,2,2)
scatter(age(NH_idx),db(FFR_SNR_sub(NH_idx)))
[rho,pval]=corr(age(NH_idx)',FFR_SNR_sub(NH_idx))
lsline
subplot(2,2,3)
%bar([1 2],[nanmean(db(FFR_sub(find(age<25)))) nanmean(db(FFR_sub(find(age>45))))])
ynh = [];onh=[];
ynh = db(FFR_sub(find(age<=25)))';
onh = db(FFR_sub(find(age>45)))';
if length(ynh)<length(onh)
    ynh(length(ynh)+1:length(onh)) = nan;
else
    onh(length(onh)+1:length(ynh)) = nan;
end
boxplot([ynh;onh]')
subplot(2,2,4)
bar([1 2],[nanmean(FFR_SNR_sub(find(age<25))) nanmean(FFR_SNR_sub(find(age>45)))])
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

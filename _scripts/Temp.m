%% plot audiogram

close all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
cd ..

UHEAL_startup
cd(datadir)
cd([datadir '/scraped']);
freq_aud = [250 500 1000 2000 4000 8000 9000 10000 11200 12500 14000 16000];
d=dir('*.mat')
% colormap
cmap = flip(cbrewer('div','RdYlBu',70-17));
figure(1)
for s=1:length(d)
    
    %load relevant file
    %load(['UH' num2str(subid(s))])
    
    load([d(s).name])
    if ~isempty(dataalm.aud)
        %get audiogram and age
        [age,aud_L,aud_R,aud_freq,gender] = get_aud(dataalm);
        
        %get stim ear
        stimear = dataalm.stim.abr.ear(1);
        %plot stim ear
        if stimear ==1 %left ear
            
            p1 = semilogx(aud_freq,aud_L','-','color',[cmap(age-17,:) 0.8]);
            plot_aud_param(p1,aud_freq)
            hold on
        else % right ear
            p1 = semilogx(aud_freq,aud_R','-','color',[cmap(age-17,:) 0.8]);
            plot_aud_param(p1,aud_freq)
            hold on
        end
        
        age_sub(s) = age; %log age
        gender_sub(s) = gender;
    else
        age_sub(s) = nan;
        gender_sub(s) = nan;
    end
   
    
end
cb=colorbar

cb.FontSize = 12;
cb.Limits = [0 1]
cb.Ticks = [linspace(0,1,5)];
cb.TickLabels = {linspace(18,70,5)};
cb.Label.String = 'Age';
cb.Label.Rotation = 90;
cb.Label.FontSize = 16;
cb.Label.FontName = 'Arial';
%colormap(cmap)

set(gcf,'position',[305 412 432 299])

%% age plot
figure(2)
hist(age_sub,max(age_sub)-min(age_sub))
xlabel('Age (years)')
ylabel('nr. of participants')
hold on
plot([25 25],[0 4],'k--')
plot([60 60],[0 4],'k--')
set(gca,'fontsize',16)
xlim([16 70])
set(gcf,'position',[305 412 432 299])

%% grouped male female
figure(3)
Y   = gender_sub(find(age_sub<=25)); 
O1  = gender_sub(find(age_sub>25 & age_sub<=60));
O2  = gender_sub(find(age_sub>60));

b1=bar([1:3],[length(find(Y==1)) length(find(O1==1)) length(find(O2==1));...
    length(find(Y==2)) length(find(O1==2)) length(find(O2==2))]','stacked')

hold on
xlabel('Group')
ylabel('Nr. of participants')
set(gca,'xtick',[1:3],'xticklabels',{'Y','O1','O2'},'fontsize',16)
legend(b1,'Female','Male')
b1(1).FaceColor = [0.2 0.2 0.5]
b1(2).FaceColor = [0.5 0.7 0.5]
set(gcf,'position',[305 412 432 299])


%%
cd(fileparts(matlab.desktop.editor.getActiveFilename))
cd ..

UHEAL_startup
cd(datadir)
cd([datadir '/scraped']);
d=dir('*.mat')


for s=1:length(d)
  
    load([d(s).name])
    if ~isempty(dataalm.nesi)
        %get audiogram and age
        nesi(s)=load(dataalm.nesi)
   
    
end


%%
%cd('/Users/jmarcher/Documents/git')
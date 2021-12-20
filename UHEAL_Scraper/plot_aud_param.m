function plot_aud_param(plotid,aud_freq)
fsize = 16;
pw = [680   523   325   282];
set(gcf,'position',pw)
set(gca,'Ydir','reverse','xtick',aud_freq([1:6 end]),'xticklabel',aud_freq([1:6 end])/1000);
set(gca,'ytick',[0 20 50 100])
ylabel('Hearing level (dB)')
xl=xlabel('Frequency (kHz)')
ax = gca;
ax.YGrid= 'on'

set(gca,'Fontsize',fsize);

ax = ancestor(plotid, 'axes');
xrule = ax.XAxis;

xlim([200 20000])
xtickangle(45)
set(gca,'Fontsize',fsize);
box off
xrule.FontSize = 14;
xl.FontSize = 16;
hold on
plot([8e3 8e3],[120 -30],'k--')

ylim([-20 110])

% cb=colorbar
% 
% cb.FontSize = 12;
% cb.Limits = [0 1]
% cb.Ticks = [linspace(0,1,5)];
% cb.TickLabels = {linspace(18,70,5)};
% cb.Label.String = 'Age';
% cb.Label.Rotation = 90;
% cb.Label.FontSize = 16;
% cb.Label.FontName = 'Arial';
%colormap(cmap)
end
%% data plotting

%% MEMR
%filename = 'JM_R_MEMR.mat';
%load(filename);
%MEMR response
close all
ipsi = data.reflex_ipsi.response;
levels = data.reflex_ipsi.labels;
target_p = data.data_3DT.TPP;
freq = data.reflex_ipsi.freq;
f_center = data.reflex_ipsi.f_center;
figure(3)
        for i=1:length(levels)
            % Plot
            semilogx(f_center, ipsi(:,i))
            axis([freq(1) freq(end) -0.1 0.1])
            a=gca;
            a.XTick = [250,500,1000,2000,4000,8000];
            a.XTickLabel = [{'250'},{'500'},{'1000'},{'2000'},{'4000'},{'8000'}];
            title(['\Delta absorbance = contracted - baseline. Target pressure: ', num2str(target_p), ' daPa'])
            xlabel('Frequency [Hz]')
            ylabel('\Delta Absorbance')
            grid on
            hold on
            for x = 1:length(levels)
                txt(x) = {[num2str(levels(x)),' dB SPL']};
            end
            legend(txt)
        end
        
        
        
        %% 3D tymp
            figure(4)
            surf(log(data.data_3DT.freq), data.data_3DT.pressure, 100*data.data_3DT.absorbance)
            title(['3DT plot - TTP: ',num2str(data.data_3DT.TPP),' daPa'])
            g1=gca;
            set(g1,'YDir','reverse')
            g1.XTick = log([250,500,1000,2000,4000,8000]);
            g1.XTickLabel = [{'250'},{'500'},{'1k'},{'2k'},{'4k'},{'8k'}];
            xlabel('Frequency [Hz]')
            ylabel('Pressure [daPa]')
            zlabel('Absorbance [%]')
            axis([log(226) log(8000) -300 200 0 100])
            
            
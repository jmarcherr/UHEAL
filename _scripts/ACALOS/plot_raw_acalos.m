%find raw data
close all
d = dir('*raw.dat')

acalos_raw = load(d.name)
% find processed data
d = dir('*nbr.dat')

acalos_proc = load(d.name)

fc = unique(acalos_raw(:,1));
%%
for fcs = 1:length(fc)
    idx = find(acalos_raw(:,1)==fc(fcs));
figure(1)
plot(acalos_raw(idx,2,:),acalos_raw(idx,3,:),'o')
hold on
figure(2)
for vars = 1:3
subplot(1,3,vars)
plot(acalos_proc(fcs,vars+2),'o')
hold on
end

end
xlabel('level (dB SPL)')
ylabel('CU')
legend(num2str(fc(:)))

% find processed data
d = dir('*nbr.dat')

acalos_proc = load(d.name)
fc = unique
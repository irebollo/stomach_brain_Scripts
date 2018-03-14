% heart_getHRVts

%% parameters
cfgMain = global_getcfgmain;

%%

% subjects = global_subjectList


subj_idx = 44


data_ECG_cutted_filename = global_filename(subj_idx,cfgMain,'data_ECG_cutted')
load(data_ECG_cutted_filename)

%  figure;plot(dataECG_cutted(3,:))
[HeartBeats, R_time] = heart_peak_detect(dataECG_cutted(3,:),5000);

% figure
% hist(IBI)

%% Correct 


timeResampled =0:1/5000:898;

Rpeak_sample = 0
for i=1:length(HeartBeats)
Rpeak_sample(i) = HeartBeats(i).R_sample;
end

for iBeat = 2:length(R_time)
IBI(iBeat-1)=R_time(iBeat)-R_time(iBeat-1);
end
% 
figure
hist(IBI)

figure
plot(Rpeak_sample)
% 
% Rpeak_sample(find(IBI>1.1))
% 
% 
% figure
% plot(dataECG_cutted)
% hold on
% plot(Rpeak_sample,2000,'o')
% 
% xlim([4344860 4364860])
% 
% % Rpeak_sample = [Rpeak_sample,hb1.DataIndex,hb2.DataIndex,hb3.DataIndex,hb4.DataIndex,hb5.DataIndex,hb6.DataIndex,hb7.DataIndex,hb8.DataIndex,hb9.DataIndex,hb10.DataIndex,hb11.DataIndex,]
% Rpeak_sample = [Rpeak_sample,hb1.DataIndex]
% Rpeak_sample_corrected = sort(Rpeak_sample) 

[t,ibi,t_int,ibi_int, F, PSD] = IBI_PSD(Rpeak_sample,5000,1,30,50,1024)

figure
plot(ibi_int)

HRV_timeseries_filename = global_filename(subj_idx,cfgMain,'HRV_timeseries')
save(HRV_timeseries_filename,'dataECG_cutted','ibi_int','t_int')


%%

% 
% close all
% 
% 
% % data_ECG_cutted_filename = global_filename(subj_idx,cfgMain,'data_ECG_cutted')
% 
% subj_idx = 36
% data_ECG_cutted_filename = global_filename(subj_idx,cfgMain,'data_ECG_cutted')
% load(data_ECG_cutted_filename)
% 
% % figure
% % plot(dataECG_cutted(1,:))
% 
% data =  dataECG_cutted(1,:);
% [HeartBeats, R_time] = heart_peak_detect(data,5000);
% 
% 
% 
% % Correct 
% 
% timeResampled =0:1/5000:898;
% 
% Rpeak_sample = []
% for i=1:length(HeartBeats)
% Rpeak_sample(i) = HeartBeats(i).R_sample;
% end
% 
% IBI = []
% for iBeat = 2:length(R_time)
% IBI(iBeat-1)=R_time(iBeat)-R_time(iBeat-1);
% end
% 
% 
% [t,ibi,t_int,ibi_int, F, PSD] = IBI_PSD(Rpeak_sample,5000,1,30,50,1024)
% 
% figure
% plot(ibi_int)
% 
% % HRV_timeseries_filename = global_filename(subj_idx,cfgMain,'HRV_timeseries')
% 
% HRV_timeseries_filename = ['E:\HRV\dataHeart\IBIts_S' ,sprintf('%.2d',subj_idx), '.mat']
% 
% save(HRV_timeseries_filename,'dataECG_cutted','ibi_int','t_int')







subjects = [8 9 10 13 15 16 17 18 19 20 21 22 23 24 25 26 29 31 32 33 34 35 36 39 40 41 43]



spectrumHRV_All_subjects = zeros(length(subjects),243);


centerFrequencyHRVAllSubjects = zeros(length(subjects),2);

for iSubj = 1:length(subjects)
    
%% Loop
subj_idx = subjects(iSubj);

rootDir= strcat(global_path2root,'\Subjects\');



HRV_timeseries_filename = global_filename(subj_idx,cfgMain,'HRV_timeseries');
load(HRV_timeseries_filename)


timeseries2Coherence = ibi_int';
load(strcat(global_path2root,'\scripts4paper\files\sampleFieldtripStruc.mat'))
labelsChannelsMAIN = {'HRV'};
labelsChannels = labelsChannelsMAIN;
clusterRegionsComparisons = timeseries2Coherence;
dataStructure.hdr = EGG_downsampled.hdr;
dataStructure.fsample = 1;
dataStructure.time{1,1}  = [0:1/dataStructure.fsample:(size(clusterRegionsComparisons,2))-1];
dataStructure.label = labelsChannels;%channelStr;
dataStructure.cfg = EGG_downsampled.cfg;
dataStructure.trial{1,1} = clusterRegionsComparisons;
dataStructure.sampleinfo = [1 length(timeseries2Coherence)]

% a data.sampleinfo field in your data structure
% > that is consistent with the actual data. In general, data.sampleinfo
% > is an Nx2 matrix containing, for each trial in the data, the indices
% > of the begin and end sample of that trial, with respect to the
% > original data set (on disk)

disp('Resampling...')
cfg = [];  %initialize configuration structure
cfg.detrend = 'no'; % remove linear trend from the data (done per trial)
cfg.demean = 'yes';
cfg.resamplefs= 0.5; % 4 x top-freq (15 cpm = 0.25 Hz) - Nyquist = 30 cpm  frequency at which the data will be resampled
HRV_downsampled = ft_resampledata(cfg,dataStructure); % This procedure also lowpass filter the data at half the new sr 


HRV_downsampled.trial{1,1} = ft_preproc_polyremoval (HRV_downsampled. trial{1,1}, 2);
HRV_downsampled.trial{1,1} = ft_preproc_standardize (HRV_downsampled. trial{1,1});

data = [HRV_downsampled.trial{1,1}(cfgMain.beginCut:cfgMain.endCut)];


nVoxels = size(data,1);

%load base fieldtrip data structure
load(strcat(global_path2root,'\scripts4paper\files\sampleFieldtripStruc.mat'))

% Define fieldtrip structure
channelStr=cell(nVoxels,1);
for iVoxel = 1:nVoxels
    channelList(iVoxel,1) = iVoxel;
    channelStr(iVoxel) = cellstr(mat2str(iVoxel));
end

dataStructure.hdr = EGG_downsampled.hdr;
dataStructure.fsample = 0.5;
dataStructure.time{1,1}  = [0:2:(size(data,2)*2)-1];
dataStructure.label = channelStr;
dataStructure.cfg = EGG_downsampled.cfg;
dataStructure.trial{1,1} = data;

cfgWelch = [];
cfgWelch.keeptrials = 'no';
cfgWelch.lengthWindow = 120;
cfgWelch.overlap = 6;

len = dataStructure.fsample*cfgWelch.lengthWindow; % length of subtrials cfg.length s in samples
dataStructure.sampleinfo=[1 max(dataStructure.time{1,1})*dataStructure.fsample];
cfg = [];
cfg.trl(:,1) = dataStructure.sampleinfo(1):(len/cfgWelch.overlap):dataStructure.sampleinfo(2)-len+1;%trial start in samples from begining of raw data
cfg.trl(:,2) = dataStructure.sampleinfo(1)+len-1:(len/cfgWelch.overlap):dataStructure.sampleinfo(2);%trial ends in samples from begining of raw data
cfg.trl(:,3) = 0; %offset of the trigger with respect to the trial
data_trials = ft_redefinetrial(cfg,dataStructure);


% Estimate spectrum
cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.output = 'pow';
cfg.pad = 1000;
cfg.foilim = [1/cfgWelch.lengthWindow 0.25]; % 0 - 6 cpm
cfg.keeptrials = cfgWelch.keeptrials;




frequencyWelch = ft_freqanalysis(cfg,data_trials);

spectrumHRV_All_subjects(iSubj,:) = frequencyWelch.powspctrm;

filterWidth = 0.03
indexFrequencies = find (frequencyWelch.freq >= 0.0666 & frequencyWelch.freq <= 0.1333); % narrow LF HRV
maxPower= max(frequencyWelch.powspctrm(indexFrequencies));% from 0.033 to 0.066 hz
        maxPowerLocation = frequencyWelch.powspctrm==maxPower;

maxFrequency = frequencyWelch.freq(find(maxPowerLocation));
centerFrequencyHRVAllSubjects(iSubj,2) = maxFrequency;
centerFrequencyHRVAllSubjects(iSubj,1) = subj_idx;


        indexFrequenciesFilter = find (frequencyWelch.freq >= maxFrequency-filterWidth & frequencyWelch.freq <= maxFrequency+filterWidth);
        filterPlot=zeros(1,length(frequencyWelch.freq));
        filterPlot(indexFrequenciesFilter)= maxPower;
                   
spectrumHeartPlot = figure
plot (frequencyWelch.freq,filterPlot,'b','lineWidth',3)
hold on
plot(frequencyWelch.freq,frequencyWelch.powspctrm,'r','lineWidth',3)


plotDir = strcat (global_path2subject(subj_idx),'PreprocessingLog',filesep);
plotFilename = strcat(plotDir,'S_',sprintf('%.2d',subj_idx),'_spectrumHeartPlot');


xlabel('Frequency')
ylabel('Power')
title(['S',sprintf('%.2d',subj_idx),32,'HRV power'],'fontsize',18)

set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gcf, 'PaperPositionMode', 'auto');

print ('-dpng', '-painters', eval('plotFilename'))
print ('-depsc2', '-painters', eval('plotFilename'))
saveas(spectrumHeartPlot,strcat(plotFilename,'.fig'))   


end

filenameAllPeaksHRV =[global_path2root filesep 'scripts4paper' filesep 'files' filesep 'HRV_peaks_info.mat']

save(filenameAllPeaksHRV,'centerFrequencyHRVAllSubjects')

%% Load and plot rotations all subjects

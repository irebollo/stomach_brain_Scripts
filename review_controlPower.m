% review_controlPower

insideBrain= tools_getIndexBrain('inside');

subjects = global_subjectList

cfgMain = global_getcfgmain;


peak = global_getEGGpeaks

SpectrumXSignalXSubject=zeros(30,2,93); % subjects, cominations freqsbins
SpectrumXSignalXSubjectPeak = zeros (30,2,61);

%%
DMNcenter_tal = [-5 -49 40];% tal Fox 2005 PNAS
DMNcenter = tal2mni(DMNcenter_tal');  % transform to MNI
DMN_center_vox = tools_mni2vox(DMNcenter,1); % transform to voxel coordinates 
emptyBrain = zeros (53,63,46); % empty 3d Volume for storing empirical PLV
% take cube of 3 voxels centered at coodinates
emptyBrain(DMN_center_vox(1)-1:DMN_center_vox(1)+1,DMN_center_vox(2)-1:DMN_center_vox(2)+1,DMN_center_vox(3)-1:DMN_center_vox(3)+1 ) = 1;
DMNCoreMask = emptyBrain(insideBrain);
DMNCoreMask = find(DMNCoreMask);
% tools_writeMri(emptyBrain,'DMN_seed')

S1_center_MNI = [45 -28 46];
S1_center_vox = tools_mni2vox(S1_center_MNI',1); % transform to voxel coordinates 
emptyBrain = zeros (53,63,46); % empty 3d Volume for storing empirical PLV
emptyBrain(S1_center_vox(1)-1:S1_center_vox(1)+1,S1_center_vox(2)-1:S1_center_vox(2)+1,S1_center_vox(3)-1:S1_center_vox(3)+1 ) = 1;
SICoreMask = emptyBrain(insideBrain);
SICoreMask = find(SICoreMask);

% tools_writeMri(emptyBrain,'SI_seed')


for iSubj = 1:length(subjects)
    subj_idx = subjects(iSubj)
%% Load BOLD

boldname = global_filename(subj_idx,cfgMain,'BOLD_filtered_fullbandFilename');
load(boldname)


BOLD_filtered_zscored_cutted_inside = BOLD_filtered_zscored(insideBrain,cfgMain.beginCut:cfgMain.endCut);
mean_precuneus_timeseries = mean(BOLD_filtered_zscored_cutted_inside(DMNCoreMask,:));
mean_ssi_timeseries = mean(BOLD_filtered_zscored_cutted_inside(SICoreMask,:));

figure
plot(mean_precuneus_timeseries)
hold on
plot(mean_ssi_timeseries,'r')

%% get fourier

timeseries2Coherence(1,:) = mean_ssi_timeseries;
timeseries2Coherence(2,:) = mean_precuneus_timeseries;
nVoxels = size(timeseries2Coherence,1);

%load base fieldtrip data structure
load(strcat(global_path2root,'\scripts4paper\files\sampleFieldtripStruc.mat'))

labelsChannelsMAIN = {'S1','DMN'};
labelsChannels = labelsChannelsMAIN;
clusterRegionsComparisons = timeseries2Coherence;
dataStructure.hdr = EGG_downsampled.hdr;
dataStructure.fsample = 0.5;
dataStructure.time{1,1}  = [0:2:(size(clusterRegionsComparisons,2)*2)-1];dataStructure.label = labelsChannels;%channelStr;
dataStructure.cfg = EGG_downsampled.cfg;
dataStructure.trial{1,1} = clusterRegionsComparisons;

%computing fourier
cfgWelch = [];
cfgWelch.lengthWindow = 120; %seconds
cfgWelch.overlap = 6;% propotion i.e 1/6
len = dataStructure.fsample*cfgWelch.lengthWindow; % length of subtrials cfg.length s in samples
dataStructure.sampleinfo=[1 max(dataStructure.time{1,1})*dataStructure.fsample];

% Cut data into trials
cfg = [];
cfg.trl(:,1) = dataStructure.sampleinfo(1):(len/cfgWelch.overlap):dataStructure.sampleinfo(2)-len+1;%trial start in samples from begining of raw data
cfg.trl(:,2) = dataStructure.sampleinfo(1)+len-1:(len/cfgWelch.overlap):dataStructure.sampleinfo(2);%trial ends in samples from begining of raw data
cfg.trl(:,3) = 0; %offset of the trigger with respect to the trial
data_trials = ft_redefinetrial(cfg,dataStructure);

% Perform frequency analysis
cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.output = 'powandcsd';
cfg.pad = 1000;
cfg.foilim = [1/cfgWelch.lengthWindow 0.1]; % 0 - 6 cpm
cfg.keeptrials = 'no';

frequencyWelch = ft_freqanalysis(cfg,data_trials);

figure
plot(frequencyWelch.freq,frequencyWelch.powspctrm)
% labels('si','dmn')


indPeaks(iSubj) = find(frequencyWelch.freq == peak(iSubj,3));
SpectrumXSignalXSubject(iSubj,:,:) = frequencyWelch.powspctrm ; % subjects, cominations freqsbins
SpectrumXSignalXSubjectPeak(iSubj,:,:) = frequencyWelch.powspctrm(:,indPeaks(iSubj)-30:indPeaks(iSubj)+30);


%% Get power in BOLD

end

figure
plot(squeeze(mean(SpectrumXSignalXSubjectPeak))')
legend('s1','Prec')
% check cube is inside the brain


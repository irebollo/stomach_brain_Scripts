pupil_coherenceEGGpupil

% initialize coh all subjects

subjects = [8 9 10 13 16 17 18 19 20 21 22 23 26 29 32 33 34 35 36 39]
% load egg
% load pupil
% put them into ft
% compute coherence


%% 

peaksAllsubjects = global_getEGGpeaks;
cfgMain = global_getcfgmain;
%% Execute

CoherenceXSubjectsXCluster_pupil = zeros(length(subjects),93); % 13 =12 clusters + global, 93 = freqs bins
CoherenceXSubjectsXCluster_peakEGG_pupil = zeros(length(subjects),61); % 13 = 12 clusters + global, 60 = freqs bins
proportionNanWholeScan = zeros(length(subjects),1);

for iSubj = 1:length(subjects)
    
    
subj_idx = subjects(iSubj)

    
filename_EGG_amplitude = global_filename(subj_idx,cfgMain,'FilenameamplitudeXVolumeBestChannel_FULLBAND');
load(filename_EGG_amplitude)

pupil_data_bp_filename = ['Y:\Subjects\Subject' sprintf('%.2d',subj_idx) '\Timeseries\Pupil\Pupil_bp_MRI_'  sprintf('%.2d',subj_idx)];

load(pupil_data_bp_filename)

proportionNanWholeScan(iSubj) = propotion_nan_wholeScan;

% clusterTimeseries = clusterTimeseries(cfgMain.beginCut:cfgMain.endCut,:);

timeseries2Coherence = zeros(2,420); % 12 clusters plus EGG plus global
timeseries2Coherence(1,:) = EGG_FullBand(cfgMain.beginCut:cfgMain.endCut);
timeseries2Coherence(2,:) = data_pupil_bp(cfgMain.beginCut:cfgMain.endCut);
% timeseries2Coherence(3:14,:)= clusterTimeseries';
nVoxels = size(timeseries2Coherence,1);

%load base fieldtrip data structure
load(strcat(global_path2root,'scripts4paper\files\sampleFieldtripStruc.mat'))

labelsChannelsMAIN = {'EGG','Pupil'};
labelsChannels = labelsChannelsMAIN;
clusterRegionsComparisons = timeseries2Coherence;
dataStructure.hdr = EGG_downsampled.hdr;
dataStructure.fsample = 0.5;
dataStructure.time{1,1}  = [0:2:(size(clusterRegionsComparisons,2)*2)-1];
dataStructure.label = labelsChannels;%channelStr;
dataStructure.cfg = EGG_downsampled.cfg;
dataStructure.trial{1,1} = clusterRegionsComparisons;




%% computing fourier
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

% Get coherence
cfg            = [];
cfg.method     = 'coh';
coherence      = ft_connectivityanalysis(cfg, frequencyWelch);
    

    CoherenceXSubjectsXCluster_pupil(iSubj,:) = coherence.cohspctrm(1,:);
    
    indPeaks(iSubj) = find(coherence.freq == peaksAllsubjects(iSubj,3)); % find the frequency bin corrdesponding to the EGG peak of the subject
    CoherenceXSubjectsXCluster_peakEGG_pupil(iSubj,:) = coherence.cohspctrm(indPeaks(iSubj)-30:indPeaks(iSubj)+30);
end

figure
nhist(proportionNanWholeScan)

% subjects_clean = proportionNanWholeScan>0 & proportionNanWholeScan<0.16

figure
plot(mean(CoherenceXSubjectsXCluster_peakEGG_pupil))

figure
plot(coherence.freq,mean(CoherenceXSubjectsXCluster_pupil))

mean_squared_coherence_peak = mean(CoherenceXSubjectsXCluster_peakEGG(:,31))*mean(CoherenceXSubjectsXCluster_peakEGG(:,31))
stf_squared_coherence_peak = std(CoherenceXSubjectsXCluster_peakEGG(:,31))*std(CoherenceXSubjectsXCluster_peakEGG(:,31))

figure
nhist(CoherenceXSubjectsXCluster_peakEGG(:,31).*CoherenceXSubjectsXCluster_peakEGG(:,31))



%% coherence clusters

%% Execute

CoherenceXSubjectsXCluster = zeros(length(subjects),13,93); % 13 =12 clusters + global, 93 = freqs bins
CoherenceXSubjectsXCluster_peakEGG = zeros(length(subjects),13,61); % 13 = 12 clusters + global, 60 = freqs bins

for iSubj = 1:length(subjects)
    subj_idx = subjects(iSubj)
    %load
    clusterTimeseries_coherence_filename = global_filename(subj_idx,cfgMain,'clusterTimeseries_coherence_filename');
    load(clusterTimeseries_coherence_filename)
    CoherenceXSubjectsXCluster(iSubj,:,:) = coherence.cohspctrm(1:13,:);
    
    indPeaks(iSubj) = find(coherence.freq == peaksAllsubjects(iSubj,3)); % find the frequency bin corrdesponding to the EGG peak of the subject
    CoherenceXSubjectsXCluster_peakEGG(iSubj,:,:) = coherence.cohspctrm(1:13,indPeaks(iSubj)-30:indPeaks(iSubj)+30);
end

%% Average them separatly in absoulte frequency or relative frequency with respecto to EGG peak (frequency bin 31)

meanCoherenceXCluster_peak = squeeze(mean(CoherenceXSubjectsXCluster_peakEGG(:,:,31)));
serCoherenceXCluster_peak = squeeze(std (CoherenceXSubjectsXCluster_peakEGG(:,:,31))./sqrt(length(subjects)));

meanSquaredCoherenceXCluster_peak = squeeze(mean(CoherenceXSubjectsXCluster_peakEGG(:,:,31).*CoherenceXSubjectsXCluster_peakEGG(:,:,31)));
serSquaredCoherenceXCluster_peak = squeeze(std (CoherenceXSubjectsXCluster_peakEGG(:,:,31).*CoherenceXSubjectsXCluster_peakEGG(:,:,31))./sqrt(length(subjects)));

minSquaredCoherenceXCluster_peak = squeeze(min(CoherenceXSubjectsXCluster_peakEGG(:,:,31).*CoherenceXSubjectsXCluster_peakEGG(:,:,31)));

maxSquaredCoherenceXCluster_peak = squeeze(max(CoherenceXSubjectsXCluster_peakEGG(:,:,31).*CoherenceXSubjectsXCluster_peakEGG(:,:,31)));

SquaredCoherenceXClusterXSubject_peak = squeeze(CoherenceXSubjectsXCluster_peakEGG(:,:,31).*CoherenceXSubjectsXCluster_peakEGG(:,:,31));
meanSquaredCoherenceXSubject_peak = squeeze(mean(CoherenceXSubjectsXCluster_peakEGG(:,2:end,31).*CoherenceXSubjectsXCluster_peakEGG(:,2:end,31),2));% 1 is global signal

%% load mean CS in Si, SIIr, MWM and EBA in all participants
insideBrain = tools_getIndexBrain('inside');


gastricNetwork_overlapPupil = ft_read_mri('Y:\ClusterResults\kw3\CA0050\Cluster_nR10000_CA0050_kw3_fir2_fspread_015_fOrder_5_tw_15csfr_ClusterMap.nii');
gastricNetwork_overlapPupil = gastricNetwork_overlapPupil.anatomy(:);
gastricNetwork_overlapPupil= gastricNetwork_overlapPupil(insideBrain);

nhist(gastricNetwork_overlapPupil(find(gastricNetwork_overlapPupil)))

clustersOI = [4 10 3 9];
clustersOI_matches = ismember(gastricNetwork_overlapPupil, clustersOI);
gastricNetwork_overlapPupil_indexes = find(clustersOI_matches);


subjects = [8 9 10 13 16 17 18 19 20 21 22 23 26 29 32 33 34 35 36 39]


emp_overlap = zeros(length(subjects),length(gastricNetwork_overlapPupil_indexes));
chance_overlap = zeros(length(subjects),length(gastricNetwork_overlapPupil_indexes));


for iSubject=1:length(subjects)
    subj_idx = subjects(iSubject)
    
    filename_PLV= global_filename(subj_idx,cfgMain,'PLVXVoxelFilename_csfr');
    empPLV = ft_read_mri([filename_PLV,'.nii']);
    empPLV= empPLV.anatomy(insideBrain);
    filename_surrogate= global_filename(subj_idx,cfgMain,'medianRotationFilename_csfr');
    chancePLV = ft_read_mri([filename_surrogate,'.nii']);
        chancePLV= chancePLV.anatomy(insideBrain);

    emp_overlap(iSubject,:)= empPLV(gastricNetwork_overlapPupil_indexes);
    chance_overlap(iSubject,:)= chancePLV(gastricNetwork_overlapPupil_indexes);
end


cs_overlap = emp_overlap-chance_overlap;
cs_overlap_mean = mean(cs_overlap,2);


coherence_peak_pupil = CoherenceXSubjectsXCluster_peakEGG(:,31).*CoherenceXSubjectsXCluster_peakEGG(:,31);

figure
plot(coherence_peak_pupil,cs_overlap_mean,'ok')

[r p]=corrcoef(coherence_peak_pupil,cs_overlap_mean)

corrbf(r(3),20)

%{

This scripts compares estimates of shared variance between all nodes of the
gastric network using the squared coherence and squared correlation
coefficient. It then performs a ttest to check for significant differences.
The resulting values are reported in the text of the paper, section 'Temporal advances and lags in the gastric network '

Input:
Cluster timeseries
    Y:\Subjects\Subject13\Timeseries\PhasesAnalysis\ClusterTimeseries_csfr_S_13_kw3_CA0050

Output:
command line output from section %% ttest correlation and coherence


IR 30/06/2017
%}

%% Set parameters

peaksAllsubjects = global_getEGGpeaks;
subjects = global_subjectList;
cfgMain = global_getcfgmain;

% initialize coherence at peak
squaredCoherence_peak = zeros(30,66);

% iterate through subjects
for iSubj = 1:length(subjects)

    
%% 1 timeseries gasnet 
    subj_idx = subjects(iSubj)
    clusterTimeseries_filename = global_filename(subj_idx,cfgMain,'clusterTimeseries_filename');
    load(clusterTimeseries_filename)

timeseries2Coherence = clusterTimeseries(cfgMain.beginCut:cfgMain.endCut,:)'; % cut first and last 15 volumes due to filter artifact
    
%load base fieldtrip data structure
load(strcat(global_path2root,'scripts4paper\files\sampleFieldtripStruc.mat'))

labelsChannelsMAIN = {'vOcc','dOCC','MWM','SSr','pCS','dPrec','RSC','SIIl','EBAr','SIIr','dPRECla','sPOS'};

labelsChannels = labelsChannelsMAIN;
clusterRegionsComparisons = timeseries2Coherence;
dataStructure.hdr = EGG_downsampled.hdr;
dataStructure.fsample = 0.5;
dataStructure.time{1,1}  = [0:2:(size(clusterRegionsComparisons,2)*2)-1];dataStructure.label = labelsChannels;%channelStr;
dataStructure.cfg = EGG_downsampled.cfg;
dataStructure.trial{1,1} = clusterRegionsComparisons;

% Cut data into 120 s long trials with 20 seconds overlap
cfgWelch = [];
cfgWelch.lengthWindow = 120;
cfgWelch.overlap = 6;
len = dataStructure.fsample*cfgWelch.lengthWindow; % length of subtrials cfg.length s in samples
dataStructure.sampleinfo=[1 max(dataStructure.time{1,1})*dataStructure.fsample];
cfg = [];
cfg.trl(:,1) = dataStructure.sampleinfo(1):(len/cfgWelch.overlap):dataStructure.sampleinfo(2)-len+1;%trial start in samples from begining of raw data
cfg.trl(:,2) = dataStructure.sampleinfo(1)+len-1:(len/cfgWelch.overlap):dataStructure.sampleinfo(2);%trial ends in samples from begining of raw data
cfg.trl(:,3) = 0; %offset of the trigger with respect to the trial
data_trials = ft_redefinetrial(cfg,dataStructure);

% computing fourier
cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.output = 'powandcsd';
cfg.pad = 1000;
cfg.foilim = [1/cfgWelch.lengthWindow 0.1]; % 0 - 6 cpm
cfg.keeptrials = 'no';
frequencyWelch = ft_freqanalysis(cfg,data_trials);
cfg            = [];

% Compuite coherence
cfg.method     = 'coh';
coherence      = ft_connectivityanalysis(cfg, frequencyWelch);

% find frequency bin corresponding to EGG peak
indPeak_subj = find (peaksAllsubjects(:,1) == subj_idx);
indPeak = find(coherence.freq == peaksAllsubjects(indPeak_subj,3));

% square it
squaredCoherence_peak(iSubj,:) = coherence.cohspctrm(:,indPeak).*coherence.cohspctrm(:,indPeak);
    
end

% average across clusters
squaredCoherence_peak_mean = mean(squaredCoherence_peak)'

mean(squaredCoherence_peak_mean)
std(squaredCoherence_peak_mean)
[c,i]= min(squaredCoherence_peak_mean)
[c,i]=max(squaredCoherence_peak_mean)

% transform vector to triangular matrix (for visualization only)
indexesMatrix = triu(ones(12),1)
matrixCorr = ones(12);
matrixCorr = matrixCorr(:)
matrixCorr(logical(indexesMatrix(:))) = squaredCoherence_peak_mean;
matrixCorr = reshape(matrixCorr,12,12)
matrixCorr(triu(matrixCorr,1)) = squaredCoherence_peak_mean

figure
imagesc(matrixCorr,[0 0.7])
colormap('hot')
title('coherence')

%% Compute correlation

timeseries_allsubjects_corr_allgasnet = zeros(30,12,450); 
for iSubj = 1:length(subjects)

    
% Load timeseries gasnet 
    subj_idx = subjects(iSubj)
    clusterTimeseries_filename = global_filename(subj_idx,cfgMain,'clusterTimeseries_filename');
    load(clusterTimeseries_filename)
    timeseries_allsubjects_corr_allgasnet(iSubj,:,:) = clusterTimeseries';
   
   
% Filter everything for pearson

% find subjects EGG peak
indPeak = find (peaksAllsubjects(:,1) == subj_idx);
mostPowerfullFrequency = peaksAllsubjects(indPeak,3);

centerFrequency = mostPowerfullFrequency; %
filter_frequency_spread=cfgMain.frequencySpread/1000; % In hz
sr = 0.5; % 1 TR = 2s
filterOrder=(cfgMain.fOrder*fix(sr/(centerFrequency-filter_frequency_spread))-1);%in nsamples
transition_width= cfgMain.transitionWidth/100; % in normalised units
% filteredMRI=tools_bpFilter(squeeze(timeseries_allsubjects_coh(iSubj,:,:))',sr,filterOrder,centerFrequency,filter_frequency_spread,transition_width,cfgMain.filterType);
filteredMRIgasnetALL=tools_bpFilter(squeeze(timeseries_allsubjects_corr_allgasnet(iSubj,:,:))',sr,filterOrder,centerFrequency,filter_frequency_spread,transition_width,cfgMain.filterType);
timeseries_allsubjects_corr_allgasnet(iSubj,:,:) =filteredMRIgasnetALL';
end

% cut 15 first and last volumes due to filter artefact
timeseries_allsubjects_corr_allgasnet = timeseries_allsubjects_corr_allgasnet(:,:,cfgMain.beginCut:cfgMain.endCut);


squaredPearsonAllgasnet= zeros(30,66);

for iSubj = 1:length(subjects)

data4correlationAllgasnet = squeeze(timeseries_allsubjects_corr_allgasnet(iSubj,:,:))';
correlationMatrixSquaredAllgasnet = corrcoef(data4correlationAllgasnet).*corrcoef(data4correlationAllgasnet);
% the correlation matrix is simmetrical, so it will include all the values twice, by taking the
% upper triangle of the matrix (without the identity diagonal line), we can recover only the first occurance of 
% each coherence value
triu_notDiagAllgasnet = triu(correlationMatrixSquaredAllgasnet,1)';
% the resulting triangular matrix is transforemd into vector
% The order in which the resulting values correspond to the order of the 
% output of fieldtrip coherence analysis
triu_notDiag_vectorAllgasnet = triu_notDiagAllgasnet(:);
% zeros correspond to the elements of the lower triangle and the identity
% line
triu_notDiag_vector_nozerosAllgasnet = triu_notDiag_vectorAllgasnet(find(triu_notDiag_vectorAllgasnet));
squaredPearsonAllgasnet(iSubj,:) = triu_notDiag_vector_nozerosAllgasnet;
end

squaredPearsonAllgasnetMatrix= zeros(30,12,12);


%% visualize the values of each node in correlation matrix 
for iSubj = 1:length(subjects)
data4correlationAllgasnet = squeeze(timeseries_allsubjects_corr_allgasnet(iSubj,:,:))';
squaredPearsonAllgasnetMatrix(iSubj,:,:) = corrcoef(data4correlationAllgasnet).*corrcoef(data4correlationAllgasnet);
end
squaredPearsonAllgasnetMatrix_mean = squeeze(mean(squaredPearsonAllgasnetMatrix));
figure
imagesc(triu(squaredPearsonAllgasnetMatrix_mean,1),[0 0.7])
colormap('hot')
title('correlation')


%% Get mean and range of correlation across nodes
mean(mean(squaredPearsonAllgasnet))
std(mean(squaredPearsonAllgasnet))

[c,i]= min(mean(squaredPearsonAllgasnet))
[c,i]=max(mean(squaredPearsonAllgasnet))

%% ttest correlation and coherence
mean_subject_squaredCoherence_peak = mean(squaredCoherence_peak,2)
mean_subject_squaredPearsonAllgasnet = mean(squaredPearsonAllgasnet,2)
[h p ci stats]=ttest(mean_subject_squaredCoherence_peak,mean_subject_squaredPearsonAllgasnet)


mean(squaredCoherence_peak,2)>mean(squaredPearsonAllgasnet,2)
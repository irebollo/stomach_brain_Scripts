%{

This script produces the analysis presented in figure 4D
It compares the estimates of functional connectivity (shared variance)
between gastric, saliency and default regions using coherence and correlations
It thens plots the timeseries of an example subject

Inputs:

Cluster timeseries
    Y:\Subjects\Subject13\Timeseries\PhasesAnalysis\ClusterTimeseries_csfr_S_13_kw3_CA0050
rVAI mask
    Y:\MasksFromOtherPapers\Deen_Insula\RVAI_3mm.nii

Output:
Figures from section %% Plot figures output
used in figure 4 of the paper

IR 20/06/2017
%}

%% Set parameters
peaksAllsubjects = global_getEGGpeaks;
subjects = global_subjectList;
cfgMain = global_getcfgmain;

insideBrain = tools_getIndexBrain('inside');

% right ventral anterior insula ROI
rVAI = ft_read_mri([global_path2root 'MasksFromOtherPapers' filesep 'Deen_Insula' filesep 'RVAI_3mm.nii']);
rVAI = logical(rVAI.anatomy(:)); % 3d to vector
rVAI = rVAI(insideBrain);

% DMN ROI
DMNcenter_tal = [-5 -49 40];% tal Fox 2005 PNAS
DMNcenter = tal2mni(DMNcenter_tal')  % transform to MNI
DMN_center_vox = tools_mni2vox(DMNcenter,1) % transform to voxel coordinates 
emptyBrain = zeros (53,63,46); % empty 3d Volume for storing empirical PLV
% take cube of 3 voxels centered at coodinates
emptyBrain(DMN_center_vox(1)-1:DMN_center_vox(1)+1,DMN_center_vox(2)-1:DMN_center_vox(2)+1,DMN_center_vox(3)-1:DMN_center_vox(3)+1 ) = 1
DMNCoreMask = emptyBrain(insideBrain);
DMNCoreMask = find(DMNCoreMask);


%% Retrieve timeseries 

timeseries_allsubjects_coh = zeros(30,6,450); % 6 timeseries
timeseries_allsubjects_corr = zeros(30,6,450); 
for iSubj = 1:length(subjects)

    
%% 1 timeseries gasnet 
    subj_idx = subjects(iSubj)
    clusterTimeseries_filename = global_filename(subj_idx,cfgMain,'clusterTimeseries_filename');
    load(clusterTimeseries_filename)
    timeseries_allsubjects_coh(iSubj,1:4,:) = clusterTimeseries(:,[1 3 4 9])';
% clusters timeseries being loaded
% SSr EBA RCZp Cun
   %% 1 timeseries insula and pcc

% whole brain timeseries
BOLDTimeseriesFilename = global_filename(subj_idx,cfgMain,strcat('filename_',cfgMain.Timeseries2Regress,'_Residuals_FB'));
timeseries = load(BOLDTimeseriesFilename);

% rVAI
timeseries_allsubjects_coh (iSubj,5,:)= nanmean(timeseries.error_csf_z(:,rVAI),2)';
   
% prec
timeseries_allsubjects_coh (iSubj,6,:)= nanmean(timeseries.error_csf_z(:,DMNCoreMask),2)';


   
   %% Filter everything for pearson
   
indPeak = find (peaksAllsubjects(:,1) == subj_idx);
mostPowerfullFrequency = peaksAllsubjects(indPeak,3);
centerFrequency = mostPowerfullFrequency; %
filter_frequency_spread=cfgMain.frequencySpread/1000; % In hz
sr = 0.5; % 1 TR = 2s
filterOrder=(cfgMain.fOrder*fix(sr/(centerFrequency-filter_frequency_spread))-1);%in nsamples
transition_width= cfgMain.transitionWidth/100; % in normalised units
filteredMRI=tools_bpFilter(squeeze(timeseries_allsubjects_coh(iSubj,:,:))',sr,filterOrder,centerFrequency,filter_frequency_spread,transition_width,cfgMain.filterType);

timeseries_allsubjects_corr(iSubj,:,:) = filteredMRI';
end

% cut first and last 15 volumes to get rid of edge artifact of the filter
timeseries_allsubjects_coh = timeseries_allsubjects_coh(:,:,cfgMain.beginCut:cfgMain.endCut);
timeseries_allsubjects_corr = timeseries_allsubjects_corr(:,:,cfgMain.beginCut:cfgMain.endCut);

%% get coherence squared
squaredCoherence_peak = zeros(30,15);
for iSubj = 1:length(subjects)

timeseries2Coherence = squeeze(timeseries_allsubjects_coh(iSubj,:,:)); % get rid of singleton_dimension
  

nVoxels = size(timeseries2Coherence,1); % number of timeseries

%load base fieldtrip data structure
load([global_path2root,'scripts4paper' filesep 'files' filesep 'sampleFieldtripStruc.mat'])

labelsChannelsMAIN = {'dOcc','rCZP','rSSI','EBA','vAI','Prec'};
labelsChannels = labelsChannelsMAIN;
clusterRegionsComparisons = timeseries2Coherence;
dataStructure.hdr = EGG_downsampled.hdr;
dataStructure.fsample = 0.5; % hz
dataStructure.time{1,1}  = [0:2:(size(clusterRegionsComparisons,2)*2)-1];dataStructure.label = labelsChannels;%channelStr;
dataStructure.cfg = EGG_downsampled.cfg;
dataStructure.trial{1,1} = clusterRegionsComparisons;

% computing fourier

% cut data into trials
cfgWelch = [];
cfgWelch.lengthWindow = 120; % in seconds
cfgWelch.overlap = 6; % 1/6 i.e. 20 seconds
len = dataStructure.fsample*cfgWelch.lengthWindow; % length of subtrials cfg.length s in samples
dataStructure.sampleinfo=[1 max(dataStructure.time{1,1})*dataStructure.fsample];

cfg = [];
cfg.trl(:,1) = dataStructure.sampleinfo(1):(len/cfgWelch.overlap):dataStructure.sampleinfo(2)-len+1;%trial start in samples from begining of raw data
cfg.trl(:,2) = dataStructure.sampleinfo(1)+len-1:(len/cfgWelch.overlap):dataStructure.sampleinfo(2);%trial ends in samples from begining of raw data
cfg.trl(:,3) = 0; %offset of the trigger with respect to the trial
data_trials = ft_redefinetrial(cfg,dataStructure); % 

% compute spectrum using welch

cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.output = 'powandcsd';
cfg.pad = 1000;
cfg.foilim = [1/cfgWelch.lengthWindow 0.1]; 
cfg.keeptrials = 'no';
frequencyWelch = ft_freqanalysis(cfg,data_trials);

% compute coherence
cfg            = [];
cfg.method     = 'coh';
coherence      = ft_connectivityanalysis(cfg, frequencyWelch);


% get coherence at the frequency bin corresponding to the peak fo the EGG
indPeak_subj = find (peaksAllsubjects(:,1) == subj_idx);
indPeak = find(coherence.freq == peaksAllsubjects(indPeak_subj,3));
squaredCoherence_peak(iSubj,:) = coherence.cohspctrm(:,indPeak).*coherence.cohspctrm(:,indPeak); % square it to get from coherence to sharted variance
    
end



%% get correlation squared
squaredPearson= zeros(30,15); % Initialize
squaredPearsonAllgasnet= zeros(30,66);

for iSubj = 1:length(subjects)

data4correlation = squeeze(timeseries_allsubjects_corr(iSubj,:,:))'; % remove the first singleton dimension

correlationMatrixSquared = corrcoef(data4correlation).*corrcoef(data4correlation); % get the correlation matrix of all timeseries

% the correlation matrix is simmetrical, so it will include all the values twice, by taking the
% upper triangle of the matrix (without the identity diagonal line), we can recover only the first occurance of 
% each coherence value
triu_notDiag = triu(correlationMatrixSquared,1)';
% the resulting triangular matrix is transforemd into vector
% The order in which the resulting values correspond to the order of the 
% output of fieldtrip coherence analysis
triu_notDiag_vector = triu_notDiag(:);
% zeros correspond to the elements of the lower triangle and the identity
% line
triu_notDiag_vector_nozeros = triu_notDiag_vector(find(triu_notDiag_vector));

squaredPearson(iSubj,:) = triu_notDiag_vector_nozeros;

end
%% get mean and SER across subjects
squaredPearson_mean = mean(squaredPearson);
squaredPearson_SER = std(squaredPearson)./sqrt(30);

squaredCoherence_peak_mean = mean(squaredCoherence_peak);
squaredCoherence_peak_SER = std(squaredCoherence_peak)./sqrt(30);

%% Plot figures output

data2plot =[squaredPearson_mean([6 7 8 9]),squaredCoherence_peak_mean([6 7 8 9])]% RCZP-rSS,RCZP-OT, RCZP-vAI RCZP-DMN
data2plot_error = [squaredPearson_SER([6 7 8 9]),squaredCoherence_peak_SER([6 7 8 9])]% RCZP-rSS,RCZP-OT, RCZP-vAI RCZP-DMN
figure
bar(data2plot,'grouped')
hold on
errorbar(data2plot(:),data2plot_error(:),'.k','LineWidth',2)



data2plot =[squaredPearson_mean([2 3 4 5]),squaredCoherence_peak_mean([2 3 4 5])]% vOCC-rSS,vOCC-OT, vOCC-vAI vOCC-DMN
data2plot_error = [squaredPearson_SER([2 3 4 5]),squaredCoherence_peak_SER([2 3 4 5])]% vOCC-rSS,vOCC-OT, vOCC-vAI vOCC-DMN
figure
bar(data2plot,'grouped')
hold on
errorbar(data2plot(:),data2plot_error(:),'.k','LineWidth',2)

%% Boxplots
figure
boxplot([squaredPearson(:,[6 7 8 9]) , squaredCoherence_peak(:,[6 7 8 9]) ,squaredPearson(:,[2 3 4 5]), squaredCoherence_peak(:,[2 3 4 5]) ])



%% plot representative subject timeseries

% cut data of all subjects into trials

data_trials_allsubjects = zeros(30,7,3,100);
for iSubj = 1:length(subjects)
    timeseries2trials = squeeze(timeseries_allsubjects_corr(iSubj,[2 4 6],:));
    

nVoxels = size(timeseries2trials,1);
%load base fieldtrip data structure
load(strcat(global_path2root,'scripts4paper\files\sampleFieldtripStruc.mat'))

labelsChannelsMAIN = {'rCZP','EBA','Prec'};
labelsChannels = labelsChannelsMAIN;
clusterRegionsComparisons = timeseries2trials;
dataStructure.hdr = EGG_downsampled.hdr;
dataStructure.fsample = 0.5;
dataStructure.time{1,1}  = [0:2:(size(clusterRegionsComparisons,2)*2)-1];dataStructure.label = labelsChannels;%channelStr;
dataStructure.cfg = EGG_downsampled.cfg;
dataStructure.trial{1,1} = clusterRegionsComparisons;


cfgWelch = [];
cfgWelch.lengthWindow = 200;
cfgWelch.overlap = 2;
len = dataStructure.fsample*cfgWelch.lengthWindow; % length of subtrials cfg.length s in samples
dataStructure.sampleinfo=[1 max(dataStructure.time{1,1})*dataStructure.fsample];
cfg = [];
cfg.trl(:,1) = dataStructure.sampleinfo(1):(len/cfgWelch.overlap):dataStructure.sampleinfo(2)-len+1;%trial start in samples from begining of raw data
cfg.trl(:,2) = dataStructure.sampleinfo(1)+len-1:(len/cfgWelch.overlap):dataStructure.sampleinfo(2);%trial ends in samples from begining of raw data
cfg.trl(:,3) = 0; %offset of the trigger with respect to the trial
data_trials = ft_redefinetrial(cfg,dataStructure);

for iTrial = 1:7
data_trials_allsubjects(iSubj,iTrial,:,:)= data_trials.trial{iTrial};
end
end

crosscorr_allsubjectes_trials = zeros(30,7,3,5);
for iSubj = 1:length(subjects)
for iTrial = 1:7

[crosscorr_allsubjectes_trials(iSubj,iTrial,1,:),lags]=xcorr(squeeze(squeeze(data_trials_allsubjects(iSubj,iTrial,1,:))),squeeze(squeeze(data_trials_allsubjects(iSubj,iTrial,2,:))),[2],'coeff');
[crosscorr_allsubjectes_trials(iSubj,iTrial,2,:),lags]=xcorr(squeeze(squeeze(data_trials_allsubjects(iSubj,iTrial,2,:))),squeeze(squeeze(data_trials_allsubjects(iSubj,iTrial,3,:))),[2],'coeff');
[crosscorr_allsubjectes_trials(iSubj,iTrial,3,:),lags]=xcorr(squeeze(squeeze(data_trials_allsubjects(iSubj,iTrial,3,:))),squeeze(squeeze(data_trials_allsubjects(iSubj,iTrial,1,:))),[2],'coeff');
end
end


comparisonsOfInterest = [7 8 9]

diffDoCCSS = ...
squaredPearson(:,comparisonsOfInterest) - squaredCoherence_peak(:,comparisonsOfInterest)
figure;pcolor(diffDoCCSS)


iSubj = 30
time2plot = [1:2:200]; % this one works quite well

data2plot = squeeze(squeeze(data_trials_allsubjects(30,5,:,:)));

% time2plot = [2:2:840];
figure
plot(time2plot,data2plot(3,:),'color',[0.7 0.7 0.7],'linewidth',3)
hold on
plot(time2plot,data2plot(1,:),'color',[1 0 0],'linewidth',3)
plot(time2plot,data2plot(2,:),'color',[0.7 0 0],'linewidth',3)
legend('Prec','rCZP','OT')
ylabel('Z amplitude')
xlabel('Seconds')


% corrcoef(squeeze(timeseries_allsubjects_corr(iSubj,1,time2plot/2)),squeeze(timeseries_allsubjects_corr(iSubj,3,time2plot/2)))
[r,lags]=xcorr(squeeze(timeseries_allsubjects_corr(iSubj,2,time2plot/2)),squeeze(timeseries_allsubjects_corr(iSubj,4,time2plot/2)),[3],'coeff')
[rOTinsula,lags]=xcorr(squeeze(timeseries_allsubjects_corr(iSubj,2,time2plot/2)),squeeze(timeseries_allsubjects_corr(iSubj,6,time2plot/2)),[3],'coeff')
[rrCZPinsula,lags]=xcorr(squeeze(timeseries_allsubjects_corr(iSubj,4,time2plot/2)),squeeze(timeseries_allsubjects_corr(iSubj,6,time2plot/2)),[3],'coeff')
figure
hold on
plot(lags*2,(r.*r)*100,'r','linewidth',3)
plot(lags*2,(rOTinsula.*rOTinsula)*100,'color',[0.7 0.7 0.7],'linewidth',3)
plot(lags*2,(rrCZPinsula.*rrCZPinsula)*100,'color',[0.5 0.5 0.5],'linewidth',3)
ylabel('%Shared Variance')
grid on
xlabel('time lag in seconds')



%% Get coherence in example subject
timeseries2Coherence = squeeze(squeeze(data_trials_allsubjects(30,5,:,:)));
subj_idx = subjects(30)    

nVoxels = size(timeseries2Coherence,1);
%load base fieldtrip data structure
load(strcat(global_path2root,'scripts4paper\files\sampleFieldtripStruc.mat'))

labelsChannelsMAIN = {'rCZP','OT','vPrec'};
labelsChannels = labelsChannelsMAIN;
clusterRegionsComparisons = timeseries2Coherence;
dataStructure.hdr = EGG_downsampled.hdr;
dataStructure.fsample = 0.5;
dataStructure.time{1,1}  = [0:2:(size(clusterRegionsComparisons,2)*2)-1];dataStructure.label = labelsChannels;%channelStr;
dataStructure.cfg = EGG_downsampled.cfg;
dataStructure.trial{1,1} = clusterRegionsComparisons;

% computing fourier
cfgWelch = [];
cfgWelch.lengthWindow = 60;
cfgWelch.overlap = 2;
len = dataStructure.fsample*cfgWelch.lengthWindow; % length of subtrials cfg.length s in samples
dataStructure.sampleinfo=[1 max(dataStructure.time{1,1})*dataStructure.fsample];
cfg = [];
cfg.trl(:,1) = dataStructure.sampleinfo(1):(len/cfgWelch.overlap):dataStructure.sampleinfo(2)-len+1;%trial start in samples from begining of raw data
cfg.trl(:,2) = dataStructure.sampleinfo(1)+len-1:(len/cfgWelch.overlap):dataStructure.sampleinfo(2);%trial ends in samples from begining of raw data
cfg.trl(:,3) = 0; %offset of the trigger with respect to the trial
data_trials = ft_redefinetrial(cfg,dataStructure);
cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.output = 'powandcsd';
cfg.pad = 1000;
cfg.foilim = [1/cfgWelch.lengthWindow 0.1]; % 0 - 6 cpm
cfg.keeptrials = 'no';
frequencyWelch = ft_freqanalysis(cfg,data_trials);
cfg            = [];
cfg.method     = 'coh';
coherence      = ft_connectivityanalysis(cfg, frequencyWelch);


indPeak_subj = find (peaksAllsubjects(:,1) == subj_idx);
indPeak = find(coherence.freq == peaksAllsubjects(indPeak_subj,3));
squaredCoherence_peak_exampleSubject = coherence.cohspctrm(:,indPeak).*coherence.cohspctrm(:,indPeak);

data_crosscor = squeeze(squeeze(crosscorr_allsubjectes_trials(30,5,:,:)));
    
figure
bar([data_crosscor(1,3).*data_crosscor(1,3);data_crosscor(1,2).*data_crosscor(1,2);squaredCoherence_peak_exampleSubject(1) ])


plot(lags,(data_crosscor(1,:).*data_crosscor(1,:))*100,'r','linewidth',3)

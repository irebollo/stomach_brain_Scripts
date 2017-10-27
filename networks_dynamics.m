paper_network_dynamics
% Perform the EGG-BOLD dPLV analysis

%% CFG giles
cfgMain = global_getcfgmain;
subjects = global_subjectList;
labelsChannelsMAIN = {'vOcc','dOCC','MWM','SSr','pCS','dPrec','RSC','SIIl','EBAr','SIIr','dPRECla','sPOS'}; % Cluster Names


%% Load cluster timeseries of each participant

clusTS_alls = zeros(30,12,420); % cluster timeseries for all subjects
EnvClusTS_alls = zeros(30,12,420); % amplitude envelope of cluster timeseries

clusTSphases_alls = zeros(30,12,420);% Phases of all clusters
EGGPhases_alls =zeros(30,420);% Phase of EGG

dPLV_allClus_alls = zeros(30,12,78); % dinamic PLV with EGG all clusters


EnvClusTS_alls_trials = zeros(30,12,78); % average amplitude envelope of cluster timeseries per trial

%% Looop
for iSubj = 1:length(subjects)

% Load
subj_idx = subjects(iSubj)

clusterTimeseries_filename = global_filename(subj_idx,cfgMain,'clusterTimeseries_filename');
load(clusterTimeseries_filename)
EGGPhaseXVolumeFilename = global_filename(subj_idx,cfgMain,'EGGPhaseXVolumeFilename');
load(EGGPhaseXVolumeFilename)

filename_EGG_amplitude = global_filename(subj_idx,cfgMain,'FilenameamplitudeXVolumeBestChannel_FULLBAND');
load(filename_EGG_amplitude)

% EGG peak info
peaksAllsubjects = global_getEGGpeaks;
indPeak = find (peaksAllsubjects(:,1) == subj_idx);
mostPowerfullFrequency = peaksAllsubjects(indPeak,3);

% Filter BOLD and EGG
centerFrequency = mostPowerfullFrequency; %
filter_frequency_spread=cfgMain.frequencySpread/1000; % In hz
sr = 0.5; % 1 TR = 2s
filterOrder=(cfgMain.fOrder*fix(sr/(centerFrequency-filter_frequency_spread))-1);%in nsamples
transition_width= cfgMain.transitionWidth/100; % in normalised units
filteredClusters=tools_bpFilter(clusterTimeseries,sr,filterOrder,centerFrequency,filter_frequency_spread,transition_width,cfgMain.filterType);
filteredEGG = tools_bpFilter(EGG_FullBand,sr,filterOrder,centerFrequency,filter_frequency_spread,transition_width,cfgMain.filterType);

% Hilbert
phaseClusters = hilbert(filteredClusters(cfgMain.beginCut:cfgMain.endCut,:));


% Cut data into trials for dPLV
timeseries2Coherence = [angle(phaseXVolume);angle(phaseClusters)';filteredEGG(cfgMain.beginCut:cfgMain.endCut);filteredClusters(cfgMain.beginCut:cfgMain.endCut,:)';abs(phaseClusters)']; 
% 1, phase EGG, 2-13 phase clusters, 14 amp/volume EGG 15:26 amp clusters,
% 27-38 amplitude envelop of clusters, 39 Order timeseries

nVoxels = size(timeseries2Coherence,1);

%load base fieldtrip data structure
load(strcat(global_path2root,'scripts4paper\files\sampleFieldtripStruc.mat'))


% Fieldtrip require every voxel to have a channnel name
channelStr=cell(nVoxels,1);
for iVoxel = 1:nVoxels
    channelList(iVoxel,1) = iVoxel;
    channelStr(iVoxel) = cellstr(mat2str(iVoxel)); 
end

labelsChannels = channelStr;
clusterRegionsComparisons = timeseries2Coherence;
dataStructure.hdr = EGG_downsampled.hdr;
dataStructure.fsample = 0.5;
dataStructure.time{1,1}  = [0:2:(size(clusterRegionsComparisons,2)*2)-1];dataStructure.label = labelsChannels;%channelStr;
dataStructure.cfg = EGG_downsampled.cfg;
dataStructure.trial{1,1} = clusterRegionsComparisons;


% 60 seconds 10 seconds overlap
% Cut data into trials
cfgWelch = [];
cfgWelch.lengthWindow = 60; %seconds
cfgWelch.overlap = 6;% propotion i.e 1/4 i.e. 10 seconds

len = dataStructure.fsample*cfgWelch.lengthWindow; % length of subtrials cfg.length s in samples
dataStructure.sampleinfo=[1 max(dataStructure.time{1,1})*dataStructure.fsample];

cfg = [];
cfg.trl(:,1) = dataStructure.sampleinfo(1):(len/cfgWelch.overlap):dataStructure.sampleinfo(2)-len+1;%trial start in samples from begining of raw data
cfg.trl(:,2) = dataStructure.sampleinfo(1)+len-1:(len/cfgWelch.overlap):dataStructure.sampleinfo(2);%trial ends in samples from begining of raw data
cfg.trl(:,3) = 0; %offset of the trigger with respect to the trial
data_trials = ft_redefinetrial(cfg,dataStructure);


% Compute dPLV
dPLV_allclusters_60_10 = zeros(12,length(data_trials.trial));
for iStep = 1:length(data_trials.trial)
dPLV_allclusters_60_10(:,iStep) = abs (mean (exp (1i* (bsxfun (@minus , data_trials.trial{1,iStep}(2:13,:)', data_trials.trial{1,iStep}(1,:)')))))'; % get PLV
end


for iStep=1:length(data_trials.trial)
time2plot(iStep) = cfgWelch.lengthWindow/2+ (iStep-1)*10
end


amp_trial = zeros(12,length(data_trials.trial));
for iStep = 1:length(data_trials.trial)
for iCluster = 1:12
amp_trial(iCluster,iStep)= mean(data_trials.trial{1,iStep}(26+iCluster,:));
end %iCluster
end


% save data
clusTS_alls(iSubj,:,:) = filteredClusters(cfgMain.beginCut:cfgMain.endCut,:)';
clusTSphases_alls(iSubj,:,:) = angle(phaseClusters)';
EGGPhases_alls(iSubj,:) = angle(phaseXVolume);
dPLV_allClus_alls(iSubj,:,:) = dPLV_allclusters_60_10;
EnvClusTS_alls(iSubj,:,:) = abs(phaseClusters)';
EnvClusTS_alls_trials(iSubj,:,:) = amp_trial;

end


%% Get corr dPLV bold Amplitude and perform stats
r_dPLV_Amp = zeros(30,12);

for iSubj =1:30
for iCluster = 1:12

[r,p]= corrcoef(squeeze(dPLV_allClus_alls(iSubj,iCluster,:)),...
    squeeze(EnvClusTS_alls_trials(iSubj,iCluster,:)));

r_dPLV_Amp(iSubj,iCluster) = r(3)

end
end


z_dPLV_Amp=fisherz(r_dPLV_Amp)
z_dPLV_Amp = reshape(z_dPLV_Amp,30,12)
[h p ci stats ]= ttest(z_dPLV_Amp)

[pthr pcor padj] = fdr(p)

tablep=[p', (p*12)', padj']



figure
bar(mean(r_dPLV_Amp))
hold on
errorbar(mean(r_dPLV_Amp),(std(r_dPLV_Amp)./sqrt(30)),'or')
ylabel('r')

set(gca,'XtickLabel',labelsChannelsMAIN)
set(gca,'TickLength',[0 0],'FontSize',12);
% shg
% end

figure
boxplot(r_dPLV_Amp,'colors','k','boxstyle','filled','plotstyle','compact','labels',labelsChannelsMAIN);shg
hold on
plot(mean(r_dPLV_Amp),'+w','linewidth',3)
plot(mean(r_dPLV_Amp),'+k','linewidth',1)

title('r (dPLV - Amplitude envelope) across subjects')
set(findall(gcf,'-property','FontSize'),'FontSize',16)

mean(mean(r_dPLV_Amp))
std(mean(r_dPLV_Amp))



%% Which regiosn fluctuate together with the EGG
for iSubj = 1:30
corrMatrix_alls(iSubj,:,:) = corr(squeeze(dPLV_allClus_alls(iSubj,:,:))');
end

corrMatrix_alls_mean = squeeze(mean(corrMatrix_alls))
mean(corrMatrix_alls(:))
figure;imagesc(triu(corrMatrix_alls_mean,1),[0 0.5]);
colormap('hot')
sum(diag(corrMatrix_alls_mean));
corrMatrix_alls_mean_nodiag = corrMatrix_alls_mean

corrMatrix_alls_mean_nodiag(find(corrMatrix_alls_mean_nodiag==1)) = nan
sum((corrMatrix_alls_mean))
mean((corrMatrix_alls_mean))

upper = triu(corrMatrix_alls_mean,1)
[b,i]=sort(upper(:),'descend')
tiedrank(upper(:))

reshape(tiedrank(upper(:)),12,12)


[h p ci stats] = ttest(corrMatrix_alls)

p_relevant = p(find(triu(squeeze(p),1)))

[pthr pcor padj] = fdr(p_relevant)

p_bonfe = squeeze(p.*66)

https://brainder.org/2011/09/05/fdr-corrected-fdr-adjusted-p-values/


%% Gastric network coupling global

global_dPLV_allClus_alls = zeros(30,1);
for iSubj = 1:30
corrmatrix_isubj = squeeze(corrMatrix_alls(iSubj,:,:));
corrmatrix_isubj = triu(corrmatrix_isubj,1);
global_dPLV_allClus_alls(iSubj,:) = mean(corrmatrix_isubj(:))
end

mean(global_dPLV_allClus_alls)
std(global_dPLV_allClus_alls)
global_dPLV_allClus_alls_fisher = fisherz(global_dPLV_allClus_alls)
[h p ci stats]=ttest(global_dPLV_allClus_alls_fisher)

figure
nhist(global_dPLV_allClus_alls)

% control 
mean(triu(corrMatrix_alls,1),2)

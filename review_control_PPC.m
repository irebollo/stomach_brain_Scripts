% pairwisePhaseconsistency

subjects = global_subjectList;
cfgMain = global_getcfgmain;
insideBrain = tools_getIndexBrain('inside');
indNoBrain = tools_getIndexBrain('outside');

peaksAllsubjects = global_getEGGpeaks;
rootDir= strcat(global_path2root,'\Subjects\');


for iSubject = 1:30

   subj_idx = subjects(iSubject)


% subjects which are not this subject

% surrogatePLVmatrix = zeros(29,length(insideBrain));

%% Load BOLD
% load BOLD of subject unfiltered
filename_bold_input = global_filename(subj_idx,cfgMain,'filename_csfr_Residuals_FB');
load(filename_bold_input)


% subjListForAnalysis = subjects~= subj_idx
% SurrogateSubjectList = subjects(subjListForAnalysis)

% for iSurrogateSubject = 1:29
    
%     tic
    
% currentSurrogateSubject = SurrogateSubjectList(iSurrogateSubject)

% load EGG first other subject
EGGAmplitudeXVolumeFilename = global_filename(subj_idx,cfgMain,strcat('EGGAmplitudeXVolumeFilename'));
load(EGGAmplitudeXVolumeFilename)

% Frequency info of surrogate EGG
% indPeak = find (peaksAllsubjects(:,1) == iSubject);
mostPowerfullFrequency = peaksAllsubjects(iSubject,3);


%Filter
centerFrequency = mostPowerfullFrequency; %
filter_frequency_spread=cfgMain.frequencySpread/1000; % In hz
sr = 0.5; % 1 TR = 2s
filterOrder=(cfgMain.fOrder*fix(sr/(centerFrequency-filter_frequency_spread))-1);%in nsamples
transition_width= cfgMain.transitionWidth/100; % in normalised units
filteredMRI=tools_bpFilter(error_csf_z,sr,filterOrder,centerFrequency,filter_frequency_spread,transition_width,cfgMain.filterType);
% phaseMRI = hilbert(filteredMRI);
% phaseMRI = angle(phaseMRI(cfgMain.beginCut:cfgMain.endCut,:)); % Cut data to have the same length as EGG (cut this way to get rid of fmri edge artifact on EGG)
% error_csf_z

% clear filteredMRI

% filteredMRI = filteredMRI(cfgMain.beginCut:cfgMain.endCut,:);

% Compute sPLV for surrogate
% ONLY IN VOXELS INSIDE BRAIN
% empPLV = zeros (1,length(insideBrain))'; % empty 3d Volume for storing empirical PLV


% empPLV = abs (mean (exp (1i* (bsxfun (@minus , phaseMRI, angle (phaseXVolume)'))))); % get PLV


% bsxfun applies the operation in @minus from the vector of the third input
% to each column of the matrix of the second input

% clear other variables
clear  phaseMRI

%%


% BOLDtimeseries.trialVector = BOLDtimeseries.trialVector(cfgMain.beginCut:cfgMain.endCut,insideBrain);
% BOLDtimeseries.BOLDtsNOPOLI = ft_preproc_polyremoval (BOLDtimeseries.trialVector', 2);
% BOLD_zscored = ft_preproc_standardize (BOLDtimeseries.BOLDtsNOPOLI);



data = [EGGTimeseries;filteredMRI(cfgMain.beginCut:cfgMain.endCut,:)'];


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
cfgWelch.keeptrials = 'yes';
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
cfg.output = 'fourier';
cfg.pad = 1000;
cfg.foilim = [mostPowerfullFrequency-0.001, mostPowerfullFrequency+0.001]; % 0 - 6 cpm
% cfg.foilim = [1/cfgWelch.lengthWindow 0.25]; % 0 - 6 cpm

cfg.keeptrials = cfgWelch.keeptrials;

cfg.channelcmb=cell(71759,2);
for iVoxel = 1:71759

   cfg.channelcmb{iVoxel,1} = '1';
   cfg.channelcmb{iVoxel,2} = char(channelStr(iVoxel+1));
end

frequencyWelch = ft_freqanalysis(cfg,data_trials);



cfg            = [];
cfg.channelcmb=cell(71759,2);
% cfg.channelcmb=cell(717,2);

for iVoxel = 1:71759
%         for iVoxel = 1:717

   cfg.channelcmb{iVoxel,1} = '1';
   cfg.channelcmb{iVoxel,2} = char(channelStr(iVoxel+1));
end
cfg.method     = 'ppc';
coherence             = ft_connectivityanalysis(cfg, frequencyWelch);



figure
nhist(coherence.ppcspctrm(:,coherence.freq==mostPowerfullFrequency))



emptyBrain = zeros (53,63,46); % empty 3d Volume for storing empirical PLV
emptyBrain = emptyBrain(:);
emptyBrain(insideBrain) = mean(coherence.ppcspctrm,2) ;
emptyBrain = reshape(emptyBrain,53,63,46);
PPC_filename =  global_filename(subj_idx,cfgMain,'PPC_filename');

tools_writeMri(emptyBrain,PPC_filename)


end

%% Create surrogate PPC by inversing the EGG
for iSubject = 1:30

   subj_idx = subjects(iSubject)


% subjects which are not this subject

% surrogatePLVmatrix = zeros(29,length(insideBrain));

%% Load BOLD
% load BOLD of subject unfiltered
filename_bold_input = global_filename(subj_idx,cfgMain,'filename_csfr_Residuals_FB');
load(filename_bold_input)


% load EGG 
EGGAmplitudeXVolumeFilename = global_filename(subj_idx,cfgMain,strcat('EGGAmplitudeXVolumeFilename'));
load(EGGAmplitudeXVolumeFilename)

EGGTimeseries = flip(EGGTimeseries)
% figure;plot(EGGTimeseries)
% hold on
% plot(EGGTimeseries,'--b')


% Frequency info of surrogate EGG
% indPeak = find (peaksAllsubjects(:,1) == iSubject);
mostPowerfullFrequency = peaksAllsubjects(iSubject,3);


%Filter
centerFrequency = mostPowerfullFrequency; %
filter_frequency_spread=cfgMain.frequencySpread/1000; % In hz
sr = 0.5; % 1 TR = 2s
filterOrder=(cfgMain.fOrder*fix(sr/(centerFrequency-filter_frequency_spread))-1);%in nsamples
transition_width= cfgMain.transitionWidth/100; % in normalised units
filteredMRI=tools_bpFilter(error_csf_z,sr,filterOrder,centerFrequency,filter_frequency_spread,transition_width,cfgMain.filterType);

clear  phaseMRI

%%


data = [EGGTimeseries;filteredMRI(cfgMain.beginCut:cfgMain.endCut,:)'];

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
cfgWelch.keeptrials = 'yes';
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
cfg.output = 'fourier';
cfg.pad = 1000;
cfg.foilim = [mostPowerfullFrequency-0.001, mostPowerfullFrequency+0.001]; % 0 - 6 cpm
% cfg.foilim = [1/cfgWelch.lengthWindow 0.25]; % 0 - 6 cpm

cfg.keeptrials = cfgWelch.keeptrials;

cfg.channelcmb=cell(71759,2);
for iVoxel = 1:71759

   cfg.channelcmb{iVoxel,1} = '1';
   cfg.channelcmb{iVoxel,2} = char(channelStr(iVoxel+1));
end

frequencyWelch = ft_freqanalysis(cfg,data_trials);



cfg            = [];
cfg.channelcmb=cell(71759,2);
% cfg.channelcmb=cell(717,2);

for iVoxel = 1:71759
%         for iVoxel = 1:717

   cfg.channelcmb{iVoxel,1} = '1';
   cfg.channelcmb{iVoxel,2} = char(channelStr(iVoxel+1));
end
cfg.method     = 'ppc';
coherence             = ft_connectivityanalysis(cfg, frequencyWelch);


figure
nhist(coherence.ppcspctrm(:,coherence.freq==mostPowerfullFrequency))


emptyBrain = zeros (53,63,46); % empty 3d Volume for storing empirical PLV
emptyBrain = emptyBrain(:);
emptyBrain(insideBrain) = mean(coherence.ppcspctrm,2) ;
emptyBrain = reshape(emptyBrain,53,63,46);
PPC_filename =  global_filename(subj_idx,cfgMain,'PPC_surrogate_filename');

tools_writeMri(emptyBrain,PPC_filename)


end

%% do stats
cfgMain.numberofrandomizations = 1000
cfgMain.clusterAlpha=0.005
timeseries_statsCluster_Regression_PPC(subjects,cfgMain)

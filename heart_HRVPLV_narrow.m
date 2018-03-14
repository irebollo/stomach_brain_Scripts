% heart_HRVPLV


% review_control_EGG_othersubjects



%% Parameters
subjects = [8 9 10 13 15 16 17 18 19 20 21 22 23 24 25 26 29 31 32 33 34 35 36 39 40 41 43]
cfgMain = global_getcfgmain;
insideBrain = tools_getIndexBrain('inside');
indNoBrain = tools_getIndexBrain('outside');

load([global_path2root filesep 'scripts4paper' filesep 'files' filesep 'HRV_peaks_info'])
rootDir= strcat(global_path2root,'\Subjects\');


frequencySpread =0.02

for iSubject = 1:length(subjects)

   subj_idx = subjects(iSubject)
   
    tic
    

%% load and preprocess HRV_ts

centerFrequency = centerFrequencyHRVAllSubjects(iSubject)
plotDir = strcat (global_path2subject(subj_idx),'PreprocessingLog',filesep);



HRV_timeseries_filename = global_filename(subj_idx,cfgMain,'HRV_timeseries')
load(HRV_timeseries_filename)


HRV_notdownsampled = ibi_int';
load(strcat(global_path2root,'\scripts4paper\files\sampleFieldtripStruc.mat'))
labelsChannelsMAIN = {'HRV'};
labelsChannels = labelsChannelsMAIN;
clusterRegionsComparisons = HRV_notdownsampled;
dataStructure.hdr = EGG_downsampled.hdr;
dataStructure.fsample = 1;
dataStructure.time{1,1}  = [0:1/dataStructure.fsample:(size(clusterRegionsComparisons,2))-1];
dataStructure.label = labelsChannels;%channelStr;
dataStructure.cfg = EGG_downsampled.cfg;
dataStructure.trial{1,1} = clusterRegionsComparisons;
dataStructure.sampleinfo = [1 length(HRV_notdownsampled)]

disp('Resampling...')
cfg = [];  %initialize configuration structure
cfg.detrend = 'no'; % remove linear trend from the data (done per trial)
cfg.demean = 'yes';
cfg.resamplefs= 0.5; % 4 x top-freq (15 cpm = 0.25 Hz) - Nyquist = 30 cpm  frequency at which the data will be resampled
HRV_downsampled = ft_resampledata(cfg,dataStructure); % This procedure also lowpass filter the data at half the new sr 


HRV_downsampled.trial{1,1} = ft_preproc_polyremoval (HRV_downsampled. trial{1,1}, 2);
HRV_downsampled.trial{1,1} = ft_preproc_standardize (HRV_downsampled. trial{1,1});



filter_frequency_spread=frequencySpread; % In hz Half widfth or fullwidth????
sr = 0.5; % 1 TR = 2s
filterOrder=(cfgMain.fOrder*fix(sr/(centerFrequency-filter_frequency_spread))-1);%in nsamples
transition_width= cfgMain.transitionWidth/100; % in normalised units
filteredHeart=tools_bpFilter(HRV_downsampled.trial{1,1},sr,filterOrder,centerFrequency,filter_frequency_spread,transition_width,cfgMain.filterType);
phaseHeart = hilbert(filteredHeart);

phaseHeart = angle(phaseHeart(cfgMain.beginCut:cfgMain.endCut)); % Cut data to have the same length as EGG (cut this way to get rid of fmri edge artifact on EGG)

HRV_ts_plot = figure
subplot(2,1,1)
plot(filteredHeart)
title(['S',sprintf('%.2d',subj_idx),32,'filtered HRV narrow'],'fontsize',18)

subplot(2,1,2)
plot(phaseHeart)

xlabel('volumes')
plotFilename = strcat(plotDir,'S_',sprintf('%.2d',subj_idx),'_Heart_HRVts_Narrow');
set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gcf, 'PaperPositionMode', 'auto');
print ('-dpng', '-painters', eval('plotFilename'))
print ('-depsc2', '-painters', eval('plotFilename'))
saveas(HRV_ts_plot,strcat(plotFilename,'.fig'))  
%% BOLD

% Load BOLD
% load BOLD of subject unfiltered
filename_bold_input = global_filename(subj_idx,cfgMain,'filename_csfr_Residuals_FB');
load(filename_bold_input)

%Filter

filter_frequency_spread=frequencySpread; % In hz Half widfth or fullwidth????
sr = 0.5; % 1 TR = 2s
filterOrder=(cfgMain.fOrder*fix(sr/(centerFrequency-filter_frequency_spread))-1);%in nsamples
transition_width= cfgMain.transitionWidth/100; % in normalised units
filteredMRI=tools_bpFilter(error_csf_z,sr,filterOrder,centerFrequency,filter_frequency_spread,transition_width,cfgMain.filterType);
phaseMRI = hilbert(filteredMRI);
phaseMRI = angle(phaseMRI(cfgMain.beginCut:cfgMain.endCut,:)); % Cut data to have the same length as EGG (cut this way to get rid of fmri edge artifact on EGG)
% error_csf_z

% check filter output
meanBOLDtsPlot =  figure
plot(nanmean(filteredMRI,2))
xlabel('volumes')
title(['S',sprintf('%.2d',subj_idx),32,'mean narrow filtered BOLD across brain'],'fontsize',18)
plotFilename = strcat(plotDir,'S_',sprintf('%.2d',subj_idx),'_Heart_MeanBOLDts_Narrow');
set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gcf, 'PaperPositionMode', 'auto');
print ('-dpng', '-painters', eval('plotFilename'))
print ('-depsc2', '-painters', eval('plotFilename'))
saveas(meanBOLDtsPlot,strcat(plotFilename,'.fig'))    


clear filteredMRI

% filteredMRI = filteredMRI(cfgMain.beginCut:cfgMain.endCut,:);

% Compute sPLV for surrogate
% ONLY IN VOXELS INSIDE BRAIN
% empPLV = zeros (1,length(insideBrain))'; % empty 3d Volume for storing empirical PLV
empPLV = abs (mean (exp (1i* (bsxfun (@minus , phaseMRI, phaseHeart'))))); % get PLV
% bsxfun applies the operation in @minus from the vector of the third input
% to each column of the matrix of the second input

plotempPLVdist = figure
nhist(empPLV)
xlabel('PLV')
title(['S',sprintf('%.2d',subj_idx),32,'empiricalPLV BrainHeart  BrainMean:' num2str(mean(empPLV))],'fontsize',18)
plotFilename = strcat(plotDir,'S_',sprintf('%.2d',subj_idx),'_Heart_EMPPLV_dist_Narrow');
set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gcf, 'PaperPositionMode', 'auto');
print ('-dpng', '-painters', eval('plotFilename'))
print ('-depsc2', '-painters', eval('plotFilename'))
saveas(plotempPLVdist,strcat(plotFilename,'.fig'))    





% clear other variables
% clear  phaseMRI

%store EMPPLV
filename_empPLVHeart = strcat(rootDir,'Subject',sprintf('%.2d',subj_idx),'\Timeseries','\Heart\empPLV_heart_narrow_subject_',sprintf('%.2d',subj_idx))

empPLV2write = zeros (53,63,46); % empty 3d Volume for storing empirical PLV
empPLV2write = empPLV2write(:); % transformed into a vector
empPLV2write(insideBrain) = empPLV; % get PLV
% bsxfun applies the operation in @minus from the vector of the third input
% to each column of the matrix of the second input
PLV3D = reshape(empPLV2write,53,63,46); % reshape it from vector to matrix
tools_writeMri(PLV3D,filename_empPLVHeart)

%% 
% store PLV in surrogate PLV matrix


surrogatePLVmatrix = zeros(360,length(insideBrain));



indexRotations=31:390; % Rotating at least one minute (30 TR = 60s) at the beggining or end
rotatedPhaseHeart = zeros(length(indexRotations),length(phaseHeart));
for iRotation = 1 : length(indexRotations)
rotatedPhaseHeart(iRotation,:) = circshift(phaseHeart,[0 indexRotations(iRotation)]);
end


for iRotation = 1 : length(indexRotations)
currentPhaseHeart = rotatedPhaseHeart (iRotation,:) ;
phaseDifference = bsxfun (@minus , phaseMRI, currentPhaseHeart');
PLV =    abs (mean (exp (1i* phaseDifference ) ) ); % 
surrogatePLVmatrix(iRotation,:) = PLV;% timeseries_get_PLV(phaseMRI,rotatedPhaseEGG(iRotation,:)'); 

disp('Rotation number for subject:')
disp(iRotation)
disp (subj_idx)
end



filename_allSurrogates = strcat(rootDir,'Subject',sprintf('%.2d',subj_idx),'\Timeseries','\Heart\allsurrogatePLV_heart_narrow_subject_',sprintf('%.2d',subj_idx),'.mat')
save(filename_allSurrogates,'surrogatePLVmatrix')


median_sPLV = median(surrogatePLVmatrix);

plotSurrPLVdist = figure
nhist(median_sPLV)
xlabel('PLV')
title(['S',sprintf('%.2d',subj_idx),32,'surrogatePLV across bain. Mean:' num2str(mean(median_sPLV))],'fontsize',18)
plotFilename = strcat(plotDir,'S_',sprintf('%.2d',subj_idx),'_Heart_surrPLV_dist_Narrow');
set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gcf, 'PaperPositionMode', 'auto');
print ('-dpng', '-painters', eval('plotFilename'))
print ('-depsc2', '-painters', eval('plotFilename'))
saveas(plotSurrPLVdist,strcat(plotFilename,'.fig'))    



surrPLV = zeros (53,63,46); % empty 3d Volume for storing empirical PLV
surrPLV = surrPLV(:); % transformed into a vector
surrPLV(insideBrain) = median_sPLV; % get PLV
% bsxfun applies the operation in @minus from the vector of the third input
% to each column of the matrix of the second input
PLV3D = reshape(surrPLV,53,63,46); % reshape it from vector to matrix

filename_medianSurrogate  = strcat(rootDir,'Subject',sprintf('%.2d',subj_idx),'\Timeseries','\Heart\surrPLV_heart_narrow_subject_',sprintf('%.2d',subj_idx),'.nii')

tools_writeMri(PLV3D,filename_medianSurrogate)
end


% once loop ends,store
subjects = [8 9 10 13 15 16 17 18 19 20 21 22 23 24 25 26 29 31 32 33 34 35 36 39 40 41 43]

cfgMain =global_getcfgmain;
cfgMain.clusterAlpha = 0.001; % exact one-sided p value that will be used for the cluster forming threshold XXXX check
cfgMain.numberofrandomizations = 1000
timeseries_statsCluster_Heart_narrow(subjects,cfgMain)
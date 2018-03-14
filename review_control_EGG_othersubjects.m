
% review_control_EGG_othersubjects



%% Parameters
subjects = global_subjectList;
cfgMain = global_getcfgmain;
insideBrain = tools_getIndexBrain('inside');
indNoBrain = tools_getIndexBrain('outside');

peaksAllsubjects = global_getEGGpeaks;
rootDir= strcat(global_path2root,'\Subjects\');

for iRealSubject = 2:30

   subj_idx = subjects(iRealSubject)


% subjects which are not this subject

surrogatePLVmatrix = zeros(29,length(insideBrain));

%% Load BOLD
% load BOLD of subject unfiltered
filename_bold_input = global_filename(subj_idx,cfgMain,'filename_csfr_Residuals_FB');
load(filename_bold_input)

%%
% Start loop

subjListForAnalysis = subjects~= subj_idx
SurrogateSubjectList = subjects(subjListForAnalysis)

for iSurrogateSubject = 1:29
    
    tic
    
currentSurrogateSubject = SurrogateSubjectList(iSurrogateSubject)

% load EGG first other subject
EGGPhaseXVolumeFilename = global_filename(currentSurrogateSubject,cfgMain,strcat('EGGPhaseXVolumeFilename'));
load(EGGPhaseXVolumeFilename)

% Frequency info of surrogate EGG
indPeak = find (peaksAllsubjects(:,1) == currentSurrogateSubject);
mostPowerfullFrequency = peaksAllsubjects(indPeak,3);


%Filter
centerFrequency = mostPowerfullFrequency; %
filter_frequency_spread=cfgMain.frequencySpread/1000; % In hz
sr = 0.5; % 1 TR = 2s
filterOrder=(cfgMain.fOrder*fix(sr/(centerFrequency-filter_frequency_spread))-1);%in nsamples
transition_width= cfgMain.transitionWidth/100; % in normalised units
filteredMRI=tools_bpFilter(error_csf_z,sr,filterOrder,centerFrequency,filter_frequency_spread,transition_width,cfgMain.filterType);
phaseMRI = hilbert(filteredMRI);
phaseMRI = angle(phaseMRI(cfgMain.beginCut:cfgMain.endCut,:)); % Cut data to have the same length as EGG (cut this way to get rid of fmri edge artifact on EGG)
% error_csf_z

clear filteredMRI

% filteredMRI = filteredMRI(cfgMain.beginCut:cfgMain.endCut,:);

% Compute sPLV for surrogate
% ONLY IN VOXELS INSIDE BRAIN
% empPLV = zeros (1,length(insideBrain))'; % empty 3d Volume for storing empirical PLV
empPLV = abs (mean (exp (1i* (bsxfun (@minus , phaseMRI, angle (phaseXVolume)'))))); % get PLV
% bsxfun applies the operation in @minus from the vector of the third input
% to each column of the matrix of the second input

% clear other variables
clear  phaseMRI

% store PLV in surrogate PLV matrix
surrogatePLVmatrix(iSurrogateSubject,:) = empPLV;

toc

end


filename_allSurrogates = strcat(rootDir,'Subject',sprintf('%.2d',subj_idx),'\Timeseries','\ControlEGG_othesubject\distribution_sPLV_s',sprintf('%.2d',subj_idx),'.mat')
save(filename_allSurrogates,'surrogatePLVmatrix')

median_sPLV = median(surrogatePLVmatrix);


surrPLV = zeros (53,63,46); % empty 3d Volume for storing empirical PLV
surrPLV = surrPLV(:); % transformed into a vector
surrPLV(insideBrain) = median_sPLV; % get PLV
% bsxfun applies the operation in @minus from the vector of the third input
% to each column of the matrix of the second input
PLV3D = reshape(surrPLV,53,63,46); % reshape it from vector to matrix

filename_medianSurrogate  = strcat(rootDir,'Subject',sprintf('%.2d',subj_idx),'\Timeseries','\ControlEGG_othesubject\median_sPLV_s',sprintf('%.2d',subj_idx),'.nii')

tools_writeMri(PLV3D,filename_medianSurrogate)
end
% once loop ends,store

cfgMain.numberofrandomizations = 1000
timeseries_statsCluster_Regression_surrogateSubject(subjects,cfgMain)
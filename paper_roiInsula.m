%{
paper_ROIinsula
Perform analysis presneted in section: Sub-threshold coupling in the right posterior insula
Load Deen 2011 parcelletion of right and left anterior, posterior and
middle insula. 
Compute for each subject average BOLD time-series in each of these ROIs, calculate empirical
and chance and perform a ttest comparing them

Inputs are stored in; Y:\MasksFromOtherPapers\Deen_Insula\
Output, results of the ttest in the command line, bonferroni corrected for multiple
comparisons

IR 04/07/2017
%}

%% Load masks
insideBrain = tools_getIndexBrain('inside');

% Load masks
rPI = ft_read_mri([global_path2root 'MasksFromOtherPapers' filesep 'Deen_Insula' filesep 'RPI_3mm.nii']);
rPI = logical(rPI.anatomy(:));
rPI = rPI(insideBrain);

rDAI = ft_read_mri([global_path2root 'MasksFromOtherPapers' filesep 'Deen_Insula' filesep 'RDAI_3mm.nii']);
rDAI = logical(rDAI.anatomy(:));
rDAI = rDAI(insideBrain);

rVAI = ft_read_mri([global_path2root 'MasksFromOtherPapers' filesep 'Deen_Insula' filesep 'RVAI_3mm.nii']);
rVAI = logical(rVAI.anatomy(:));
rVAI = rVAI(insideBrain);

lPI = ft_read_mri([global_path2root 'MasksFromOtherPapers' filesep 'Deen_Insula' filesep 'LPI_3mm.nii']);
lPI = logical(lPI.anatomy(:));
lPI = lPI(insideBrain);

lDAI = ft_read_mri([global_path2root filesep 'MasksFromOtherPapers' filesep 'Deen_Insula' filesep 'LDAI_3mm.nii']);
lDAI = logical(lDAI.anatomy(:));
lDAI = lDAI(insideBrain);

lVAI = ft_read_mri([global_path2root filesep 'MasksFromOtherPapers' filesep 'Deen_Insula' filesep 'LVAI_3mm.nii']);
lVAI = logical(lVAI.anatomy(:));
lVAI = lVAI(insideBrain);

% Load the information about the peaks of the EGG
peaksAllsubjects = global_getEGGpeaks;

%initialize
empPLV_ROI = zeros(6,length(subjects));
chancePLV_ROI= zeros(6,length(subjects));

%% Loop subject
for iSubj = 1:length(subjects)
    subj_idx = subjects(iSubj)
%% Load and filter BOLD timeseries from subject 

EGGPhaseXVolumeFilename = global_filename(subj_idx,cfgMain,strcat('EGGPhaseXVolumeFilename'));
load(EGGPhaseXVolumeFilename)

BOLDTimeseriesFilename = global_filename(subj_idx,cfgMain,strcat('filename_',cfgMain.Timeseries2Regress,'_Residuals_FB'));
timeseries = load(BOLDTimeseriesFilename);

indPeak = find (peaksAllsubjects(:,1) == subj_idx);
mostPowerfullFrequency = peaksAllsubjects(indPeak,3);


centerFrequency = mostPowerfullFrequency; %
filter_frequency_spread=cfgMain.frequencySpread/1000; % In hz
sr = 0.5; % 1 TR = 2s
filterOrder=(cfgMain.fOrder*fix(sr/(centerFrequency-filter_frequency_spread))-1);%in nsamples
transition_width= cfgMain.transitionWidth/100; % in normalised units
filteredMRI=tools_bpFilter(eval(strcat('timeseries.',cell2mat(fields(timeseries)))),sr,filterOrder,centerFrequency,filter_frequency_spread,transition_width,cfgMain.filterType);

filteredMRI = filteredMRI(cfgMain.beginCut:cfgMain.endCut,:);

%% Timeseries of ROI

timeseriesROI = zeros(420,6);
timeseriesROI(:,1) = nanmean(filteredMRI(:,rPI),2);
timeseriesROI(:,2) = nanmean(filteredMRI(:,rDAI),2);
timeseriesROI(:,3) = nanmean(filteredMRI(:,rVAI),2);
timeseriesROI(:,4) = nanmean(filteredMRI(:,lPI),2);
timeseriesROI(:,5) = nanmean(filteredMRI(:,lDAI),2);
timeseriesROI(:,6) = nanmean(filteredMRI(:,lVAI),2);

phaseROI = hilbert (timeseriesROI);

%% EmpPLV

    phaseDifferenceEmp = bsxfun (@minus , angle(phaseROI), angle(phaseXVolume )');

    PLV_emp =  abs (mean (exp (1i* phaseDifferenceEmp ) ) ); %
    
%% rotations PLV
indexRotations=31:390; % Rotating at least one minute (30 TR = 60s) at the beggining or end
rotatedPhaseEGG = zeros(length(indexRotations),length(phaseXVolume));
for iRotation = 1 : length(indexRotations)
rotatedPhaseEGG(iRotation,:) = circshift(phaseXVolume,[0 indexRotations(iRotation)]);
end

RPLV = zeros(length(indexRotations),6); % R from rotated

for iRotation = 1 : length(indexRotations)
    currentPhaseEGG = rotatedPhaseEGG (iRotation,:) ;
    phaseDifference = bsxfun (@minus , angle(phaseROI), angle(currentPhaseEGG )');
    PLV =  abs (mean (exp (1i* phaseDifference ) ) ); %
    RPLV(iRotation,:) = PLV;% timeseries_get_PLV(phaseMRI,rotatedPhaseEGG(iRotation,:)');
end
medianPLV = median(RPLV,1);

empPLV_ROI (:,iSubj)= PLV_emp
chancePLV_ROI(:,iSubj)= medianPLV;


end

%% ttest
[h p ci stats] = ttest2(empPLV_ROI',chancePLV_ROI');

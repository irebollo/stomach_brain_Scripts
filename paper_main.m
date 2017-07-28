%{
% paper_main

Scripts that reproduces the analysis on physiens experiment for
publishing

0 - Set global variables
1 - MRI preprocessing with SPM
2 - preprocessing of EGG data
2B - Control 1: EGG Peaks outside the scanner
3 - timeseries analisis
4 - Group level statistics of Gastric-BOLD coupling
4b -Table1
5 - Offset control
6 - Randomization of rotations control
7 - Variance explained in gastric network clusters
8 - Control FWD displacement
9 - Bayes stats
10 - Insula ROI post-hoc and autonomic overlap
11 - Phase organization in the gastric network

Ignacio Rebollo 05/07/2017
%}
%% 0 set global variables
cfgMain = global_getcfgmain;
subjects= global_subjectList;

%% 1 MRI preprocessing 
% ATTENTION THIS CELL SHOULD BE RUN WITH MATLAB 2011, WHILE THE REST IS RUN
% WITH MATLAB 2013
% Calls SPM 8 

% slice timing correction,  normalization to MNI template, motion correction spatial smoothing

for iSubject=1:length(subjects)
    mri_group_preproc(subjects(iSubject),cfgMain.kernelWidth)
end

%% 2 - preprocessing of EGG data
% from now on MATLAB 2013
cfgMain.automaticChannelSelection = 1 % the first time is set to 0 then after calling prepro_egg_saveChannelInfo(cfgMain) is turn to 1
cfgMain.plotFigures= 1
for iSubj = 1: length(subjects)
    close all
    subj_idx = subjects(iSubj)
        prepro_egg(subj_idx,cfgMain)
end

% After calling one time and choosing and storing the best channel and most power frequency
% of each subject manually, automatic channel selection is turned on to use the
% manually selected info for future EGG analysis (E.g. filtering BOLD and 
% offset control)

prepro_egg_saveChannelInfo(cfgMain) 
% cfgMain.automaticChannelSelection = 1

%% 2b - 1st control peaks outside the scaner are the same inside the scaner ?

paper_controlPeaksOustideScanner(subjects,cfgMain)


%% 3 - timeseries analisis
cfgMain.plotFigures= 0

for iSubj = 1:length(subjects)
    subj_idx = subjects(iSubj)
timeseries_prepare_import2matlab(subj_idx,cfgMain) 
% Concatenate the swaf image into a matrix containing each voxel timeseries, done once per
%  kernelwidth
timeseries_preprocessBOLD(subj_idx,cfgMain) % preprocess the timeseries (remove polinomial and filter)
timeseries_csfSignal_obtainAndRegress(subj_idx,cfgMain) % obtain the wm and csf signals from the BOLD preprocessed timeseries and stores them
timeseries_preparePhases_Regression(subj_idx,cfgMain) % Filter and hilbert transform residuals of csf regression
timeseries_mapPLV_Regression(subj_idx,cfgMain) % Obtain PLV per voxel
timeseries_medianRotation_Regression(subj_idx,cfgMain) % Obtain chance PLV per voxel  

end

%% 4 - Group level statistics of Gastric-BOLD coupling

timeseries_statsCluster_Regression(subjects,cfgMain)


%% 4 b Anatomical characterizaion of gastric-BOLD network(table 1)

resultsCharacterization = paper_MakeTable1(cfgMain)

%% 5- offset control

paper_controlOffset_MAIN

%% 6 - Randomization of rotations control 

% calculate PLV for all rotations
for iSubj = 1:length(subjects)
    subj_idx = subjects(iSubj)
timeseries_AllRotations_Regression(subj_idx,cfgMain) 
end
paper_randomization_randomize
paper_randomization_figure

%% 7 - Variance explained in gastric network clusters
   
networks_coherence_main(cfgMain)

%% 8- Control on FWD displacement

for iSubj = 1:length(subjects)
    subj_idx = subjects(iSubj)
networks_couplingStrenght_obtain(subj_idx,cfgMain) 
timeseries_FWD_BOLD_GLM(subj_idx,cfgMain) 
end

paper_control_corrMovement_CS

%% 9- Bayes stats
paper_bayes_MAIN

%% 10 -  Insula ROI post-hoc and autonomic overlap

paper_ROIinsula
paper_quantifyAutonomicOverlap

%% 11- Phase organziation in gastric network

% obtain phases angles of each subject

for iSubj = 1:length(subjects)
    subj_idx = subjects(iSubj)
networks_phaseAngle_obtain(subj_idx,cfgMain)
end

Paper_ComparisonSharedVariance
Paper_ComparisonSharedVarianceAllNodes

function timeseries_mapPLV_Regression(subj_idx,cfgMain)
%{
Computes and store in timeseries/data folder a 3D nifti volume with Phase locking value between the EGG and each
voxel timeseries


inputs:
subj_idx = s number
cfgMain must contain fields
    kernelWidth,Timeseries2Regress,frequencySpread ,fOrder,beginCut,endCut
kernelWidth: with of the smoothing kernel from preprocessing, paper  = % 3mm
cfgMain.Timeseries2Regress should be 'csf' to load residuals of csf regression
fOrder : multiplicative factor of the filter order
frequencySpread: spead of the time domain filter in hz * 1000, paper = 0.015 hz = 15,
begin and end cut are the voulmes that are discarded to avoid the filter
ringing artifact
cfgMain.transitionWidth is the transition width of the filter, paper is 15
offset is with respect to EGG peaking filter, only for control analysis.
offset is in hz x 1000 e.g. and offset of 0.006 hz is a value of 6

It's input is the output of the script timeseries_preparePhases_Regression 
Y:\Subjects\Subject13\Timeseries\MRItimeseries\csfResiduals_FB_phases_s13_kw3_fir2_fspread_015
and EGG phases
Y:\Subjects\Subject13\Timeseries\EGGtimeseries\PhaseXvolume_S_13_fir2_fspread_015_ord_5_tw_15

Output: saves data in subject timeseries folder as a 3D .nii image
Y:\Subjects\Subject13\Timeseries\PhasesAnalysis\PLVxVoxel_csfr_S_13_kw3_fir2_fspread_015_fOrder_5_tw_15

It uses mike cohen matlab implementation of PLV (2013)

IR commented on 28/06/2017
%}

%% Pass parameters of cfgMain to function, load data and set output filenames

fOrder = cfgMain.fOrder;
frequencySpread = cfgMain.frequencySpread;
kernelWidth= cfgMain.kernelWidth;

% Get EGG phase x volume
EGGPhaseXVolumeFilename = global_filename(subj_idx,cfgMain,strcat('EGGPhaseXVolumeFilename'));
load(EGGPhaseXVolumeFilename)

PLVXVoxelFilename = global_filename(subj_idx,cfgMain,strcat('PLVXVoxelFilename_',cfgMain.Timeseries2Regress)); % output filename
BOLDPhasesTimeseriesFilename = global_filename(subj_idx,cfgMain,strcat('filename_',cfgMain.Timeseries2Regress,'_Residuals_FB_phases')); % input filename

load (BOLDPhasesTimeseriesFilename)

plotDir = strcat (global_path2subject(subj_idx),'PreprocessingLog',filesep);
plotFilename = strcat(plotDir,'S_',sprintf('%.2d',subj_idx),'_PLVhistogram');

[indNoBrain] = tools_getIndexBrain('outside');
[indBrain] = tools_getIndexBrain('inside');
%% Calculate PLV

empPLV = zeros (53,63,46); % empty 3d Volume for storing empirical PLV
empPLV = empPLV(:); % transformed into a vector
empPLV(indBrain) = abs (mean (exp (1i* (bsxfun (@minus , angle(phaseMRI), angle (phaseXVolume)'))))); % get PLV
% bsxfun applies the operation in @minus from the vector of the third input
% to each column of the matrix of the second input
empPLV(indNoBrain) = 0;
PLV3D = reshape(empPLV,53,63,46); % reshape it from vector to matrix

%% Save data

tools_writeMri(PLV3D,PLVXVoxelFilename)

%% Sanity plot

if cfgMain.savePlots == 1
    
    insideBrain = tools_getIndexBrain('inside');
    
    voxelCoordinates = sub2ind([53,63,46],11,30,37); % Right somatosensory cortex
    voxelCoordinates_inside = zeros(153594,1);
    voxelCoordinates_inside(voxelCoordinates)=1;
    voxelCoordinates_inside = voxelCoordinates_inside(insideBrain);
    ind_voxelCoordinates_inside = find(voxelCoordinates_inside);
    
    if cfgMain.plotFigures == 0;
        SanityPlot = figure('visible','off');
    else
        SanityPlot = figure('visible','on');
    end
    
    % Plot histogram of PLV across the brain
    
    nhist(empPLV(indBrain))
    xlabel('PLV')
    title(['S',sprintf('%.2d',subj_idx),32,'PLV across bain. Mean:' num2str(mean(empPLV(indBrain))) ' rSS voxel:' 32 num2str(empPLV(indBrain(ind_voxelCoordinates_inside)))],'fontsize',18)
    
    
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    set(gcf, 'PaperPositionMode', 'auto');
    
    print ('-dpng', '-painters', eval('plotFilename'))
    print ('-depsc2', '-painters', eval('plotFilename'))
    saveas(SanityPlot,strcat(plotFilename,'.fig'))
    
end

end
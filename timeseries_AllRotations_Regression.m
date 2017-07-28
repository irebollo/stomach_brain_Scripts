function timeseries_AllRotations_Regression(subj_idx,cfg)
%{

Stores a matlab file with the PLV obtained from all 360 possible rotations
of EGG signal. The file is later used in the false positive rate control
script paper_control_randomziationTimeShifts

inputs:
subj_idx = s number
cfg must contain fields
    kernelWidth,Timeseries2Regress,frequencySpread ,fOrder,beginCut,endCut
kernelWidth: with of the smoothing kernel from preprocessing, paper  =
3mm
cfg.Timeseries2Regress should be 'csf' to load residuals of csf
regression
fOrder : multiplicative factor of the order of the filter 
frequencySpread: spead of the time domain filter in hz * 1000, paper = 0.015 hz = 15,
begin and end cut are the voulmes that are discarded to avoid the filter
ringing artifact
cfg.transitionWidth is the transition width of the filter, paper is 15 
offset is with respect to EGG peaking filter, only for control analysis.
offset is in hz x 1000 e.g. and offset of 0.006 hz is a value of 6

Input file
BOLD
Y:\Subjects\Subject13\Timeseries\MRItimeseries\csfResiduals_FB_phases_s13_kw3_fir2_fspread_015
and EGG phases
Y:\Subjects\Subject13\Timeseries\EGGtimeseries\PhaseXvolume_S_13_fir2_fspread_015_ord_5_tw_15

Output: saves PLV of each rotation  in subject timeseries folder as a roationXvoxel matlab
file
Y:\Subjects\Subject13\Timeseries\PhasesAnalysis\AllRotationsPLV_csfr_S_13_kw3_fir2_fspread_015_fOrder_5_tw_15

IR commented 28/06/2017

%}

%% Import cfg parameters
fOrder = cfg.fOrder;
frequencySpread = cfg.frequencySpread;
kernelWidth= cfg.kernelWidth;
offset = cfg.offset;

outsideBrain = tools_getIndexBrain('outside');
insideBrain = tools_getIndexBrain('inside');

%% Load data and outputfilename

AllRotationFilename = global_filename(subj_idx,cfg,strcat('AllRotationsFilename_',cfg.Timeseries2Regress)); % output
BOLDPhasesTimeseriesFilename = global_filename(subj_idx,cfg,strcat('filename_',cfg.Timeseries2Regress,'_Residuals_FB_phases')); % input filename

load (BOLDPhasesTimeseriesFilename)

EGGPhaseXVolumeFilename = global_filename(subj_idx,cfg,strcat('EGGPhaseXVolumeFilename')); 
load(EGGPhaseXVolumeFilename)

%% Rotate EGG

indexRotations=31:390; % Rotating at least one minute (30 TR = 60s) at the beggining or end
rotatedPhaseEGG = zeros(length(indexRotations),length(phaseXVolume));
for iRotation = 1 : length(indexRotations)
rotatedPhaseEGG(iRotation,:) = circshift(phaseXVolume,[0 indexRotations(iRotation)]);
end


%% Calculate PLV
disp('+++++++++++++++++++++++++++++++ RPLV')

% initialize structure for distribution of rotated PLV values
RPLV = zeros(length(indexRotations),153594); % R from rotated

%iterate through all rotation and calculate PLV

for iRotation = 1 : length(indexRotations)
    

PLV = zeros (53,63,46);
PLV = PLV(:);

currentPhaseEGG = rotatedPhaseEGG (iRotation,:) ;
phaseDifference = bsxfun (@minus , angle(phaseMRI), angle(currentPhaseEGG )');
PLV(insideBrain) =    abs (mean (exp (1i* phaseDifference ) ) ); % 

RPLV(iRotation,:) = PLV;% timeseries_get_PLV(phaseMRI,rotatedPhaseEGG(iRotation,:)'); 

disp('Rotation number for subject:')
disp(iRotation)
disp (subj_idx)
end



save(AllRotationFilename,'RPLV')

end
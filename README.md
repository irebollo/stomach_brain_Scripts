# stomach_brain_Scripts
Scripts to reproduce the analysis described by Rebollo et al on stomach brain synchronisation

General organization of scripts:

The scripts are organized according to prefixes

Global scripts: Global variables, like subject indexes, path, filenames
MRI scripts: Preprocessing of mri antomical and functional images using SPM8
Prepro scripts: preprocessing of the EGG inside and outside the scanner
timeseries scripts: Analysis of BOLD and EGG timeseries, estimation of chance and empirical PLV, group level statistics
Tools scripts: bandpass filter, brain atlas lookup, write mri, vox2mni,mni2vox
Network scripts: Coherence and phase angle analysis of gastric network
Paper scripts: Scripts to make specific analysis for figures for the paper


The main script is paper_main, from which all the analysis can be launched. All analysis were performed with Matlab 2013 (and SPM 8) except one step of fMRI preprocessing, documented in paper_main

Global scripts
global_path2root
path to Physiens root folder

global_path2subject
path to subject specific folder

global_filename
Function that defines the path of all files in PHYSIENS study by concatenating the different elements composing the path and filename. The output of the function is a string containing the complete path to the requested file that can be used to load the file

Global_subjectList
List of all artifact free subjects included on PHYSIENS final analysis

global_getcfgmain
Function defining a matlab structure containing the final parameters used through the preprocessing and analysis pipeline. This structure is a mandatory input to most PHYSIENS functions, either to know which parameter apply or which filename retrieve/store using global_filename. In order to retrieve the path to a file, the different parts of the cfgMain
structure are concatenated when calling global_filename

global_getEGGpeaks
Load the file containing the channel number and exact frequency for all subjects during the fMRI acquisition, which was manually estimated in the first call of prepro_egg and stored in the log file together with the EGG phase timeseries in each subject timeseries folder. The file to be loaded was created by the function prepro_egg_saveChannelInfo

global_getEGGpeaksOutsideScanner
Retrieve the frequency of the EGG peak recorded when subjects where outside the scanner in order to do the analysis performed in paper_controlPeaksOustideScanner

global_getGastricNetwork
Retrieves and load the mask of the gastric network obtained from 30 subjects. 

MRI
MRI scripts are run with matlab 2011 and SPM8. These scripts needs the functions in the folder MRITools

MRI_group_preproc

Toplevel batch for preprocessing with SPM8 by Anne-Dominique Lodeho. To use this batch, please edit the beginning of this file and run it. 
you need to specify which particular tasks to perform, through the TODO structure with binary fields (1|0):
    todo.slice_timing = 1;
    todo.realign      = 1;
    todo.normalize    = 1;
    todo.coregister   = 1;
    todo.apply_norm   = 1;
    todo.smooth       = 1;
    todo.run          = 1;
This functions iterates through subjects and calls MRI_single_preproc that does the actual call to SPM

MRI_single_preproc
Create an SPM8 job for single subject's preprocessing by Anne-Dominique Lodeho with the parameters specified in MRI_group_preproc

Prepro scripts 
Prepro_egg
Basic preprocessing steps for obtaining and storing the EGG phase per mri volume

prepro_egg_loadData
Function to retrieve subject EGG and IRM volume markers recorded inside or outside the scanner. This function is called by prepro_egg

prepro_egg_path2data
This function get the full path of different kind of files, e.g. brainamp markers, for a particular subject. This function is called by prepro_egg_loadData

prepro_egg_downsampleMarkers
Get closest sample number of new marker based on the timestamps of the not downsampled data

prepro_egg_saveChannelInfo
Create a file in the folder scripts/files called EGG_peaks_info based on the log of the preprocessing of each subject, it retrieves the mostpowerfull frequency and channel of each subject and put it on a table. The original selection of channels and frequency peaks is done on the first call to prepro_egg, with the paramter cfgMain.automaticChannelSelection == 0
	
prepro_egg_outsideScanner
Similar to prepro_egg buy preprocess the brainamp files containing the EGG recorded outside the scanner.

timeseries scripts
timeseries_prepare_import2matlab
Construct fMRI timeseries (time,voxels) from swaf images located at each subject folder and saves them as a structure in the HDD.

timeseries_preprocessBOLD
Takes the output of timeseries_prepare_import2matlab that is stored in each subject timeseries\data folder and preprocess it so it can be further analyzed Preprocessing steps, implemented in fieldtrip, includes remove polinomial trends of second degree, bandpass filter between 0.01 and 0.1 Hz and standarize to Z units the output is stored in the tiemseries folder of each subject

timeseries_csfSignal_obtainAndRegress
This function retrieves the timeseries of a 3x3 cube of voxels located in the CSF and regress it out the of all the brain voxels, the timeseries, the residuals of the regression and the betas coefficient are all stores in the timeseries folder of each subject.

timeseries_preparePhases_Regression
Load EGG and FMRI timeseries, filter fMRI timeseries at gastric peak, extract phases of fmri timeseries(CSFregressed) and save them into disk in timeseries folder.

timeseries_mapPLV_Regression
Computes and store in timeseries/data folder a 3D nifti volume with phase locking value between the EGG and each voxel timeseries

timeseries_medianRotation_Regression
Stores a .nii image of the median PLV obtained with a timeshifted EGG signals concatenated and rotated). It does all 360 possible rotations of EGG signal and calculate PLV at each
rotations and takes the median PLV obtained across timeshift across voxels to be used as surrogate PLV in group level statistics

timeseries_AllRotations_Regression
Stores a matlab file with the PLV obtained from all 360 possible rotations of EGG signal. The file is later used in the false positive rate control script paper_control_randomziationTimeShifts

timeseries_statsCluster_Regression
Perform group level statistics comparing empirical versus chance PLV by using the clustering randomization procedure provided by fieldtrip (Maris and Oostenveld 2007)

timeseries_FWD_BOLD_GLM
Store the timeseries of absolute derivative of movement (separately for rotation and translations) i.e. Framewise displacement. Then compute the GLM between the BOLD timeseries and the framewise displacement. To be later used by the script paper_control_corrMovement_CS.Store figure into disk for each subject.
Networks scripts
networks_coherence_getClustertimeseries
Compute and stores the average BOLD time series across each significant cluster of the gastric network. Later used for coherence and angle analysis.

networks_coherence_getFullBandEGG
Compute and stores the fullband (0.01-0.10 EGG time series in order to compute coherence with the BOLD signal in significant clusters.

networks_coherence_getCoherence
This function computes the coherence between the EGG and all gastric network regions and between all gastric network regions for a given participant and store it into disk
networks_coherence_group
This Script loads the BOLD-EGG coherence values of all subjects in a loop and average them at gastric peaking frequency separately for each subject. The resulting values are used in table 1 of the paper

networks_coherence_main
Simple script that calls the relevant functions to perform the coherence analysis between gastric network and the EGG at the group level.

networks_couplingStrenght_obtain
Computes difference between median and empirical PLV (coupling strength) and stores in the HDD. These values are later used in the control with head micro movements

networks_phaseAngle_obtain
This function 
1 Loads cluster timeseries for each subject and EGG
2 Get angle and value of phase locking for each cluster
3 Subtract mean phase across clusters to the angle of each cluster
4 Plots previous steps and store them into HDD
5 Stores the angle and magnitude of PLV for each cluster

networks_phaseAngle
This scripts loads the angle of phase locking of each cluster and subject, average them, perform the watson-william statistical test, plot them, and store them in nifti format (Figure 4 of the paper).

Paper scripts
paper_controlPeaksOustideScanner
Control to check if the EGG peaks are the same inside and outside the scanner. First calls prepro_egg_outsideScanner to preprocess EGG outside the scaner. Manual corrections in 4 subjects is applied. Then perform a ttest between peaks inside and outside the scanner. 

paper_MakeTable1
The output of this function is used in table 1 of the paper to identify the regions belonging to the gastric network using automatic anatomical labeling(AAL).
Input:  cfgMain structure used to obtain results, which allows to load the group
level masks. The output is the area number, the number and volume of voxels
in that area, the percentage of the area in the cluster, the tvalue at peak
coordinates and the x y z MNI coordinates of the peak

paper_controlOffset_MAIN
Perform filter offset control. Filters BOLD and EGG timeseries with frequency offsets from the EGGpeak in steps of 0.001 Hz, calculates empirical and chance PLV and performs group
level statistics. It then loads the group level stats of all offsets and plots them to be
used in figure 2B Modifying the parameter cfgMain.offset changes the filenames of all relevant files when calling global_filename in order to reflect the desired filter offset.

paper_ROIinsula
Perform analysis presented in section: Sub-threshold coupling in the right posterior insula
Load Deen 2011 parcellation of right and left anterior, posterior and
middle insula. Compute for each subject average BOLD time-series in each of these ROIs, calculate empirical and chance and perform a ttest comparing them.

Paper_ComparisonSharedVariance
This script produces the analysis presented in figure 4D It compares the estimates of functional connectivity (shared variance) between gastric, saliency and default regions using coherence and correlations It thens plots the timeseries of an exempe subject.

Paper_ComparisonSharedVarianceAllNodes
This scripts compares estimates of shared variance between all nodes of the gastric network using the squared coherence and squared correlation coefficient. It then performs a ttest to check for significant differences. The resulting values are reported in the text of the paper, section 'Temporal advances and lags in the gastric network '.

paper_control_corrMovement_CS
Perform correlation between the framewise displacement betas and coupling
Strength in all voxels of the gastric network. Reported in section Gastric-coupling in gastric network is not related to head micromovements

paper_bayes_MAIN
Compute Bayes Factor for the influence of different variables on
gastric-BOLD coupling strength
1- Framewise displacement, one sample ttest
2- EGG power, frequency, traitAnxiety, BMI. Correlation
3- Gender, two sample ttest

paper_quantifyAutonomicOverlap
 This script quantify the overlap between voxels of the gastric network and the electrodermal activity related regions provided by Beisnner 2013.

paper_randomization_createRandomizationMatrix
Creates and saves a matrix defining the time-shifts that will be used for each subject at each iteration for the randomization control. To call only once, and then load the output matrix stored in Y:\scripts\files\randomizationRotationMatrix

paper_randomization_randomize
This script performs the group level statistics with PLV obtained with random timeshifts to the EGG timeseries across subjects.
The PLV obtained with all possible timeshifts is obtained with randomization of EGG phases using the script timeseries_AllRotations_Regression 
The randomization matrix is loaded, which contains for each iteration (1:1000), which PLV of time-shifted EGG (out of 360), will be assigned to each of the 30 subjects
It then calls paper_randomization_statsCluster to perform group level statistics on these surragate data.

paper_randomization_statsCluster
This script performs the group level statistics for the randomization
control. Instead of comparing the median rotation PLV against the empirical
PLV, it compares the median rotation PLV against a random rotation, which
is previously determined in the file Y:\scripts\files\randomizationRotationMatrix

paper_randomization_figure

Tools scripts
tools_get_indexBrain
Get index of voxel outside/inside brain in vector format based on SPM apriori mask

tools_bpFilter
Forward and reverse filter the signal using filtfilt. This corrects for phase distortion introduced by a one-pass filter. filterType has to be specified as string input e.g. ’fir2’

tools_vox2mni
Voxel space to MNI space

tools_mni2vox
MNI space to Voxel space 

tools_writeMri
Write 3d or 4d data in nifti file using fieldtrip

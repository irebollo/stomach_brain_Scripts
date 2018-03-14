function filename =  global_filename (subj_idx,cfg,fileType)
% (subj_idx,cfg,fileType)
%{
Define the path of all files in PHYSIENS study by concatenating the
different elements composing the path and filename.
The output of the function is a string containing the complete path to the requested
file that can be used to load the file

The input consist of:
subj_idx: subject number, if a file that is not specific to a subject is
used this input should be 0
cfg: configuration structure, should include all parameters that are
going to be in the filename e.g. cfg.frequencySpread

fileType: a string containign the type of file the function has to search

example call
 path2medianRotation =  global_filename (13,cfg,'medianRotationFilename')
 path2GroupLevelStatisticMap =  global_filename (0,cfg,'clusterOutputFilename')

IR, commented on the 15/09/2016
And checked it on the 28/06/2017

%}
%% Define the nature of the requested filed

SubjectDataRoot = strcat(global_path2subject(subj_idx),'Timeseries',filesep);

fileisCluster = strcmp (fileType,'clusterOutputFilename'); % group level statistics cluster files are treated separatly
fileisPPCCluster = strcmp (fileType,'clusterOutputFilenamePPC'); % group level statistics cluster files are treated separatly

fileisHeartCluster = strcmp (fileType,'HeartclusterOutputFilename'); % group level statistics cluster files are treated separatly
%% Once it knows what it lookin


if fileisCluster == 1
    
    if cfg.randomized == 1 % randomized refers to false positive control
        
        clusterOutputFilename = strcat(global_path2root,'ClusterResults',filesep,'Randomizations',filesep,'RandomCluster_N',num2str(cfg.currentRandomizationIteration),'_nR',num2str(cfg.numberofrandomizations),'_CA',sprintf('%.3d',cfg.clusterAlpha*1000),...
            '_kw',num2str(cfg.kernelWidth),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),'_fOrder_',num2str(cfg.fOrder),'_tw_',sprintf('%.2d',cfg.transitionWidth),'_offset_',sprintf('%.4d',cfg.offset*10));
        
    elseif cfg.offset ~= 0 % filter offset control
        
        clusterOutputFilename = strcat(global_path2root,'ClusterResults',filesep,'OffsetControl',filesep,'Cluster_','nR',num2str(cfg.numberofrandomizations),'_CA',sprintf('%.4d',cfg.clusterAlpha*10000),...
            '_kw',num2str(cfg.kernelWidth),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),'_fOrder_',num2str(cfg.fOrder),'_tw_',sprintf('%.2d',cfg.transitionWidth),cfg.Timeseries2Regress,'_offset_',sprintf('%.4d',cfg.offset*10));
        
    else
        if cfg.hrfconvolved == 1
            clusterOutputFilename = strcat(global_path2root,'ClusterResults',filesep,'kw',num2str(cfg.kernelWidth),filesep,'CA',sprintf('%.4d',cfg.clusterAlpha*10000),filesep,'Cluster_convolvedEGG_','nR',num2str(cfg.numberofrandomizations),'_CA',sprintf('%.4d',cfg.clusterAlpha*10000),...
                '_kw',num2str(cfg.kernelWidth),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),'_fOrder_',num2str(cfg.fOrder),'_tw_',sprintf('%.2d',cfg.transitionWidth),cfg.Timeseries2Regress);
        elseif cfg.scrubbed == 1
            clusterOutputFilename = strcat(global_path2root,'ClusterResults',filesep,'kw',num2str(cfg.kernelWidth),filesep,'CA',sprintf('%.4d',cfg.clusterAlpha*10000),filesep,'Cluster_Scrubbed_','nR',num2str(cfg.numberofrandomizations),'_CA',sprintf('%.4d',cfg.clusterAlpha*10000),...
                '_kw',num2str(cfg.kernelWidth),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),'_fOrder_',num2str(cfg.fOrder),'_tw_',sprintf('%.2d',cfg.transitionWidth),cfg.Timeseries2Regress);
            
        else
            clusterOutputFilename = strcat(global_path2root,'ClusterResults',filesep,'kw',num2str(cfg.kernelWidth),filesep,'CA',sprintf('%.4d',cfg.clusterAlpha*10000),filesep,'Cluster_','nR',num2str(cfg.numberofrandomizations),'_CA',sprintf('%.4d',cfg.clusterAlpha*10000),...
                '_kw',num2str(cfg.kernelWidth),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),'_fOrder_',num2str(cfg.fOrder),'_tw_',sprintf('%.2d',cfg.transitionWidth),cfg.Timeseries2Regress);
        end
    end
    %
    
elseif fileisPPCCluster == 1
    clusterOutputFilenamePPC = strcat(global_path2root,filesep,'ClusterResults',filesep,'PPC',filesep,'CA',sprintf('%.4d',cfg.clusterAlpha*10000),filesep,'Cluster_','nR',num2str(cfg.numberofrandomizations),'_CA',sprintf('%.4d',cfg.clusterAlpha*10000),...
        '_kw',num2str(cfg.kernelWidth),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),'_fOrder_',num2str(cfg.fOrder),'_tw_',sprintf('%.2d',cfg.transitionWidth),cfg.Timeseries2Regress);    
    
elseif fileisHeartCluster == 1
    HeartclusterOutputFilename = strcat(global_path2root,filesep,'ClusterResults',filesep,'Heart',filesep,'CA',sprintf('%.4d',cfg.clusterAlpha*10000),filesep,'Cluster_','nR',num2str(cfg.numberofrandomizations),'_CA',sprintf('%.4d',cfg.clusterAlpha*10000),...
        '_kw',num2str(cfg.kernelWidth),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),'_fOrder_',num2str(cfg.fOrder),'_tw_',sprintf('%.2d',cfg.transitionWidth),cfg.Timeseries2Regress);
    
else
    
    
    if cfg.offset ~= 0 % if it's looking to the files corresponding to the filter offset control
        
        
        EGGPhaseXVolumeFilename = strcat(SubjectDataRoot,'EGGtimeseries',filesep,...
            'PhaseXvolume_S_',sprintf('%.2d',subj_idx),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),...
            '_ord_',num2str(cfg.fOrder),'_tw_',sprintf('%.2d',cfg.transitionWidth),'_offset_',sprintf('%.4d',cfg.offset*10));
        EGGAmplitudeXVolumeFilename = strcat(SubjectDataRoot,'EGGtimeseries',filesep,...
            'EGGtimeseries_S_',sprintf('%.2d',subj_idx),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),...
            '_ord_',num2str(cfg.fOrder),'_tw_',sprintf('%.2d',cfg.transitionWidth),'_offset_',sprintf('%.4d',cfg.offset*10));
        
        %% BOLD
        
        BOLDTimeseriesFilename = strcat(SubjectDataRoot,'MRItimeseries',filesep,'fMRItimeseries_S',sprintf('%.2d',subj_idx),'_kw',num2str(cfg.kernelWidth),'_offset_',sprintf('%.4d',cfg.offset*10),'.mat');
        
        BOLDBPAmplitudeTimeseriesFilename = strcat(SubjectDataRoot,'MRItimeseries',filesep,'BOLDFiltered_S',sprintf('%.2d',subj_idx),'_kw',...
            num2str(cfg.kernelWidth),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),'_fOrder_',num2str(cfg.fOrder),'_tw_',sprintf('%.2d',cfg.transitionWidth),'_offset_',sprintf('%.4d',cfg.offset*10));
        BOLD_filtered_fullbandFilename = strcat(SubjectDataRoot,'MRItimeseries',filesep,'BOLDFULLBAND_S',sprintf('%.2d',subj_idx),'_kw',...
            num2str(cfg.kernelWidth),'_offset_',sprintf('%.4d',cfg.offset*10));
        BOLD_filtered_narrowbandFilename = strcat(SubjectDataRoot,'MRItimeseries',filesep,'BOLDnarrowBAND_S',sprintf('%.2d',subj_idx),'_kw',...
            num2str(cfg.kernelWidth),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),'_offset_',sprintf('%.4d',cfg.offset*10));
        
        filename_csf_Signal_FB = strcat(SubjectDataRoot,'MRItimeseries',filesep,'filename_csfSignal_FB_S_',sprintf('%.2d',subj_idx),'_kw',...
            num2str(cfg.kernelWidth),'_offset_',sprintf('%.4d',cfg.offset*10));
        filename_csfr_Residuals_FB =strcat(SubjectDataRoot,'MRItimeseries',filesep,'csfRegressionResiduals_FB_S_',sprintf('%.2d',subj_idx),'_kw',...
            num2str(cfg.kernelWidth));
        filename_csfr_Betas_FB = strcat(SubjectDataRoot,'MRItimeseries',filesep,'filename_csfRegressionBetas_FB_S_',sprintf('%.2d',subj_idx),'_kw',...
            num2str(cfg.kernelWidth));
        filename_csfr_Residuals_FB_phases = strcat(SubjectDataRoot,'MRItimeseries',filesep,'csfResiduals_FB_phases_s',sprintf('%.2d',subj_idx),'_kw',...
            num2str(cfg.kernelWidth),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),'_offset_',sprintf('%.4d',cfg.offset*10));
        
        
        %% Phase analysis
        
        PLVXVoxelFilename_csfr = strcat(SubjectDataRoot,'PhasesAnalysis',filesep,'PLVxVoxel_csfr_S_',sprintf('%.2d',subj_idx),...
            '_kw',num2str(cfg.kernelWidth),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),...
            '_fOrder_',num2str(cfg.fOrder),'_tw_',sprintf('%.2d',cfg.transitionWidth),'_offset_',sprintf('%.4d',cfg.offset*10));
        medianRotationFilename_csfr = strcat(SubjectDataRoot,'PhasesAnalysis',filesep,'medianRotation_csfr_S_',sprintf('%.2d',subj_idx),...
            '_kw',num2str(cfg.kernelWidth),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),'_fOrder_',num2str(cfg.fOrder),'_tw_',sprintf('%.2d',cfg.transitionWidth),'_offset_',sprintf('%.4d',cfg.offset*10));
        
    else % Normal case no offset
        
        %% EGG
        EGGPhaseXVolumeFilename = strcat(SubjectDataRoot,'EGGtimeseries',filesep,...
            'PhaseXvolume_S_',sprintf('%.2d',subj_idx),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),...
            '_ord_',num2str(cfg.fOrder),'_tw_',sprintf('%.2d',cfg.transitionWidth));
        EGGAmplitudeXVolumeFilename = strcat(SubjectDataRoot,'EGGtimeseries',filesep,...
            'EGGtimeseries_S_',sprintf('%.2d',subj_idx),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),...
            '_ord_',num2str(cfg.fOrder),'_tw_',sprintf('%.2d',cfg.transitionWidth));
        
        FilenameamplitudeXVolumeBestChannel_FULLBAND = strcat(SubjectDataRoot,'EGGtimeseries',filesep,...
            'EGGTimeseriesFullband_S_',sprintf('%.2d',subj_idx));
        
        EGGoutsideScannerFilename = strcat(SubjectDataRoot,'EGGtimeseries',filesep,...
            'EGGoutsideScannerInfo',sprintf('%.2d',subj_idx));
        
        
        EGGConvolved_phase_Filename = strcat(SubjectDataRoot,'EGGtimeseries',filesep,...
            'PhaseXvolume_Convolved_S_',sprintf('%.2d',subj_idx),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),...
            '_ord_',num2str(cfg.fOrder),'_tw_',sprintf('%.2d',cfg.transitionWidth));
        EGGConvolved_amplitude_Filename = strcat(SubjectDataRoot,'EGGtimeseries',filesep,...
            'AmplitudeXvolume_Convolved_S_',sprintf('%.2d',subj_idx),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),...
            '_ord_',num2str(cfg.fOrder),'_tw_',sprintf('%.2d',cfg.transitionWidth));
        
        
        %% BOLD
        
        BOLDTimeseriesFilename = strcat(SubjectDataRoot,'MRItimeseries',filesep,'fMRItimeseries_S',sprintf('%.2d',subj_idx),'_kw',num2str(cfg.kernelWidth),'.mat');
        
        BOLDBPAmplitudeTimeseriesFilename = strcat(SubjectDataRoot,'MRItimeseries',filesep,'BOLDFiltered_S',sprintf('%.2d',subj_idx),'_kw',...
            num2str(cfg.kernelWidth),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),'_fOrder_',num2str(cfg.fOrder),'_tw_',sprintf('%.2d',cfg.transitionWidth));
        BOLD_filtered_fullbandFilename = strcat(SubjectDataRoot,'MRItimeseries',filesep,'BOLDFULLBAND_S',sprintf('%.2d',subj_idx),'_kw',...
            num2str(cfg.kernelWidth));
        BOLD_filtered_narrowbandFilename = strcat(SubjectDataRoot,'MRItimeseries',filesep,'BOLDnarrowBAND_S',sprintf('%.2d',subj_idx),'_kw',...
            num2str(cfg.kernelWidth),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread));
        
        filename_csf_Signal_FB = strcat(SubjectDataRoot,'MRItimeseries',filesep,'filename_csfSignal_FB_S_',sprintf('%.2d',subj_idx),'_kw',...
            num2str(cfg.kernelWidth));
        filename_csfr_Residuals_FB =strcat(SubjectDataRoot,'MRItimeseries',filesep,'csfRegressionResiduals_FB_S_',sprintf('%.2d',subj_idx),'_kw',...
            num2str(cfg.kernelWidth));
        filename_csfr_Betas_FB = strcat(SubjectDataRoot,'MRItimeseries',filesep,'filename_csfRegressionBetas_FB_S_',sprintf('%.2d',subj_idx),'_kw',...
            num2str(cfg.kernelWidth));
        filename_csfr_Residuals_FB_phases = strcat(SubjectDataRoot,'MRItimeseries',filesep,'csfResiduals_FB_phases_s',sprintf('%.2d',subj_idx),'_kw',...
            num2str(cfg.kernelWidth),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread));
         filename_BOLDpoweratGastric = strcat(SubjectDataRoot,'MRItimeseries',filesep,'BOLDpoweratGastric_s',sprintf('%.2d',subj_idx),'_kw',...
            num2str(cfg.kernelWidth),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread));
        
        %% Phase analysis
        
        if cfg.scrubbed == 1
            PLVXVoxelFilename_csfr = strcat(SubjectDataRoot,'PhasesAnalysis',filesep,'PLVxVoxel_scrubbed_csfr_S_',sprintf('%.2d',subj_idx),...
                '_kw',num2str(cfg.kernelWidth),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),...
                '_fOrder_',num2str(cfg.fOrder),'_tw_',sprintf('%.2d',cfg.transitionWidth));
            medianRotationFilename_csfr = strcat(SubjectDataRoot,'PhasesAnalysis',filesep,'medianRotation_scrubbed_csfr_S_',sprintf('%.2d',subj_idx),...
                '_kw',num2str(cfg.kernelWidth),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),'_fOrder_',num2str(cfg.fOrder),'_tw_',sprintf('%.2d',cfg.transitionWidth));
            
        elseif cfg.hrfconvolved == 1
            PLVXVoxelFilename_csfr = strcat(SubjectDataRoot,'PhasesAnalysis',filesep,'PLVxVoxel_csfr_convolvedEGG_S_',sprintf('%.2d',subj_idx),...
                '_kw',num2str(cfg.kernelWidth),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),...
                '_fOrder_',num2str(cfg.fOrder),'_tw_',sprintf('%.2d',cfg.transitionWidth));
            medianRotationFilename_csfr = strcat(SubjectDataRoot,'PhasesAnalysis',filesep,'medianRotation_csfr_convolvedEGG_S_',sprintf('%.2d',subj_idx),...
                '_kw',num2str(cfg.kernelWidth),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),'_fOrder_',num2str(cfg.fOrder),'_tw_',sprintf('%.2d',cfg.transitionWidth));
            
        else % NORMAL CASE
            PLVXVoxelFilename_csfr = strcat(SubjectDataRoot,'PhasesAnalysis',filesep,'PLVxVoxel_csfr_S_',sprintf('%.2d',subj_idx),...
                '_kw',num2str(cfg.kernelWidth),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),...
                '_fOrder_',num2str(cfg.fOrder),'_tw_',sprintf('%.2d',cfg.transitionWidth));
            medianRotationFilename_csfr = strcat(SubjectDataRoot,'PhasesAnalysis',filesep,'medianRotation_csfr_S_',sprintf('%.2d',subj_idx),...
                '_kw',num2str(cfg.kernelWidth),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),'_fOrder_',num2str(cfg.fOrder),'_tw_',sprintf('%.2d',cfg.transitionWidth));
            
        end
        
        
        clusterTimeseries_filename = strcat(SubjectDataRoot,'PhasesAnalysis',filesep,'ClusterTimeseries_csfr_S_',sprintf('%.2d',subj_idx),...
            '_kw',num2str(cfg.kernelWidth),'_CA',sprintf('%.4d',cfg.clusterAlpha*10000));
        
        clusterTimeseries_spectrum_filename  = strcat(SubjectDataRoot,'PhasesAnalysis',filesep,'ClusterTimeseriesSpectrum_csfr_S_',sprintf('%.2d',subj_idx),...
            '_kw',num2str(cfg.kernelWidth),'_CA',sprintf('%.4d',cfg.clusterAlpha*10000));
        
        clusterTimeseries_coherence_filename  = strcat(SubjectDataRoot,'PhasesAnalysis',filesep,'ClusterTimeseriesCoherence_csfr_S_',sprintf('%.2d',subj_idx),...
            '_kw',num2str(cfg.kernelWidth),'_CA',sprintf('%.4d',cfg.clusterAlpha*10000));
        
        clusterTimeseries_phaseAngle_filename  = strcat(SubjectDataRoot,'PhasesAnalysis',filesep,'ClusterTimeseriesPhaseAngle_csfr_S_',sprintf('%.2d',subj_idx),...
            '_kw',num2str(cfg.kernelWidth),'_CA',sprintf('%.4d',cfg.clusterAlpha*10000));
        
        clusterTimeseries_insula_phaseAngle_filename  = strcat(SubjectDataRoot,'PhasesAnalysis',filesep,'ClusterTimeseriesInsulaPhaseAngle_csfr_S_',sprintf('%.2d',subj_idx),...
            '_kw',num2str(cfg.kernelWidth),'_CA',sprintf('%.4d',cfg.clusterAlpha*10000));
        
        voxelsTimeseries_phaseAngle_filename = strcat(SubjectDataRoot,'PhasesAnalysis',filesep,'VoxelTimeseriesPhaseAngle_csfr_S_',sprintf('%.2d',subj_idx),...
            '_kw',num2str(cfg.kernelWidth),'_CA',sprintf('%.4d',cfg.clusterAlpha*10000));
        
        AllRotationsFilename_csfr = strcat(SubjectDataRoot,'PhasesAnalysis',filesep,'AllRotationsPLV_csfr_S_',sprintf('%.2d',subj_idx),...
            '_kw',num2str(cfg.kernelWidth),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),'_fOrder_',num2str(cfg.fOrder),'_tw_',sprintf('%.2d',cfg.transitionWidth));
        
        couplingStrenghtFilename_csfr = strcat(SubjectDataRoot,'GlobalSignal',filesep,'couplingStrenght_csfr_S_',sprintf('%.2d',subj_idx),...
            '_kw',num2str(cfg.kernelWidth),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),'_fOrder_',num2str(cfg.fOrder),'_tw_',sprintf('%.2d',cfg.transitionWidth));
        
        PPC_filename = strcat(SubjectDataRoot,'PhasesAnalysis',filesep,'PPC_S_',sprintf('%.2d',subj_idx),...
            '_kw',num2str(cfg.kernelWidth),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),'_fOrder_',num2str(cfg.fOrder),'_tw_',sprintf('%.2d',cfg.transitionWidth));
        PPC_surrogate_filename = strcat(SubjectDataRoot,'PhasesAnalysis',filesep,'PPC_surrogate_S_',sprintf('%.2d',subj_idx),...
            '_kw',num2str(cfg.kernelWidth),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread),'_fOrder_',num2str(cfg.fOrder),'_tw_',sprintf('%.2d',cfg.transitionWidth));  
                
        % ECG
        data_ECG_cutted = strcat(SubjectDataRoot,'Heart',filesep,'ECG_noArtifact_cutted_S',sprintf('%.2d',subj_idx),'.mat');
        HRV_timeseries = strcat(SubjectDataRoot,'Heart',filesep,'IBIts_S',sprintf('%.2d',subj_idx),'.mat');
        HRV_BOLD_spectrum = strcat(SubjectDataRoot,'Heart',filesep,'HRV_BOLD_spectrum_S',sprintf('%.2d',subj_idx),'.mat');
        HRV_BOLD_coherence = strcat(SubjectDataRoot,'Heart',filesep,'HRV_BOLD_coherence_S',sprintf('%.2d',subj_idx),'.mat');
        HRV_BOLD_coherence_map = strcat(SubjectDataRoot,'Heart',filesep,'HRV_BOLD_coherence_map_S',sprintf('%.2d',subj_idx));
        
        
        %% Global Not used
        
        %         GlobalSignal_CSFr_FB_filename = strcat(SubjectDataRoot,'GlobalSignal',filesep,'GlobalSignal_CSFr_FB_S_',sprintf('%.2d',subj_idx),'_kw',...
        %             num2str(cfg.kernelWidth),'_fir2_fspread_',sprintf('%.3d',cfg.frequencySpread));
        %
        %
        %         GlobalSignal_CSFr_FB_betas_filename = strcat(SubjectDataRoot,'GlobalSignalRegressionBetas_FB_S_',sprintf('%.2d',subj_idx),'_kw',...
        %             num2str(cfg.kernelWidth));
        
        
        
        
        
        
        %% Other files
        filename_fMRImovementParameters = strcat(global_path2subject(subj_idx),'fMRI',filesep,...
            'acquisition1',filesep,'RestingState',filesep,'rp_afPHYSIENS_Sujet',sprintf('%.2d',subj_idx),...
            '-0002-00001-000001-01.txt');
        filename_InstantaneousMovementEstimates = strcat(SubjectDataRoot,'Other',filesep,...
            'InstantaneousMovementEstimates_S_',sprintf('%.2d',subj_idx));
        filename_BetasInstantaneousMovementEstimates = strcat(SubjectDataRoot,'Other',filesep,...
            'BetasInstantaneousMovementEstimates_S_',sprintf('%.2d',subj_idx));
        
        
        
    end
end


filename = eval(fileType);

end
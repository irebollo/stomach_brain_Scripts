function networks_couplingStrenght_obtain(subj_idx,cfgMain)
%{
Computes difference between median and empirical PLV (coupling strenght)
and stores in the HDD. These values are later used in the control with
head micromovements

Inputs:
empirical PLV x Voxel
    Y:\Subjects\Subject13\Timeseries\PhasesAnalysis\PLVxVoxel_csfr_S_13_kw3_fir2_fspread_015_fOrder_5_tw_15
chance PLV per voxel
    Y:\Subjects\Subject13\Timeseries\PhasesAnalysis\medianRotation_csfr_S_13_kw3_fir2_fspread_015_fOrder_5_tw_15

Outputs:
Coupling strenght x Voxel
    Y:\Subjects\Subject13\Timeseries\GlobalSignal\couplingStrenght_csfr_S_13_kw3_fir2_fspread_015_fOrder_5_tw_15


%}

%% Load files

couplingStrenghtFilename_csfr = global_filename(subj_idx,cfgMain,'couplingStrenghtFilename_csfr')
insideBrain = tools_getIndexBrain('inside');

PLVXVoxelFilename = global_filename(subj_idx,cfgMain,'PLVXVoxelFilename_csfr');
PLVXVoxel =ft_read_mri(strcat(PLVXVoxelFilename,'.nii'));
PLVXVoxel = PLVXVoxel.anatomy(:);
PLVXVoxel = PLVXVoxel(insideBrain);

medianRotationFilename = global_filename(subj_idx,cfgMain,'medianRotationFilename_csfr');
medianRotation =ft_read_mri(strcat(medianRotationFilename,'.nii'));
medianRotation = medianRotation.anatomy(:);
medianRotation = medianRotation(insideBrain);

couplingstrenght = PLVXVoxel - medianRotation;

%% Plots

plotDir = strcat (global_path2subject(subj_idx),'PreprocessingLog',filesep);
plotFilename = strcat(plotDir,'S_',sprintf('%.2d',subj_idx),'_CouplingStrenght');

if cfgMain.savePlots == 1
    
    if cfgMain.plotFigures == 0;
        SanityPlot = figure('visible','off');
    else
        SanityPlot = figure('visible','on');
    end
   
    voxelCoordinates = sub2ind([53,63,46],11,30,37); % Right somatosensory cortex
    voxelCoordinates_inside = zeros(153594,1);
    voxelCoordinates_inside(voxelCoordinates)=1;
    voxelCoordinates_inside = voxelCoordinates_inside(insideBrain);
    ind_voxelCoordinates_inside = find(voxelCoordinates_inside);
    
    nhist(couplingstrenght)
    xlabel('Coupling strenght')
    title(['S',sprintf('%.2d',subj_idx),32,'Coupling across bain. Mean:' num2str(mean(couplingstrenght)) ' rSS voxel:' 32 num2str(couplingstrenght(ind_voxelCoordinates_inside))],'fontsize',18)
    
    
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    set(gcf, 'PaperPositionMode', 'auto');
    
    print ('-dpng', '-painters', eval('plotFilename'))
    print ('-depsc2', '-painters', eval('plotFilename'))
    saveas(SanityPlot,strcat(plotFilename,'.fig'))
    
end
%% save

save(couplingStrenghtFilename_csfr,'couplingstrenght')

end
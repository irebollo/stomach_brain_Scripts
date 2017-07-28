function timeseries_FWD_BOLD_GLM(subj_idx,cfgMain)

%{
Store the timeseries of absolute derivative of movement (separately for
rotation and translations) i.e. Framewise displacement. 
Y:\Subjects\Subject13\Timeseries\Other\InstantaneousMovementEstimates_S_13

Then compute the
GLM between the BOLD timeseries and the framewise displacement and stored in:
Y:\Subjects\Subject13\Timeseries\Other\BetasInstantaneousMovementEstimates_S_13

To be later used by the script paper_control_corrMovement_CS.


Store figure into disk for each subject.
IR 04/07/2017
%}

%% Input and output filenames 


plotDir = strcat (global_path2subject(subj_idx),'PreprocessingLog',filesep);
plotFilename = strcat(plotDir,'S_',sprintf('%.2d',subj_idx),'_MovementEffectINGlobalSignal');

filename_InstantaneousMovementEstimates = global_filename(subj_idx,cfgMain,'filename_InstantaneousMovementEstimates');
filename_BetasInstantaneousMovementEstimates = global_filename(subj_idx,cfgMain,'filename_BetasInstantaneousMovementEstimates');

% Load mask brain
insideBrain = tools_getIndexBrain('inside');
gray_mask_filename = ['D:' filesep 'Physiens' filesep 'scripts4paper' filesep 'files' filesep 'gray_mask.hdr']
grayMatter = ft_read_mri(gray_mask_filename);
grayMatter = logical(grayMatter.anatomy);
grayMatter= grayMatter(:);
grayMatter_insideBrain = grayMatter(insideBrain);

% Load EGG
filenameEGG = global_filename(subj_idx,cfgMain,'EGGAmplitudeXVolumeFilename');
load(filenameEGG)

% load gs
GlobalSignal_CSFr_FB_filename = global_filename(subj_idx,cfgMain,'GlobalSignal_CSFr_FB_filename');
load(GlobalSignal_CSFr_FB_filename)

% load BOLD timeseries
BOLDTimeseriesFilename = global_filename(subj_idx,cfgMain,strcat('filename_',cfgMain.Timeseries2Regress,'_Residuals_FB'));
timeseries = load(BOLDTimeseriesFilename);

%% Do
% filter
peaksAllsubjects = global_getEGGpeaks;
indPeak = find (peaksAllsubjects(:,1) == subj_idx);
mostPowerfullFrequency = peaksAllsubjects(indPeak,3);
centerFrequency = mostPowerfullFrequency; %
filter_frequency_spread=cfgMain.frequencySpread/1000; % In hz
sr = 0.5; % 1 TR = 2s
filterOrder=(cfgMain.fOrder*fix(sr/(centerFrequency-filter_frequency_spread))-1);%in nsamples
transition_width= cfgMain.transitionWidth/100; % in normalised units
filteredMRI=tools_bpFilter(eval(strcat('timeseries.',cell2mat(fields(timeseries)))),sr,filterOrder,centerFrequency,filter_frequency_spread,transition_width,cfgMain.filterType);
filtered_GS = nanmean(filteredMRI(:,grayMatter_insideBrain),2);

mvp = global_filename(subj_idx,cfgMain,'filename_fMRImovementParameters');
Movement = load(mvp)

Movement_sum_abs_translation_Instantaneous = zeros(450,1);
Movement_sum_abs_rotations_Instantaneous = zeros(450,1);
for i=2:length(Movement)
   Movement_sum_abs_translation_Instantaneous(i,:) = sum(abs(Movement(i-1,1:3)' - Movement(i,1:3)'));
   Movement_sum_abs_rotations_Instantaneous(i,:)   = sum(abs(Movement(i-1,4:6)' - Movement(i,4:6)')); % In radians

end

% multiply by fifty (mean radius of the brain) to pass from radians to militers
framewiseDisplacement = [(Movement_sum_abs_rotations_Instantaneous).*50 + Movement_sum_abs_translation_Instantaneous]


% figure
% plot(framewiseDisplacement)

save(filename_InstantaneousMovementEstimates,'framewiseDisplacement','Movement_sum_abs_translation_Instantaneous','Movement_sum_abs_rotations_Instantaneous')

%% GLM; bold timeseries explained by framewise displacement 

toBeExplained = timeseries.error_csf_z(cfgMain.beginCut:cfgMain.endCut,:); % BOLD timeseries will be the variable to be explained out in the GLM
explainingVariables = []
explainingVariables = zscore(framewiseDisplacement(cfgMain.beginCut:cfgMain.endCut))

% figure
% plot(explainingVariables)
% plot(zscore(Movement_sum_abs_translation_Instantaneous),'r')
% plot(zscore(Movement_sum_abs_rotations_Instantaneous),'g')

betas_instantaneousMovement = EfficientGLM(toBeExplained,explainingVariables); % Obtain the betas indicating how much the predicting variable predicts the data
% 
% figure
% plot(betas_instantaneousMovement,couplingstrenght,'ok')

save(filename_BetasInstantaneousMovementEstimates,'betas_instantaneousMovement')
% this file is later used to correlate with copupling strenght

% subplot(5,1,3)
% imagesc(abs(predictedBOLD(cfgMain.beginCut:cfgMain.endCut,grayMatter_insideBrain)'))

%% Make figure

if cfgMain.savePlots == 1
       
    if cfgMain.plotFigures == 0;
        SanityPlot = figure('visible','off');
    else
        SanityPlot = figure('visible','on');
    end
   
    
subplot(4,1,1)
plot(zscore(EGGTimeseries),'LineWidth',3)
hold on
% plot(gs_timecourse,'r','LineWidth',3)
plot(zscore(filtered_GS(cfgMain.beginCut:cfgMain.endCut)),'r--','LineWidth',3)

xlim([0 420])
legend('EGG','GS')
title(['subject ' num2str(subj_idx)],'FontSize',20)


subplot(4,1,2)
plot(framewiseDisplacement(cfgMain.beginCut:cfgMain.endCut),'g','LineWidth',3)
xlim([0 420])
title('Movement','FontSize',20)
ylim([0 1])

subplot(4,1,3)

imagesc(timeseries.error_csf_z(cfgMain.beginCut:cfgMain.endCut,grayMatter_insideBrain)')
title('Unfiltered voxels','FontSize',20)


subplot(4,1,4)
imagesc(filteredMRI(cfgMain.beginCut:cfgMain.endCut,grayMatter_insideBrain)')
title('filtered voxels','FontSize',20)
    

% subplot(5,1,5)
% imagesc(abs(predictedBOLD(cfgMain.beginCut:cfgMain.endCut,grayMatter_insideBrain)'))


    colormap(hot)

    
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    set(gcf, 'PaperPositionMode', 'auto');
    
    print ('-dpng', '-painters', eval('plotFilename'))
    print ('-depsc2', '-painters', eval('plotFilename'))
    saveas(SanityPlot,strcat(plotFilename,'.fig'))
    
end

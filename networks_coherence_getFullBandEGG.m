function networks_coherence_getFullBandEGG(subj_idx,cfgMain)
%{
Compute and stores the fullband (0.01-0.1 hz) EGG time series in order to compute
coherence with the BOLD signal in significant clusters

Input: Raw brainamp files stored in :
Y:\Subjects\Subject13\brainamp\with MRI

Output:
Y:\Subjects\Subject13\Timeseries\EGGtimeseries\EGGTimeseriesFullband_S_13



IR 28/06/2017
%}
%% parameters inputs and outputs
 
FilenameamplitudeXVolumeBestChannel_FULLBAND = global_filename(subj_idx,cfgMain,'FilenameamplitudeXVolumeBestChannel_FULLBAND')

 peaksAllsubjects = global_getEGGpeaks;
 indPeak = find (peaksAllsubjects(:,1) == subj_idx);
 bestChannel = peaksAllsubjects(indPeak,2);
    

%% Load raw data

[EGG_raw,markers_raw] = prepro_egg_loadData(subj_idx,'EGG',1); % Only loads those channel with label 'EGG'

%% Downsample data

disp('Resampling...')
cfg = [];  %configuration structure
cfg.detrend = 'no'; % remove linear trend from the data (done per trial)
cfg.demean = 'yes';
cfg.resamplefs= 0.5; % 4 x top-freq (15 cpm = 0.25 Hz) - Nyquist = 30 cpm  frequency at which the data will be resampled
EGG_downsampled = ft_resampledata(cfg,EGG_raw); % This supposed to filter the data at half the new sr but we do not know where it does it

% Downsample markers

markersDownsampled = prepro_egg_downsampleMarkers(EGG_raw,EGG_downsampled,markers_raw); % Get closest sample number of new marker based on the timestamps of the not downsampled data

%% Get begginign and end of fMRI adquisition and cut file

beginIRM = markersDownsampled(1,1); % 450 volumes, 900 seconds, end IRM is at marker 450
endIRM = markersDownsampled(450,1);
cfg_resample = [];
cfg_resample.begsample = (beginIRM)*EGG_downsampled.fsample;
cfg_resample.endsample = (endIRM)*EGG_downsampled.fsample;
EGG_downsampled_cut = ft_redefinetrial(cfg_resample,EGG_downsampled);
EGG_downsampled_cut.time{1,1} = EGG_downsampled_cut.time{1,1} - EGG_downsampled_cut.time{1,1}(1);% time starts at 0

%% Preprocess and filter EGG
EGG_FullBand = ft_preproc_polyremoval (EGG_downsampled_cut.trial{1,1}(bestChannel,:), 2);
EGG_FullBand = ft_preproc_bandpassfilter (EGG_FullBand,0.5,[0.01 0.1]);

EGG_FullBand = EGG_FullBand(cfgMain.beginCut:cfgMain.endCut)
EGG_FullBand = ft_preproc_standardize (EGG_FullBand);

%% Plots

plotDir = strcat (global_path2subject(subj_idx),'PreprocessingLog',filesep);
plotFilename = strcat(plotDir,'S_',sprintf('%.2d',subj_idx),'_EGGFULLBAND');

if cfgMain.savePlots == 1
    
    if cfgMain.plotFigures == 0;
        SanityPlot = figure('visible','off');
    else
        SanityPlot = figure('visible','on');
    end
  
    plot(EGG_FullBand)
    xlabel('time')
    ylabel('Z')
    title(['S',sprintf('%.2d',subj_idx),32,'Filtered 0.1 - 1 Hz detrended'],'fontsize',18)
    
    
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    set(gcf, 'PaperPositionMode', 'auto');
    
    print ('-dpng', '-painters', eval('plotFilename'))
    print ('-depsc2', '-painters', eval('plotFilename'))
    saveas(SanityPlot,strcat(plotFilename,'.fig'))
    
end

%% save

save(FilenameamplitudeXVolumeBestChannel_FULLBAND,'EGG_FullBand')
end


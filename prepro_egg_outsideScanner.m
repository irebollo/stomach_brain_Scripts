function prepro_egg_outsideScanner(subj_idx,cfgEGG)

% Similar to prepro_egg buy preprocess the brainamp files containing the EGG recorded outside the
% scanner,

% cfgEGG.automaticChannelSelection: THIS SHOULD BE TURNED TO 1 so that it
% loads the same channel and frequency as for the 

% INPUT:
% subj_idx= number of subject
% cfgEGG must contain
% cfgEGG.fOrder FILTER ORDER
% cfgEGG.frequencySpread % HWHM of filter
% cfgEGG.plotFigures, if 1 will plot figures and store them on disk (where??)


% Output: amplitude X volume and phase per volume, stored separatly in each
% subject's tiemseires folder, together with best channel and most power
% frequency
%
% This Function first loads each subject raw brainamp data and marker files and preprocess them,
% including downsampling to 10 hz (also the markers), removoing the segments
% before and after the MRI acqusition. Next step involves analyzing the
% power spectrum of the data. For this we used welch method in  fieldtrip,
% which implies estimating the power in olveraping trials of 200 seconds and then averaging the power over trials
% This script uses the channel and power info used

% IR 14/09/2016

%% Pass cfg parameters to function variables

fOrder = cfgEGG.fOrder;
frequencySpread = cfgEGG.frequencySpread;
plotFigures = cfgEGG.plotFigures;

beginCut =cfgEGG.beginCut;
endCut = cfgEGG.endCut;
automaticChannelSelection = cfgEGG.automaticChannelSelection;
offset = cfgEGG.offset;

logEGGpreprocessing = []; % initialize log structure

%% Define paths, input and output filenames


% Filenames of Plots and logs
plotDir = strcat (global_path2subject(subj_idx),'PreprocessingLog',filesep);
downsampledAllChannelsFilename = strcat(plotDir,'S_',sprintf('%.2d',subj_idx),'timeseriesAllchannelsDownsampled_OutsideScanner');
spectrumFilename = strcat(plotDir,'S_',sprintf('%.2d',subj_idx),'_spectrum_OutsideScanner');
filteredPlotFilename= strcat(plotDir,'S_',sprintf('%.2d',subj_idx),'timeseriesBestchannel_OutsideScanner');


% Output
    EGGoutsideScannerFilename = global_filename(subj_idx,cfgEGG,'EGGoutsideScannerFilename');

%% Load raw data

[EGG_raw,markers_raw] = prepro_egg_loadData(subj_idx,'EGG',0); % Only loads those channel with label 'EGG'

%% Downsample data

disp('Resampling...')
cfg = [];  %configuration structure
cfg.detrend = 'no'; % remove linear trend from the data (done per trial)
cfg.demean = 'yes';
cfg.resamplefs= 10; % 4 x top-freq (15 cpm = 0.25 Hz) - Nyquist = 30 cpm  frequency at which the data will be resampled
EGG_downsampled = ft_resampledata(cfg,EGG_raw); % This supposed to filter the data at half the new sr but we do not know where it does it

%% Plot downsampled
    figureDownsampled = figure('visible','off');

for iChannel = 1:4
    subplot(4,1,iChannel)
    plot(EGG_downsampled.time{1,1}(1,:),EGG_downsampled.trial{1,1}(iChannel,:))
    title(strcat(('Outside scaner Downsampled EGG n°'),num2str(iChannel),32,'participant',sprintf('%.2d',subj_idx)), 'fontsize',11);
    xlabel('time in s')
    ylabel('amplitude')
    set(gca,'fontsize',12)
end

set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gcf, 'PaperPositionMode', 'auto');

print ('-dpng', '-painters', eval('downsampledAllChannelsFilename'))

print ('-depsc2', '-painters', eval('downsampledAllChannelsFilename'))
saveas(figureDownsampled,strcat(downsampledAllChannelsFilename,'.fig'))

%% Welch's spectrum

% outside
len = EGG_downsampled.fsample*200; % length of subtrials of 200s in samples
EGG_downsampled.sampleinfo=[1 max(EGG_downsampled.time{1,1})*EGG_downsampled.fsample];

cfg = [];
cfg.trl(:,1) = EGG_downsampled.sampleinfo(1):(len/4):EGG_downsampled.sampleinfo(2)-len+1;%trial start in samples from begining of raw data
cfg.trl(:,2) = EGG_downsampled.sampleinfo(1)+len-1:(len/4):EGG_downsampled.sampleinfo(2);%trial ends in samples from begining of raw data
cfg.trl(:,3) = 0; %offset of the trigger with respect to the trial

EGG_trials = ft_redefinetrial(cfg,EGG_downsampled);

%Now we will run a Hann tapered fft on each of the trials, and the resulting power spectra will be averaged for us giving a smooth power output.
% fft
cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.output = 'pow';
cfg.pad = 1000;
cfg.foilim = [0 0.1]; % 0 - 6 cpm
frequencyWelch = ft_freqanalysis(cfg,EGG_trials);


% Retrieve information about best channel from the disk

    peaksAllsubjects = global_getEGGpeaks;
    indPeak = find (peaksAllsubjects(:,1) == subj_idx);
    bestChannel = peaksAllsubjects(indPeak,2);

           
    indexFrequencies = find (frequencyWelch.freq >= 0.03333 & frequencyWelch.freq <= 0.06666);
    
    %Automatically get the channel with highest peak
    maxPowerXChannel = zeros(4,2); % column 1 max power, column 2 frequency location
    for iChannel=1:4
        maxPowerXChannel(iChannel,1) = max(frequencyWelch.powspctrm(iChannel,indexFrequencies));% from 0.03 to 0.07 hz %0.04 o ).7
        maxPowerLocation = frequencyWelch.powspctrm(iChannel,:)==maxPowerXChannel(iChannel,1);
        maxPowerXChannel(iChannel,2) = frequencyWelch.freq(find(maxPowerLocation));
    
    end
    
   
    mostPowerfullFrequency = maxPowerXChannel(bestChannel,2);

    filterWidth= frequencySpread/1000;% Attention for visualization in the spectrum domain only, filter shape is different
    
    filterPlot=zeros(4,100);
    for iChannel=1:4
        indexFrequenciesFilter = find (frequencyWelch.freq >= maxPowerXChannel(iChannel,2)-filterWidth & frequencyWelch.freq <= maxPowerXChannel(iChannel,2)+filterWidth);
        filterPlot(iChannel,indexFrequenciesFilter)= maxPowerXChannel(iChannel,1);
    end
    
    indexFrequenciesPloting = find (frequencyWelch.freq >= 0.0189 & frequencyWelch.freq <= 0.0698);
    % Which frequency range is going to appear in the plot
    
        

        figureSpectrum = figure('visible','off');
        for iChannel=1:4
            subplot(2,2,iChannel);
            plot(frequencyWelch.freq(indexFrequenciesPloting),frequencyWelch.powspctrm(iChannel,indexFrequenciesPloting),'-o','lineWidth',3);
% %             if iChannel == mostPowerfullChannel
                title(strcat(('Welch outside scaner EGG n°'),num2str(iChannel),32,'freq',num2str(maxPowerXChannel(iChannel,2)),32,'participant',sprintf('%.2d',subj_idx)),'fontweight','bold', 'fontsize',11);
            hold on;
            plot (frequencyWelch.freq(indexFrequenciesPloting),filterPlot(iChannel,indexFrequenciesPloting),'r','lineWidth',3)
            set(gca,'fontsize',11)
            
            set(gcf,'units','normalized','outerposition',[0 0 1 1])
            set(gcf, 'PaperPositionMode', 'auto');
            
        end
        print ('-depsc2', '-painters', eval('spectrumFilename'))
        print ('-dpng', '-painters', eval('spectrumFilename'))

        saveas(figureSpectrum,strcat(spectrumFilename,'.fig'))

%% Fill log file
% 
logEGGpreprocessing.subjectNumber = subj_idx;
logEGGpreprocessing.cfgEGG = cfgEGG;
logEGGpreprocessing.bestChannel = bestChannel;
logEGGpreprocessing.mostPowerfullFrequency = mostPowerfullFrequency;

%% Save

save(EGGoutsideScannerFilename,'logEGGpreprocessing')


end
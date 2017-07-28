function prepro_egg(subj_idx,cfgMain)

%{
Basic preprocessing steps for obtaining and storing the EGG phase per mri
volume

INPUT to function:
subj_idx= number of subject

cfgMain with the following PARAMETERS: 
cfg.fOrder filter order multiplicative factor
cfg.frequencySpread % HWHM of filter
cfg.plotFigures, if 1 will plot figures and store them on preprocessingLog folder of each subject
cfg.automaticChannelSelection: if 1 uses the information stored in the disk of which frequency and channel should be
used for each subject, if 0 user have to choose manually which is the best channel and then call prepro_egg_saveChannelInfo to store this information
on the disk

cfgMain.beginCut;cfgMain.endCut; at which volume data will be cuted to get rid of the fMRI artifact

Input file: Y:\Subjects\Subject13\brainamp\with MRI\PHYSIENS_2014-10-170001.eeg
Output: mean amplitude and phase per volume, stored separatly in each
subject's timeseries folder together with a log file containing
information about the preprocessing: frequency and channel used

Y:\Subjects\Subject13\Timeseries\EGGtimeseries\PhaseXvolume_S_13_fir2_fspread_015_ord_5_tw_15

This Function first loads each subject raw brainamp data and marker files and preprocess them,
including downsampling to 10 hz (also the markers), removoing the segments
before and after the MRI acqusition. Next step involves analyzing the
power spectrum of the data. For this we used welch method in  fieldtrip,
which implies estimating the power in olveraping trials of 200 seconds and then averaging the power over trials
The scripts ask to visually identify the channel with more power and
normogastria and this info is stored at the disk in the log file. The
signal is filtered at this frequency and channel, then resampled to the
SR of the MRI (0.5 Hz), then the hilbert transform is applied and the
resulting value of phase angle per volume is stored into disk

IR 28/07/2015
IR 12/09/2016 Fixed an error while storing the mostpowerfull frequency
and channel
IR 13/09/2016 Added sanity check, store in the HDD a plot with the phases
and amplitude timeseries of the chosen channel and the unfiltered power spectra of that channel,
it thens performs the fourier transform of the filtered timeseries and
plots them

%}
%% Pass cfg parameters to function variables

fOrder = cfgMain.fOrder;
frequencySpread = cfgMain.frequencySpread;
plotFigures = cfgMain.plotFigures;
beginCut =cfgMain.beginCut;
endCut = cfgMain.endCut;
automaticChannelSelection = cfgMain.automaticChannelSelection;
offset = cfgMain.offset;

logEGGpreprocessing = []; % initialize log structure

%% Define paths, input and output filenames

% Filenames of Plots and logs
plotDir = strcat (global_path2subject(subj_idx),'PreprocessingLog',filesep);
downsampledAllChannelsFilename = strcat(plotDir,'S_',sprintf('%.2d',subj_idx),'timeseriesAllchannelsDownsampled_InsideScanner');
spectrumFilename = strcat(plotDir,'S_',sprintf('%.2d',subj_idx),'_spectrum_InsideScanner');
filteredPlotFilename= strcat(plotDir,'S_',sprintf('%.2d',subj_idx),'timeseriesBestchannel_InsideScanner');

% Output
EGGPhaseXVolumeFilename = global_filename(subj_idx,cfgMain,'EGGPhaseXVolumeFilename');
EGGAmplitudeXVolumeFilename = global_filename(subj_idx,cfgMain,'EGGAmplitudeXVolumeFilename');

%% Load raw data

[EGG_raw,markers_raw] = prepro_egg_loadData(subj_idx,'EGG',1); % Only loads those channel with label 'EGG'

%% Downsample data

disp('Resampling...')
cfg = [];  %initialize configuration structure
cfg.detrend = 'no'; % remove linear trend from the data (done per trial)
cfg.demean = 'yes';
cfg.resamplefs= 10; % 4 x top-freq (15 cpm = 0.25 Hz) - Nyquist = 30 cpm  frequency at which the data will be resampled
EGG_downsampled = ft_resampledata(cfg,EGG_raw); % This procedure also lowpass filter the data at half the new sr 

% Downsample markers

markersDownsampled = prepro_egg_downsampleMarkers(EGG_raw,EGG_downsampled,markers_raw); 
% Get closest sample number of new marker based on the timestamps of the not downsampled data

%% Get begginign and end of fMRI adquisition and cut file for calculating the spectrum

beginIRM = markersDownsampled(1,1); % 450 volumes, 900 seconds, end IRM is at marker 450
endIRM = markersDownsampled(450,1);
cfg_resample = [];
cfg_resample.begsample = (beginIRM)*EGG_downsampled.fsample;
cfg_resample.endsample = (endIRM)*EGG_downsampled.fsample;
EGG_downsampled_cut = ft_redefinetrial(cfg_resample,EGG_downsampled);
EGG_downsampled_cut.time{1,1} = EGG_downsampled_cut.time{1,1} - EGG_downsampled_cut.time{1,1}(1);% time starts at 0

%% Plot downsampled cut file to visually check for artifacts

if plotFigures == 0;
    figureDownsampled = figure('visible','off');
else
    figureDownsampled = figure('visible','on');
end

for iChannel = 1:4
    subplot(4,1,iChannel)
    plot(EGG_downsampled_cut.time{1,1}(1,:),EGG_downsampled_cut.trial{1,1}(iChannel,:))
    title(strcat(('Downsampled EGG n°'),num2str(iChannel),32,'participant',sprintf('%.2d',subj_idx)), 'fontsize',11);
    xlabel('time in s')
    ylabel('amplitude')
    set(gca,'fontsize',12)
end

set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gcf, 'PaperPositionMode', 'auto');

if cfgMain.savePlots == 1
    print ('-dpng', '-painters', eval('downsampledAllChannelsFilename'))
    print ('-depsc2', '-painters', eval('downsampledAllChannelsFilename'))
    saveas(figureDownsampled,strcat(downsampledAllChannelsFilename,'.fig'))
end

%% Welch's spectrum
% Whelch spectrum is calculated with data including signal before and after
% the IRM recording (for visualization)


% The 900 seconds segment is divided into trials of 200s with 1/4 of
% overlap

len = EGG_downsampled_cut.fsample*200; % length of subtrials of 200s in samples

EGG_downsampled_cut.sampleinfo=[1 max(EGG_downsampled_cut.time{1,1})*EGG_downsampled_cut.fsample]; 

cfg = [];
cfg.trl(:,1) = EGG_downsampled_cut.sampleinfo(1):(len/4):EGG_downsampled_cut.sampleinfo(2)-len+1;%trial start in samples from begining of raw data
cfg.trl(:,2) = EGG_downsampled_cut.sampleinfo(1)+len-1:(len/4):EGG_downsampled_cut.sampleinfo(2);%trial ends in samples from begining of raw data
cfg.trl(:,3) = 0; %offset of the trigger with respect to the trial

EGG_trials = ft_redefinetrial(cfg,EGG_downsampled_cut);

%Now we will run a Hann tapered fft on each of the trials, and the resulting power spectra will
% be averaged for us giving a smooth power output
% fft
cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.output = 'pow';
cfg.pad = 1000; % zero padding
cfg.foilim = [0 0.1]; % 0 - 10 cpm
frequencyWelch = ft_freqanalysis(cfg,EGG_trials);


%% Get channel with more power

% Check if automatic channel selection is on, if it is, load chanel and frequency info
if cfgMain.automaticChannelSelection == 1
    peaksAllsubjects = global_getEGGpeaks;
    indPeak = find (peaksAllsubjects(:,1) == subj_idx);
    bestChannel = peaksAllsubjects(indPeak,2);
    
    mostPowerfullFrequency = peaksAllsubjects(indPeak,3);
    mostPowerfullChannel = 0;
    
else % choose manually the channel to use
    filterWidth= frequencySpread/1000;% Attention for visualization in the spectrum domain only, filter shape is different
    
    %Welch
    % Search for the largest peak in frquencies from 2 to 4 cycles per
    % second (normogastria)
    indexFrequencies = find (frequencyWelch.freq >= 0.03333 & frequencyWelch.freq <= 0.06666); % normogastria = 2-4 cpm
   
    %Automatically get the channel with highest peak
    maxPowerXChannel = zeros(4,2); % column 1 max power, column 2 frequency location
    for iChannel=1:4
        maxPowerXChannel(iChannel,1) = max(frequencyWelch.powspctrm(iChannel,indexFrequencies));% from 0.033 to 0.066 hz
        maxPowerLocation = frequencyWelch.powspctrm(iChannel,:)==maxPowerXChannel(iChannel,1);
        maxPowerXChannel(iChannel,2) = frequencyWelch.freq(find(maxPowerLocation));
    end
    
    [highestPower, mostPowerfullChannel] = max(maxPowerXChannel(:,1));
    mostPowerfullFrequency = maxPowerXChannel(mostPowerfullChannel,2);
    
    % get index where filter should be ploted
    filterPlot=zeros(4,100);
    for iChannel=1:4
        indexFrequenciesFilter = find (frequencyWelch.freq >= maxPowerXChannel(iChannel,2)-filterWidth & frequencyWelch.freq <= maxPowerXChannel(iChannel,2)+filterWidth);
        filterPlot(iChannel,indexFrequenciesFilter)= maxPowerXChannel(iChannel,1);
    end
    
    % Which frequency range is going to appear in the plot
    indexFrequenciesPloting = find (frequencyWelch.freq >= 0.0189 & frequencyWelch.freq <= 0.0698);
    
    
    if plotFigures == 0;
        figureSpectrum = figure('visible','off');
    else
        figureSpectrum = figure('visible','on');
             
        for iChannel=1:4
            subplot(2,2,iChannel);
            plot(frequencyWelch.freq(indexFrequenciesPloting),frequencyWelch.powspctrm(iChannel,indexFrequenciesPloting),'-o','lineWidth',3);
            if iChannel == mostPowerfullChannel
                title(strcat(('Welch EGG n°'),num2str(iChannel),32,'participant',sprintf('%.2d',subj_idx)),'fontweight','bold', 'fontsize',11);
            else
                title(strcat(('Welch EGG n°'),num2str(iChannel),32,'participant',sprintf('%.2d',subj_idx)), 'fontsize',11);
                
            end
            hold on;
            plot (frequencyWelch.freq(indexFrequenciesPloting),filterPlot(iChannel,indexFrequenciesPloting),'r','lineWidth',3)
            set(gca,'fontsize',11)
            
            set(gcf,'units','normalized','outerposition',[0 0 1 1])
            set(gcf, 'PaperPositionMode', 'auto');
            
        end
        
        if cfgMain.savePlots == 1
            print ('-dpng', '-painters', eval('spectrumFilename'))
            print ('-depsc2', '-painters', eval('spectrumFilename'))
            saveas(figureSpectrum,strcat(spectrumFilename,'.fig'))
        end
        
    end %plotFigures
    
    message1=strcat(' Most powerfull channel is channel',32,num2str(mostPowerfullChannel));
    
%% Manually select witch channel to use

    fprintf(message1)
    
    bestChannel = str2double(input('\nEnter the number of the channel you want to use:\n' ,'s'));
    mostPowerfullFrequency = maxPowerXChannel(bestChannel,2);
    
    message2=strcat(' Most powerfull frequency in this channel is ',32,num2str(mostPowerfullFrequency),'do you want to use that frequency?');
    fprintf(message2)
    
    useAutomaticFrequency = str2double(input('\n if yes enter 1,if Not put the number of the frequency you want to use:\n' ,'s'));
    if useAutomaticFrequency ~=1
        mostPowerfullFrequency = useAutomaticFrequency;
        
        
        message3=strcat(' The  power is this channel is ',32,num2str(maxPowerXChannel(bestChannel,1))...
            ,32,' and the peak frequency is',32, num2str(mostPowerfullFrequency) );
        fprintf('\n')
        fprintf(message3)
        
        
        presstocontinue = str2double(input('\nPress any key to continue:\n' ,'s'));
        
    end % not automatic channel selection
end 
    
    %% Filter at EGG peak frequency
    
    data=cell2mat(EGG_downsampled.trial);
    data=data(bestChannel,:);
    srate=EGG_downsampled.fsample;
    center_frequency=mostPowerfullFrequency + (offset/1000) ;
    filter_frequency_spread=frequencySpread / 1000;
    % HWHM of the filter, in order to be able to use the
    % parameter in the filename, the input to the function is 1000 bigger than
    % the actual frequency spead. e.g HWHM 0.005 Hz = frequencySpread=5
    
    lowerFilterBound = center_frequency - filter_frequency_spread; % When filtering starts
    filterOrder=fOrder*fix(srate/lowerFilterBound);%in nsamples
    transition_width= cfgMain.transitionWidth/100; % in normalised units
    [datapoints_EGG_ds_bandpass]= ...
        tools_bpFilter(data,srate,filterOrder,center_frequency,filter_frequency_spread,transition_width,cfgMain.filterType); % filtering
    EGG_ds_bandpass = EGG_downsampled; % Copy structure
    EGG_ds_bandpass.trial = mat2cell(datapoints_EGG_ds_bandpass); % update actual timeseries
    
    %% First downsample, then hilbert, then get mean phase per volume
    
    % Downsample to MRI sampling rate
    
    disp('Resampling...')
    cfg = [];  %configuration structure
    cfg.detrend = 'no'; % remove linear trend from the data (done per trial)
    cfg.demean = 'yes';
    cfg.resamplefs= 0.5; %one volume every two seconds
    EGG_ds_bp_downsampled = ft_resampledata(cfg,EGG_ds_bandpass);
    
    % downsample markers
    markers_DS_05Hz = prepro_egg_downsampleMarkers(EGG_raw,EGG_ds_bp_downsampled,markers_raw); % Downsample markwers to new SR
    
    % Hilbert transform
        phaseEGG = hilbert (EGG_ds_bp_downsampled.trial{1,1});
    
    % Get average phase value per volume
    phaseXVolume = zeros(1,450);
    nVolumes=450;
    for iTrial=1:nVolumes
        phaseXVolume(1,iTrial) = ...
            mean(phaseEGG(1,markers_DS_05Hz(iTrial,3):markers_DS_05Hz(iTrial,4)));
    end
    
    % Control, check that averaging phase is working
    % figure
    % plot(angle(phaseEGG(markers_DS_05Hz(1,3):end)),'-o')
    % hold on
    % plot(angle(phaseXVolume),'-or')
    
    
    phaseXVolume = phaseXVolume(beginCut:endCut); % cut begining and end of IRM acquisition
    EGGTimeseries = EGG_ds_bp_downsampled.trial{1,1};
    EGGTimeseries = EGGTimeseries(beginCut:endCut);
    
    %% Fill log file
    
    logEGGpreprocessing.subjectNumber = subj_idx;
    logEGGpreprocessing.cfgMain = cfgMain;
    logEGGpreprocessing.mostPowerfullChannel = mostPowerfullChannel;
    logEGGpreprocessing.bestChannel = bestChannel;
    logEGGpreprocessing.outputFilename = EGGPhaseXVolumeFilename;
    logEGGpreprocessing.automaticChannelSelection = cfgMain.automaticChannelSelection;
    logEGGpreprocessing.mostPowerfullFrequency = mostPowerfullFrequency;
%     logEGGpreprocessing.powerinChosenChannel = maxPowerXChannel(bestChannel,1);
    
    %% Save
    
    save(EGGPhaseXVolumeFilename,'phaseXVolume','logEGGpreprocessing')
    save(EGGAmplitudeXVolumeFilename,'EGGTimeseries','logEGGpreprocessing')
    save(strcat(EGGPhaseXVolumeFilename,'_log'),'logEGGpreprocessing')
    %% Plots and sanity check
    
    
    if plotFigures == 0;
        SanityPlot = figure('visible','off');
    else
        SanityPlot = figure('visible','on');
    end
    
    % Spectrum best channel unfltered timeseries
    subplot(4,1,1)
    indexFrequenciesPloting2 = find (frequencyWelch.freq >= 0.035 & frequencyWelch.freq <= 0.065);
    
    plot(frequencyWelch.freq(indexFrequenciesPloting2),frequencyWelch.powspctrm(bestChannel,indexFrequenciesPloting2),'-o','lineWidth',3);
    
    title(strcat(('Welch EGG n°'),num2str(bestChannel),32,'participant',sprintf('%.2d',subj_idx),32,'frequency ',32,num2str(mostPowerfullFrequency)),'fontweight','bold', 'fontsize',11);
    
    set(gca,'fontsize',11)
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    set(gcf, 'PaperPositionMode', 'auto');
    
    subplot(4,1,2)
    
    load(strcat(global_path2root,'scripts4paper', filesep ,'files', filesep ,'sampleFieldtripStruc.mat'))
    
    data = EGGTimeseries;
    nVoxels = size(data,1);
    
    % Define fieldtrip structure
    channelStr=cell(nVoxels,1);
    for iVoxel = 1:nVoxels
        channelList(iVoxel,1) = iVoxel;
        channelStr(iVoxel) = cellstr(mat2str(iVoxel));
    end
    
    dataStructure.hdr = EGG_downsampled.hdr;
    dataStructure.fsample = 0.5;
    dataStructure.time{1,1}  = [0:2:(size(data,2)*2)-1];
    dataStructure.label = channelStr;
    dataStructure.cfg = EGG_downsampled.cfg;
    dataStructure.trial{1,1} = data;
    
    len = dataStructure.fsample*120; % length of subtrials cfg.length s in samples
    dataStructure.sampleinfo=[1 max(dataStructure.time{1,1})*dataStructure.fsample];
    cfg = [];
    cfg.trl(:,1) = dataStructure.sampleinfo(1):(len/6):dataStructure.sampleinfo(2)-len+1;%trial start in samples from begining of raw data
    cfg.trl(:,2) = dataStructure.sampleinfo(1)+len-1:(len/6):dataStructure.sampleinfo(2);%trial ends in samples from begining of raw data
    cfg.trl(:,3) = 0; %offset of the trigger with respect to the trial
    data_trials = ft_redefinetrial(cfg,dataStructure);
    
    
    % Estimate spectrum of filtered timeseries to check if filter worked
    % properly
    cfg = [];
    cfg.method = 'mtmfft';
    cfg.taper = 'hanning';
    cfg.output = 'pow';
    cfg.pad = 1000;
    cfg.foilim = [1/120 0.1]; % 0 - 6 cpm
    cfg.keeptrials = 'no';
    
    frequencyWelchFiltered = ft_freqanalysis(cfg,data_trials);
    
    indexFrequenciesFiltered = find (frequencyWelchFiltered.freq >= 0.035 & frequencyWelchFiltered.freq <= 0.065);
    
    
    plot(frequencyWelchFiltered.freq(indexFrequenciesFiltered),frequencyWelchFiltered.powspctrm(indexFrequenciesFiltered),'-o','lineWidth',3);
    
    title(strcat('Welch EGG n°',num2str(bestChannel),32,'participant',sprintf('%.2d',subj_idx),32,'frequency ',32,num2str(frequencyWelchFiltered.freq(frequencyWelchFiltered.powspctrm == max(frequencyWelchFiltered.powspctrm)))),'fontweight','bold', 'fontsize',11);
    
    
    
    subplot(4,1,3)
    
    plot(30:2:868,angle(phaseXVolume),'-or','lineWidth',2)
    title(strcat('PhaseXVolume for subject' ,32,sprintf('%.2d',subj_idx),32,'channel ',num2str(bestChannel),32,'frequency ',32,num2str(mostPowerfullFrequency)), 'fontsize',11)
    
    xlabel('S','fontsize',11); ylabel('Angle in radians','fontsize',11);
    xlim([0,900])
    
    subplot(4,1,4)
    plot(30:2:868,EGGTimeseries,'-ob','lineWidth',2)
    title(strcat('AmplitudeXVolume for subject' ,32,sprintf('%.2d',subj_idx),32,'channel ',num2str(bestChannel),32,'frequency ',32,num2str(mostPowerfullFrequency)), 'fontsize',11)
    
    xlabel('S','fontsize',11); ylabel('Angle in radians','fontsize',11);
    xlim([0,900])
    
    if cfgMain.savePlots == 1
        
        print ('-dpng', '-painters', eval('filteredPlotFilename'))
        
        print ('-depsc2', '-painters', eval('filteredPlotFilename'))
        saveas(SanityPlot,strcat(filteredPlotFilename,'.fig'))
    end
end
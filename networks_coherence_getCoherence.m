function networks_coherence_getCoherence(subj_idx,cfgMain)


%{
This function computes the coherence between the EGG and all gastric
network regions and between all gastric network regions for a given
participant and store it into disk

Input:
Cluster timeseries
    Y:\Subjects\Subject13\Timeseries\PhasesAnalysis\ClusterTimeseries_csfr_S_13_kw3_CA0050
EGG amplitude
    Y:\Subjects\Subject13\Timeseries\EGGtimeseries\EGGTimeseriesFullband_S_13

Output:
Coherence clusters
    Y:\Subjects\Subject13\Timeseries\PhasesAnalysis\ClusterTimeseriesCoherence_csfr_S_13_kw3_CA0050
Power spectrum clusters
    Y:\Subjects\Subject13\Timeseries\PhasesAnalysis\ClusterTimeseriesSpectrum_csfr_S_13_kw3_CA0050

IR 28/06/2017

%}
%% get input - output  filename and files
clusterTimeseries_filename = global_filename(subj_idx,cfgMain,'clusterTimeseries_filename');
load(clusterTimeseries_filename)
clusterTimeseries_spectrum_filename = global_filename(subj_idx,cfgMain,'clusterTimeseries_spectrum_filename');
clusterTimeseries_coherence_filename = global_filename(subj_idx,cfgMain,'clusterTimeseries_coherence_filename');

insideBrain = tools_getIndexBrain('inside');

filename_EGG_amplitude = global_filename(subj_idx,cfgMain,'FilenameamplitudeXVolumeBestChannel_FULLBAND');
load(filename_EGG_amplitude)
filename_GlobalSignal = global_filename(subj_idx,cfgMain,'GlobalSignal_CSFr_FB_filename');
load(filename_GlobalSignal)


% EGG_FullBand = EGG_FullBand(cfgMain.beginCut:cfgMain.endCut); % this is
% now already curted in the coherence_getFULLbandEGG script

clusterTimeseries = clusterTimeseries(cfgMain.beginCut:cfgMain.endCut,:);

timeseries2Coherence = zeros(14,420) % 12 clusters plus EGG plus global
timeseries2Coherence(1,:) = EGG_FullBand;
timeseries2Coherence(2,:) = gs_timecourse';
timeseries2Coherence(3:14,:)= clusterTimeseries';
nVoxels = size(timeseries2Coherence,1);

%load base fieldtrip data structure
load(strcat(global_path2root,'scripts4paper\files\sampleFieldtripStruc.mat'))

labelsChannelsMAIN = {'EGG','Global','Ling','Cun','SMA_RCZp','SSr','CCZ','dPrec','PCC','TemParL','OccTemR','TemParR','M1SMA','OccSupR'};
labelsChannels = labelsChannelsMAIN;
clusterRegionsComparisons = timeseries2Coherence;
dataStructure.hdr = EGG_downsampled.hdr;
dataStructure.fsample = 0.5;
dataStructure.time{1,1}  = [0:2:(size(clusterRegionsComparisons,2)*2)-1];dataStructure.label = labelsChannels;%channelStr;
dataStructure.cfg = EGG_downsampled.cfg;
dataStructure.trial{1,1} = clusterRegionsComparisons;




%% computing fourier
cfgWelch = [];
cfgWelch.lengthWindow = 120; %seconds
cfgWelch.overlap = 6;% propotion i.e 1/6
len = dataStructure.fsample*cfgWelch.lengthWindow; % length of subtrials cfg.length s in samples
dataStructure.sampleinfo=[1 max(dataStructure.time{1,1})*dataStructure.fsample];

% Cut data into trials
cfg = [];
cfg.trl(:,1) = dataStructure.sampleinfo(1):(len/cfgWelch.overlap):dataStructure.sampleinfo(2)-len+1;%trial start in samples from begining of raw data
cfg.trl(:,2) = dataStructure.sampleinfo(1)+len-1:(len/cfgWelch.overlap):dataStructure.sampleinfo(2);%trial ends in samples from begining of raw data
cfg.trl(:,3) = 0; %offset of the trigger with respect to the trial
data_trials = ft_redefinetrial(cfg,dataStructure);

% Perform frequency analysis
cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.output = 'powandcsd';
cfg.pad = 1000;
cfg.foilim = [1/cfgWelch.lengthWindow 0.1]; % 0 - 6 cpm
cfg.keeptrials = 'no';

frequencyWelch = ft_freqanalysis(cfg,data_trials);

% Get coherence
cfg            = [];
cfg.method     = 'coh';
coherence      = ft_connectivityanalysis(cfg, frequencyWelch);

%% save
save(clusterTimeseries_spectrum_filename,'frequencyWelch')
save(clusterTimeseries_coherence_filename,'coherence')

%% Plots

plotDir = strcat (global_path2subject(subj_idx),'PreprocessingLog',filesep);

plotFilenameFreq = strcat(plotDir,'S_',sprintf('%.2d',subj_idx),'_Spectrum');
plotFilenameCohEGG = strcat(plotDir,'S_',sprintf('%.2d',subj_idx),'_CoherenceEGG');
plotFilenameCohGS = strcat(plotDir,'S_',sprintf('%.2d',subj_idx),'_CoherenceGS');

if cfgMain.savePlots == 1
    
    if cfgMain.plotFigures == 0;
        SanityPlotFreq = figure('visible','off');
        SanityPlotCohEGG = figure('visible','off');
        SanityPlotCohGS = figure('visible','off');


    else
        SanityPlotFreq = figure('visible','on');
        SanityPlotCohEGG = figure('visible','on');
        SanityPlotCohGS = figure('visible','on');

    end
  
    SanityPlotFreq;
    subplot(5,1,1)
    plot(frequencyWelch.freq,frequencyWelch.powspctrm(1,:)','LineWidth',3)
    legend(labelsChannelsMAIN{1})
    h = gca; set(h,'FontSize',16)
    title(['S',sprintf('%.2d',subj_idx),32,'Spectrum'],'fontsize',18)
    subplot(5,1,2)
    plot(frequencyWelch.freq,frequencyWelch.powspctrm(2:4,:)','LineWidth',3)
    legend(labelsChannelsMAIN{2:4})
    h = gca; set(h,'FontSize',16)

    subplot(5,1,3)
    plot(frequencyWelch.freq,frequencyWelch.powspctrm(5:7,:)','LineWidth',3)
    legend(labelsChannelsMAIN{5:7})
        h = gca; set(h,'FontSize',16)

    subplot(5,1,4)
    plot(frequencyWelch.freq,frequencyWelch.powspctrm(8:11,:)','LineWidth',3)
    legend(labelsChannelsMAIN{8:11})
        h = gca; set(h,'FontSize',16)

    subplot(5,1,5)
    plot(frequencyWelch.freq,frequencyWelch.powspctrm(12:14,:)','LineWidth',3)
    legend(labelsChannelsMAIN{12:14})
        h = gca; set(h,'FontSize',16)

    
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    set(gcf, 'PaperPositionMode', 'auto');
    
    print ('-dpng', '-painters', eval('plotFilenameFreq'))
    print ('-depsc2', '-painters', eval('plotFilenameFreq'))
    saveas(SanityPlotFreq,strcat(plotFilenameFreq,'.fig'))
    
    
    % coherence plot EGG
    SanityPlotCohEGG;
    subplot(4,1,1)
    plot(coherence.freq,coherence.cohspctrm(1:4,:)','LineWidth',3)
    legend(coherence.labelcmb{1:4,1})
    h = gca; set(h,'FontSize',16)
    title(['S',sprintf('%.2d',subj_idx),32,'Coherence with EGG'],'fontsize',18)
    
    subplot(4,1,2)
    plot(coherence.freq,coherence.cohspctrm(4:7,:)','LineWidth',3)
    legend(coherence.labelcmb{4:7,1})
    h = gca; set(h,'FontSize',16)

    subplot(4,1,3)
    plot(coherence.freq,coherence.cohspctrm(8:10,:)','LineWidth',3)
    legend(coherence.labelcmb{8:10,1})
    h = gca; set(h,'FontSize',16)

    subplot(4,1,4)
    plot(coherence.freq,coherence.cohspctrm(11:13,:)','LineWidth',3)
    legend(coherence.labelcmb{11:13,1})
    h = gca; set(h,'FontSize',16)

   
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    set(gcf, 'PaperPositionMode', 'auto');
    
    print ('-dpng', '-painters', eval('plotFilenameCohEGG'))
    print ('-depsc2', '-painters', eval('plotFilenameCohEGG'))
    saveas(SanityPlotCohEGG,strcat(plotFilenameCohEGG,'.fig'))
    
     
    
    % coherence plot Global Signal
    
        SanityPlotCohGS;
    subplot(4,1,1)
    plot(coherence.freq,coherence.cohspctrm(14:16,:)','LineWidth',3)
    legend(coherence.labelcmb{14:16,1})
    h = gca; set(h,'FontSize',16)
    title(['S',sprintf('%.2d',subj_idx),32,'Coherence with GS'],'fontsize',18)
    
    subplot(4,1,2)
    plot(coherence.freq,coherence.cohspctrm(17:19,:)','LineWidth',3)
    legend(coherence.labelcmb{17:19,1})
    h = gca; set(h,'FontSize',16)

    subplot(4,1,3)
    plot(coherence.freq,coherence.cohspctrm(20:22,:)','LineWidth',3)
    legend(coherence.labelcmb{20:22,1})
    h = gca; set(h,'FontSize',16)

    subplot(4,1,4)
    plot(coherence.freq,coherence.cohspctrm(23:25,:)','LineWidth',3)
    legend(coherence.labelcmb{23:25,1})
    h = gca; set(h,'FontSize',16)

   
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    set(gcf, 'PaperPositionMode', 'auto');
    
    print ('-dpng', '-painters', eval('plotFilenameCohGS'))
    print ('-depsc2', '-painters', eval('plotFilenameCohGS'))
    saveas(SanityPlotCohGS,strcat(plotFilenameCohGS,'.fig'))
    
    
end

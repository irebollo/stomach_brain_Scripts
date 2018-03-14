function heart_cohHRV_BOLD(subj_idx,cfgMain)

insideBrain = tools_getIndexBrain('inside');

%% get coherence
% 
% figure
% plot(t_int,ibi_int)
%%  load and preproces HRV ts
% downsample to 0.05 Hz

HRV_timeseries_filename = global_filename(subj_idx,cfgMain,'HRV_timeseries');
load(HRV_timeseries_filename)


timeseries2Coherence = ibi_int';
load(strcat(global_path2root,'\scripts4paper\files\sampleFieldtripStruc.mat'))
labelsChannelsMAIN = {'HRV'};
labelsChannels = labelsChannelsMAIN;
clusterRegionsComparisons = timeseries2Coherence;
dataStructure.hdr = EGG_downsampled.hdr;
dataStructure.fsample = 1;
dataStructure.time{1,1}  = [0:1/dataStructure.fsample:(size(clusterRegionsComparisons,2))-1];
dataStructure.label = labelsChannels;%channelStr;
dataStructure.cfg = EGG_downsampled.cfg;
dataStructure.trial{1,1} = clusterRegionsComparisons;
dataStructure.sampleinfo = [1 length(timeseries2Coherence)]

% a data.sampleinfo field in your data structure
% > that is consistent with the actual data. In general, data.sampleinfo
% > is an Nx2 matrix containing, for each trial in the data, the indices
% > of the begin and end sample of that trial, with respect to the
% > original data set (on disk)

disp('Resampling...')
cfg = [];  %initialize configuration structure
cfg.detrend = 'no'; % remove linear trend from the data (done per trial)
cfg.demean = 'yes';
cfg.resamplefs= 0.5; % 4 x top-freq (15 cpm = 0.25 Hz) - Nyquist = 30 cpm  frequency at which the data will be resampled
HRV_downsampled = ft_resampledata(cfg,dataStructure); % This procedure also lowpass filter the data at half the new sr 


HRV_downsampled.trial{1,1} = ft_preproc_polyremoval (HRV_downsampled. trial{1,1}, 2);
HRV_downsampled.trial{1,1} = ft_preproc_standardize (HRV_downsampled. trial{1,1});

% HRV_downsampled.trial{1,1} (cfgMain.beginCut:cfgMain.endCut)
% figure
% plot(HRV_downsampled.trial{1,1})


%% load preprocess Brain

filename_bold_input = global_filename(subj_idx,cfgMain,'BOLDTimeseriesFilename');
load(filename_bold_input)


% 
% timeseries2Coherence(1,:) = mean_ssi_timeseries;
% timeseries2Coherence(2,:) = mean_precuneus_timeseries;
% nVoxels = size(timeseries2Coherence,1);


BOLDtimeseries.trialVector = BOLDtimeseries.trialVector(cfgMain.beginCut:cfgMain.endCut,insideBrain);
BOLDtimeseries.BOLDtsNOPOLI = ft_preproc_polyremoval (BOLDtimeseries.trialVector', 2);
BOLD_zscored = ft_preproc_standardize (BOLDtimeseries.BOLDtsNOPOLI);



data = [HRV_downsampled.trial{1,1}(cfgMain.beginCut:cfgMain.endCut);BOLD_zscored];


nVoxels = size(data,1);

%load base fieldtrip data structure
load(strcat(global_path2root,'\scripts4paper\files\sampleFieldtripStruc.mat'))

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

cfgWelch = [];
cfgWelch.keeptrials = 'no';
cfgWelch.lengthWindow = 120;
cfgWelch.overlap = 6;

len = dataStructure.fsample*cfgWelch.lengthWindow; % length of subtrials cfg.length s in samples
dataStructure.sampleinfo=[1 max(dataStructure.time{1,1})*dataStructure.fsample];
cfg = [];
cfg.trl(:,1) = dataStructure.sampleinfo(1):(len/cfgWelch.overlap):dataStructure.sampleinfo(2)-len+1;%trial start in samples from begining of raw data
cfg.trl(:,2) = dataStructure.sampleinfo(1)+len-1:(len/cfgWelch.overlap):dataStructure.sampleinfo(2);%trial ends in samples from begining of raw data
cfg.trl(:,3) = 0; %offset of the trigger with respect to the trial
data_trials = ft_redefinetrial(cfg,dataStructure);


% Estimate spectrum
cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.output = 'powandcsd';
cfg.pad = 1000;
cfg.foilim = [1/cfgWelch.lengthWindow 0.25]; % 0 - 6 cpm
cfg.keeptrials = cfgWelch.keeptrials;

cfg.channelcmb=cell(71759,2);
for iVoxel = 1:71759
   cfg.channelcmb{iVoxel,1} = '1';
   cfg.channelcmb{iVoxel,2} = char(channelStr(iVoxel+1));
end

frequencyWelch = ft_freqanalysis(cfg,data_trials);



cfg            = [];
cfg.channelcmb=cell(71759,2);
for iVoxel = 1:71759
   cfg.channelcmb{iVoxel,1} = '1';
   cfg.channelcmb{iVoxel,2} = char(channelStr(iVoxel+1));
end
cfg.method     = 'coh';
coherence             = ft_connectivityanalysis(cfg, frequencyWelch);


%% plot specturm

figure
plot(frequencyWelch.freq,frequencyWelch.powspctrm(1,:))

figure
plot(frequencyWelch.freq,mean(coherence.cohspctrm).*mean(coherence.cohspctrm))


emptyBrain = zeros (53,63,46); % empty 3d Volume for storing empirical PLV
emptyBrain = emptyBrain(:);
emptyBrain(insideBrain) = mean(coherence.cohspctrm.*coherence.cohspctrm,2) ;
emptyBrain = reshape(emptyBrain,53,63,46);



%% Save

       HRV_BOLD_spectrum_filename =  global_filename(subj_idx,cfgMain,'HRV_BOLD_spectrum');
       HRV_BOLD_coherence_filename = global_filename(subj_idx,cfgMain,'HRV_BOLD_coherence');
       HRV_BOLD_coherence_map_filename = global_filename(subj_idx,cfgMain,'HRV_BOLD_coherence_map');

       save(HRV_BOLD_spectrum_filename,'frequencyWelch')
      save(HRV_BOLD_coherence_filename,'coherence')

tools_writeMri(emptyBrain,HRV_BOLD_coherence_map_filename)

end

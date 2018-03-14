function heart_buildHRVregressors(subj_idx)


% 0.15 to 0.4
% 0.04 to 0.15



% filter
% hilbert
% store envelope
% make fisrt level GLM

cfgMain = global_getcfgmain
% subj_idx = 13

load(['/media/ignacio/IREBOLLOENS/HRVts/IBIts_S' sprintf('%.2d',subj_idx) '.mat'])
%%


% 
% HRV_timeseries_filename = global_filename(subj_idx,cfgMain,'HRV_timeseries')
% load(HRV_timeseries_filename)
% 
% 
% HRV_notdownsampled = ibi_int';
% load(strcat(global_path2root,'\scripts4paper\files\sampleFieldtripStruc.mat'))
load('/home/ignacio/Matlab/scripts4paper/files/sampleFieldtripStruc.mat')


% labelsChannelsMAIN = {'HRV'};
% labelsChannels = labelsChannelsMAIN;
% clusterRegionsComparisons = HRV_notdownsampled;
% dataStructure.hdr = EGG_downsampled.hdr;
% dataStructure.fsample = 1;
% dataStructure.time{1,1}  = [0:1/dataStructure.fsample:(size(clusterRegionsComparisons,2))-1];
% dataStructure.label = labelsChannels;%channelStr;
% dataStructure.cfg = EGG_downsampled.cfg;
% dataStructure.trial{1,1} = clusterRegionsComparisons;
% dataStructure.sampleinfo = [1 length(HRV_notdownsampled)]

% disp('Resampling...')
% cfg = [];  %initialize configuration structure
% cfg.detrend = 'no'; % remove linear trend from the data (done per trial)
% cfg.demean = 'yes';
% cfg.resamplefs= 0.5; % 4 x top-freq (15 cpm = 0.25 Hz) - Nyquist = 30 cpm  frequency at which the data will be resampled
% HRV_downsampled = ft_resampledata(cfg,dataStructure); % This procedure also lowpass filter the data at half the new sr 

% 
% HRV_polyremoval= ft_preproc_polyremoval (ibi_int, 2);
% HRV_downsampled.trial{1,1} = ft_preproc_standardize (HRV_downsampled. trial{1,1});
% 

%% LF HRV

filter_frequency_spread=0.06 ; % In hz Half widfth or fullwidth????
centerFrequency = 0.1
sr = 1; % 1 TR = 2s
filterOrder=(cfgMain.fOrder*fix(sr/(centerFrequency-filter_frequency_spread))-1);%in nsamples
transition_width= cfgMain.transitionWidth/100; % in normalised units
LFfilteredHeart=tools_bpFilter(ibi_int,sr,filterOrder,centerFrequency,filter_frequency_spread,transition_width,cfgMain.filterType);
HilbertLFHRV = hilbert(LFfilteredHeart);
phaseLFHRV = angle(HilbertLFHRV); % Cut data to have the same length as EGG (cut this way to get rid of fmri edge artifact on EGG)
AmplitudeEnvelopeLFHRV = abs(HilbertLFHRV); % Cut data to have the same length as EGG (cut this way to get rid of fmri edge artifact on EGG)

figure
plot(LFfilteredHeart)
hold on
plot(AmplitudeEnvelopeLFHRV,'r')

%% HF HRV

filter_frequency_spread=0.125 ; % In hz Half widfth or fullwidth????
centerFrequency = 0.275
sr = 1; % 1 TR = 2s
filterOrder=(cfgMain.fOrder*fix(sr/(centerFrequency-filter_frequency_spread))-1);%in nsamples
transition_width= cfgMain.transitionWidth/100; % in normalised units
HFfilteredHeart=tools_bpFilter(ibi_int,sr,filterOrder,centerFrequency,filter_frequency_spread,transition_width,cfgMain.filterType);
HilbertHFHRV = hilbert(HFfilteredHeart);
phaseHFHRV = angle(HilbertHFHRV); % Cut data to have the same length as EGG (cut this way to get rid of fmri edge artifact on EGG)
AmplitudeEnvelopeHFHRV = abs(HilbertHFHRV); % Cut data to have the same length as EGG (cut this way to get rid of fmri edge artifact on EGG)

figure
plot(ibi_int)
hold on
plot(AmplitudeEnvelopeLFHRV,'r')
plot(AmplitudeEnvelopeHFHRV,'G')


% downsample lFHRV envelope

% load(strcat(global_path2root,'\scripts4paper\files\sampleFieldtripStruc.mat'))
load('/home/ignacio/Matlab/scripts4paper/files/sampleFieldtripStruc.mat')


labelsChannelsMAIN = {'LFHRV','HFHRV'};
labelsChannels = labelsChannelsMAIN;
clusterRegionsComparisons = [AmplitudeEnvelopeLFHRV';AmplitudeEnvelopeHFHRV'];
dataStructure.hdr = EGG_downsampled.hdr;
dataStructure.fsample = 1;
dataStructure.time{1,1}  = [0:1/dataStructure.fsample:(size(clusterRegionsComparisons,2))-1];
dataStructure.label = labelsChannels;%channelStr;
dataStructure.cfg = EGG_downsampled.cfg;
dataStructure.trial{1,1} = clusterRegionsComparisons;
dataStructure.sampleinfo = [1 length(AmplitudeEnvelopeLFHRV)]

disp('Resampling...')
cfg = [];  %initialize configuration structure
cfg.detrend = 'no'; % remove linear trend from the data (done per trial)
cfg.demean = 'yes';
cfg.resamplefs= 0.5; % 4 x top-freq (15 cpm = 0.25 Hz) - Nyquist = 30 cpm  frequency at which the data will be resampled
HRV_downsampled = ft_resampledata(cfg,dataStructure); % This procedure also lowpass filter the data at half the new sr 

figure
% plot(ibi_int)
hold on
plot(HRV_downsampled.trial{1,1}(1,:),'-or')
plot(HRV_downsampled.trial{1,1}(2,:),'-og')

HFHRVregressor = zscore(HRV_downsampled.trial{1,1}(2,:)')
LFHRVregressor = zscore(HRV_downsampled.trial{1,1}(1,:)')

save(['/media/ignacio/IREBOLLOENS/HRVts/HRVregressors_S' sprintf('%.2d',subj_idx) '.mat'])




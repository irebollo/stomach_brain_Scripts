function frequencyWelch = tools_estimateMRIspectrum (data,cfgWelch)


%data = VoxelxTime data
% cfg.keeptrials =1 / =0 must include keep trials
% cfg.overlap: proportion of overlap: 1/overlap
% cfg.lengthWindow of windows in seconds


%%
nVoxels = size(data,1);

%load base fieldtrip data structure
load(strcat(global_path2root,'scripts4paper\files\sampleFieldtripStruc.mat'))

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
cfg.output = 'pow';
cfg.pad = 1000;
cfg.foilim = [1/cfgWelch.lengthWindow 0.1]; % 0 - 6 cpm
cfg.keeptrials = cfgWelch.keeptrials;

frequencyWelch = ft_freqanalysis(cfg,data_trials);




end


%indexFrequencies = find (frequencyWelch.freq >= 0.0200 & frequencyWelch.freq <= 0.0800);
% 
% figure;plot(frequencyWelch.freq,frequencyWelch.powspctrm(1,:),'lineWidth',3);title(strcat('S',num2str(subj_idx),32,'V1',32,dataType),'FontSize',20);ylim([0 0.5])
% figure;plot(frequencyWelch.freq,frequencyWelch.powspctrm(2,:),'lineWidth',3);title(strcat('S',num2str(subj_idx),32,'V2',32,dataType),'FontSize',20);ylim([0 0.5])
% figure;plot(frequencyWelch.freq,frequencyWelch.powspctrm(3,:),'lineWidth',3);title(strcat('S',num2str(subj_idx),32,'V3',32,dataType),'FontSize',20);ylim([0 0.5])
% 

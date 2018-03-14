% function heart_cohHRV_BOLD(subj_idx,cfgMain)



subjects = global_subjectList
cfgMain = global_getcfgmain

insideBrain = tools_getIndexBrain('inside');

peaksAllsubjects = global_getEGGpeaks;



powerAllSubjects = zeros(30,length(insideBrain),1);


for iSubject=1:length(subjects)
    subj_idx = subjects(iSubject)
%% load preprocess Brain

filename_bold_input = global_filename(subj_idx,cfgMain,'BOLD_filtered_fullbandFilename');
load(filename_bold_input)


mostPowerfullFrequency = peaksAllsubjects(iSubject,3);

% 
% timeseries2Coherence(1,:) = mean_ssi_timeseries;
% timeseries2Coherence(2,:) = mean_precuneus_timeseries;
% nVoxels = size(timeseries2Coherence,1);



data = [BOLD_filtered_zscored(insideBrain,:)];


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
cfg.output = 'pow';
cfg.pad = 1000;
cfg.foilim = [mostPowerfullFrequency-0.001, mostPowerfullFrequency+0.001]; % 0 - 6 cpm
% cfg.foilim = [1/cfgWelch.lengthWindow 0.25]; % 0 - 6 cpm

frequencyWelch = ft_freqanalysis(cfg,data_trials);

powerAllSubjects(iSubject,:,1) = nanmean(frequencyWelch.powspctrm,2);

filename_BOLDpoweratGastric = global_filename(subj_idx,cfgMain,'filename_BOLDpoweratGastric');
save(filename_BOLDpoweratGastric,'frequencyWelch')

end


emptyBrain = zeros (53,63,46); % empty 3d Volume for storing empirical PLV
emptyBrain = emptyBrain(:);
emptyBrain(insideBrain) = nanmean(powerAllSubjects);
emptyBrain = reshape(emptyBrain,53,63,46);
tools_writeMri(emptyBrain,'Y:\ClusterResults\power\PowerAtGastric')

nhist(nanmean(powerAllSubjects))


%% load CS of all subjects and perform correlation

gastricNetwork = global_getGastricNetwork;
gastricNetwork = gastricNetwork(:);

insideBrain = tools_getIndexBrain('inside');

emp_brain = zeros(length(subjects),length(insideBrain));
emp_gastric = zeros(length(subjects),length(find(gastricNetwork)));
chance_brain = zeros(length(subjects),length(insideBrain));
chance_gastric = zeros(length(subjects),length(find(gastricNetwork)));

for iSubject=1:length(subjects)
    subj_idx = subjects(iSubject)
    
    filename_PLV= global_filename(subj_idx,cfgMain,'PLVXVoxelFilename_csfr');
    empPLV = ft_read_mri([filename_PLV,'.nii']);
    filename_surrogate= global_filename(subj_idx,cfgMain,'medianRotationFilename_csfr');
    chancePLV = ft_read_mri([filename_surrogate,'.nii']);
    
    emp_brain(iSubject,:)= empPLV.anatomy(insideBrain);
    emp_gastric(iSubject,:)= empPLV.anatomy(gastricNetwork);
    chance_brain(iSubject,:)= chancePLV.anatomy(insideBrain);
    chance_gastric(iSubject,:)= chancePLV.anatomy(gastricNetwork);
end

cs_gastric = emp_gastric-chance_gastric;
cs_brain = emp_brain-chance_brain;

mean_cs_gastric = mean(cs_gastric);
mean_cs_brain= mean(cs_brain);

figure
plot(nanmean(powerAllSubjects),mean_cs_brain,'ob')

figure
plot(nanmean(powerAllSubjects(:,gastricNetwork(insideBrain))),mean_cs_gastric,'ob')


% iterate through subjects, 
for iSubj = 1:length(subjects)
subj_idx = subjects(iSubj)

% 
% % load betas of framewise displacement
% filename_BetasInstantaneousMovementEstimates = global_filename(subj_idx,cfgMain,'filename_BetasInstantaneousMovementEstimates');
% load(filename_BetasInstantaneousMovementEstimates)
% BetasMovement(iSubj,:) = betas_instantaneousMovement;
% 
% % load coupling strenght
% couplingStrenghtFilename_csfr = global_filename(subj_idx,cfgMain,'couplingStrenghtFilename_csfr');
% load(couplingStrenghtFilename_csfr);
% cs_allsubjects(iSubj,:) = couplingstrenght;
% 
% 
% % identify and exclude NaNs
ind_nan_pow = isnan(powerAllSubjects(iSubj,:));
ind_nan_cs = isnan(cs_brain(iSubj,:));
not_nan = ~ind_nan_pow & ~ind_nan_cs ;

[r p ] = corrcoef(powerAllSubjects(iSubj,not_nan),cs_brain(iSubj,not_nan));
Pow_cs_r_wholeBrain(iSubj) = r(3)

% we are only interested in the correlation inside the gastric network
% BetasMovement_gasnet = BetasMovement(gasnet);
% couplingstrenght_gasnet = couplingstrenght(gasnet);

ind_nan_pow_gasnet = isnan(powerAllSubjects(iSubj,gastricNetwork(insideBrain)));
ind_nan_cs_gasnet = isnan(cs_gastric(iSubj,:));
not_nan = ~ind_nan_pow_gasnet & ~ind_nan_cs_gasnet ;

[r p ] = corrcoef(powerAllSubjects(iSubj,gastricNetwork(insideBrain)),cs_gastric(iSubj,:));
Pow_cs_r_gasnet(iSubj) = r(3)

end


% individual subjects
    figure
    Pow_cs_r_gasnet_fisherZ = fisherz(Pow_cs_r_gasnet)
    violin(Pow_cs_r_gasnet')
    h=gca
    set(h,'fontsize',16)
    [h p ci stats] = ttest(Pow_cs_r_gasnet_fisherZ)
    title(['Individual r values Gastric Network p ' num2str(p)])
    
    
    % individual subjects
    figure
    Pow_cs_r_brain_fisherZ = fisherz(Pow_cs_r_wholeBrain)
    violin(Pow_cs_r_wholeBrain')
    h=gca
    set(h,'fontsize',16)
    [h p ci stats] = ttest(Pow_cs_r_brain_fisherZ)
    title(['Individual r values Gastric Network p ' num2str(p)])
   

% Prior H1: corresponds to an effect differing from 0 with a p-value of 0.05
nobs = 30; % number of observations
xref = +1.699; % significant t value for one-sided ttest 29 degrees of freedom reference effect size (significant effect)
% source: http://www.ttable.org/uploads/2/1/7/9/21795380/9754276.png?852

% disp('FWD:')
% load([global_path2root 'data4paper' filesep 'rFWDxCS_subjects.mat'])
% data obtained form script paper_control_corrMovement_CS

xdat = Pow_cs_r_brain_fisherZ

% perform ttest
[bf_log10]= my_ttest_bayes(xdat, xref);
% interpret bf
[res_Bayes, bf] = interpret_Bayes(bf_log10)
% intert it to 
bf_inverted_unlog = 1/(10^bf)
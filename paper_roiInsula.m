%{
paper_ROIinsula
Perform analysis presneted in section: Sub-threshold coupling in the right posterior insula
Load Deen 2011 parcelletion of right and left anterior, posterior and
middle insula. 
Compute for each subject average BOLD time-series in each of these ROIs, calculate empirical
and chance and perform a ttest comparing them

Inputs are stored in; Y:\MasksFromOtherPapers\Deen_Insula\
Output, results of the ttest in the command line, bonferroni corrected for multiple
comparisons

IR 04/07/2017
%}
%% load cfg 

subjects = global_subjectList;
cfgMain = global_getcfgmain;

%% Load masks
insideBrain = tools_getIndexBrain('inside');

% Load masks
rPI = ft_read_mri([global_path2root 'MasksFromOtherPapers' filesep 'Deen_Insula' filesep 'RPI_3mm.nii']);
rPI = logical(rPI.anatomy(:));
rPI = rPI(insideBrain);

rDAI = ft_read_mri([global_path2root 'MasksFromOtherPapers' filesep 'Deen_Insula' filesep 'RDAI_3mm.nii']);
rDAI = logical(rDAI.anatomy(:));
rDAI = rDAI(insideBrain);

rVAI = ft_read_mri([global_path2root 'MasksFromOtherPapers' filesep 'Deen_Insula' filesep 'RVAI_3mm.nii']);
rVAI = logical(rVAI.anatomy(:));
rVAI = rVAI(insideBrain);

lPI = ft_read_mri([global_path2root 'MasksFromOtherPapers' filesep 'Deen_Insula' filesep 'LPI_3mm.nii']);
lPI = logical(lPI.anatomy(:));
lPI = lPI(insideBrain);

lDAI = ft_read_mri([global_path2root filesep 'MasksFromOtherPapers' filesep 'Deen_Insula' filesep 'LDAI_3mm.nii']);
lDAI = logical(lDAI.anatomy(:));
lDAI = lDAI(insideBrain);

lVAI = ft_read_mri([global_path2root filesep 'MasksFromOtherPapers' filesep 'Deen_Insula' filesep 'LVAI_3mm.nii']);
lVAI = logical(lVAI.anatomy(:));
lVAI = lVAI(insideBrain);


gastricNetwork = global_getGastricNetwork;
gastricNetwork= gastricNetwork(:);
gastricNetwork=gastricNetwork(insideBrain);


gray_mask_filename = 'Y:\scripts4paper\files\gray_mask.hdr'

grayMatter = ft_read_mri(gray_mask_filename);
grayMatter = logical(grayMatter.anatomy);
grayMatter=grayMatter(:);
grayMatter=grayMatter(insideBrain);

% Load the information about the peaks of the EGG
peaksAllsubjects = global_getEGGpeaks;


%initialize

empPLV_ROI = zeros(18,length(subjects));
chancePLV_ROI= zeros(18,length(subjects));

chancePLV_medianGray = zeros(1,length(subjects));
empPLV_medianGray = zeros(1,length(subjects));
empPLV_medianGastric= zeros(1,length(subjects));
chancePLV_medianGastric= zeros(1,length(subjects));

clusterMap = ft_read_mri('Y:\ClusterResults\kw3\CA0050\Cluster_nR10000_CA0050_kw3_fir2_fspread_015_fOrder_5_tw_15csfr_ClusterMap.nii')
clusterMap = clusterMap.anatomy(:);
clusterMap = clusterMap(insideBrain);
% nhist(clusterMap(find(clusterMap)));


%% Loop subject
for iSubj = 1:length(subjects)
    subj_idx = subjects(iSubj)
%% Load and filter BOLD timeseries from subject 

EGGPhaseXVolumeFilename = global_filename(subj_idx,cfgMain,strcat('EGGPhaseXVolumeFilename'));
load(EGGPhaseXVolumeFilename)

BOLDTimeseriesFilename = global_filename(subj_idx,cfgMain,strcat('filename_',cfgMain.Timeseries2Regress,'_Residuals_FB'));
timeseries = load(BOLDTimeseriesFilename);

indPeak = find (peaksAllsubjects(:,1) == subj_idx);
mostPowerfullFrequency = peaksAllsubjects(indPeak,3);


centerFrequency = mostPowerfullFrequency; %
filter_frequency_spread=cfgMain.frequencySpread/1000; % In hz
sr = 0.5; % 1 TR = 2s
filterOrder=(cfgMain.fOrder*fix(sr/(centerFrequency-filter_frequency_spread))-1);%in nsamples
transition_width= cfgMain.transitionWidth/100; % in normalised units
filteredMRI=tools_bpFilter(eval(strcat('timeseries.',cell2mat(fields(timeseries)))),sr,filterOrder,centerFrequency,filter_frequency_spread,transition_width,cfgMain.filterType);

filteredMRI = filteredMRI(cfgMain.beginCut:cfgMain.endCut,:);

%% Timeseries of ROI

timeseriesROI = zeros(420,18);
timeseriesROI(:,1) = nanmean(filteredMRI(:,rPI),2);
timeseriesROI(:,2) = nanmean(filteredMRI(:,rDAI),2);
timeseriesROI(:,3) = nanmean(filteredMRI(:,rVAI),2);
timeseriesROI(:,4) = nanmean(filteredMRI(:,lPI),2);
timeseriesROI(:,5) = nanmean(filteredMRI(:,lDAI),2);
timeseriesROI(:,6) = nanmean(filteredMRI(:,lVAI),2);

timeseriesROI(:,7) = nanmean(filteredMRI(:,clusterMap==1),2);
timeseriesROI(:,8) = nanmean(filteredMRI(:,clusterMap==2),2);
timeseriesROI(:,9) = nanmean(filteredMRI(:,clusterMap==3),2);
timeseriesROI(:,10) = nanmean(filteredMRI(:,clusterMap==4),2);
timeseriesROI(:,11) = nanmean(filteredMRI(:,clusterMap==5),2);
timeseriesROI(:,12) = nanmean(filteredMRI(:,clusterMap==6),2);
timeseriesROI(:,13) = nanmean(filteredMRI(:,clusterMap==7),2);
timeseriesROI(:,14) = nanmean(filteredMRI(:,clusterMap==8),2);
timeseriesROI(:,15) = nanmean(filteredMRI(:,clusterMap==9),2);
timeseriesROI(:,16) = nanmean(filteredMRI(:,clusterMap==10),2);
timeseriesROI(:,17) = nanmean(filteredMRI(:,clusterMap==11),2);
timeseriesROI(:,18) = nanmean(filteredMRI(:,clusterMap==12),2);

phaseROI = hilbert (timeseriesROI);

%% EmpPLV

    phaseDifferenceEmp = bsxfun (@minus , angle(phaseROI), angle(phaseXVolume )');

    PLV_emp =  abs (mean (exp (1i* phaseDifferenceEmp ) ) ); %
    
%% rotations PLV
indexRotations=31:390; % Rotating at least one minute (30 TR = 60s) at the beggining or end
rotatedPhaseEGG = zeros(length(indexRotations),length(phaseXVolume));
for iRotation = 1 : length(indexRotations)
rotatedPhaseEGG(iRotation,:) = circshift(phaseXVolume,[0 indexRotations(iRotation)]);
end

RPLV = zeros(length(indexRotations),18); % R from rotated

for iRotation = 1 : length(indexRotations)
    currentPhaseEGG = rotatedPhaseEGG (iRotation,:) ;
    phaseDifference = bsxfun (@minus , angle(phaseROI), angle(currentPhaseEGG )');
    PLV =  abs (mean (exp (1i* phaseDifference ) ) ); %
    RPLV(iRotation,:) = PLV;% timeseries_get_PLV(phaseMRI,rotatedPhaseEGG(iRotation,:)');
end
medianPLV = median(RPLV,1);

empPLV_ROI (:,iSubj)= PLV_emp
chancePLV_ROI(:,iSubj)= medianPLV;


%% Get mean EMP and Chance in regions not gastric not insula
% 
% % load median PLV all voxels
% medianRotationFilename = global_filename(subj_idx,cfgMain,strcat('medianRotationFilename_',cfgMain.Timeseries2Regress)); % output
% medianRotationWholeBrain = ft_read_mri([medianRotationFilename ,'.nii'])
% medianRotationWholeBrain = medianRotationWholeBrain.anatomy(:);
% medianRotationWholeBrain = medianRotationWholeBrain(insideBrain);
% 
% chancePLV_medianGastric=  nanmedian(medianRotationWholeBrain(gastricNetwork));
% 
% 
% medianRotationWholeBrain(gastricNetwork)=nan;
% medianRotationWholeBrain(lVAI)=nan;
% medianRotationWholeBrain(lDAI)=nan;
% medianRotationWholeBrain(rVAI)=nan;
% medianRotationWholeBrain(rDAI)=nan;
% medianRotationWholeBrain(rPI)=nan;
% 
% chancePLV_medianGray(iSubj)= nanmedian(medianRotationWholeBrain(grayMatter));
% 
% 
% 
% % medianRotationWholeBrain = medianRotationWhol
% % load empirical PLV all voxels
% PLVXVoxelFilename = global_filename(subj_idx,cfgMain,strcat('PLVXVoxelFilename_',cfgMain.Timeseries2Regress)); % output filename
% empPLVWholeBrain = ft_read_mri([PLVXVoxelFilename,'.nii'])
% empPLVWholeBrain=empPLVWholeBrain.anatomy(:);
% empPLVWholeBrain=empPLVWholeBrain(insideBrain);
% 
% 
% empPLV_medianGastric(iSubj)=nanmedian(empPLVWholeBrain(gastricNetwork))
% 
% 
% empPLVWholeBrain(gastricNetwork)=nan;
% empPLVWholeBrain(lVAI)=nan;
% empPLVWholeBrain(lDAI)=nan;
% empPLVWholeBrain(rVAI)=nan;
% empPLVWholeBrain(rDAI)=nan;
% empPLVWholeBrain(rPI)=nan;
% 
% empPLV_medianGray(iSubj)= nanmedian(empPLVWholeBrain(grayMatter));
% 
% 
% chancePLV_medianGray(iSubj)= nanmedian(medianRotationWholeBrain(grayMatter));


end




%% get effect sizes for each node

meanEMP_roi = mean(empPLV_ROI,2);
meanChance_roi = mean(chancePLV_ROI,2);

stdEMP_roi = std(empPLV_ROI,[],2);
stdChance_roi = std(chancePLV_ROI,[],2);

cohenD_roi = (meanEMP_roi-meanChance_roi)./((stdEMP_roi+stdChance_roi)/2)

mean(cohenD_roi(7:end))
std(cohenD_roi(7:end))
%% First ttest
[h p ci stats] = ttest2(empPLV_ROI',chancePLV_ROI');

%% bayes
nobs = 30; % number of observations
xref = +1.699; % significant t value for one-sided ttest 29 degrees of freedom reference effect size (significant effect)
% source: http://www.ttable.org/uploads/2/1/7/9/21795380/9754276.png?852

for iRoi =1:6
xdat = chancePLV_ROI(iRoi,:)-empPLV_ROI(iRoi,:)

[bf_log10(iRoi)]= my_ttest_bayes(xdat, xref);
[res_Bayes(iRoi,:), bf(iRoi)] = interpret_Bayes(bf_log10(iRoi))
.(iRoi) = 1/(10^bf(iRoi))

end

%% Bayes against median gray cs

median_cs_gray = empPLV_medianGray - chancePLV_medianGray;
median_cs_gastric = empPLV_medianGastric - chancePLV_medianGastric;
difference_gray_gastric = (median_cs_gastric - median_cs_gray)/2

coupling_strenght_roi = empPLV_ROI - chancePLV_ROI


for iRoi=1:6
   [h p ci stats]= ttest2(coupling_strenght_roi(iRoi,:),difference_gray_gastric)
   pp(iRoi)= p  
   tt(iRoi) = stats.tstat
end

for iRoi=1:6
t1smpbf(tt(iRoi),30)
end

for iRoi=1:6
   mean(coupling_strenght_roi(iRoi,:)-difference_gray_gastric)
end




nobs = 30; % number of observations
xref = +1.699; % significant t value for one-sided ttest 29 degrees of freedom reference effect size (significant effect)
% xref = +2.042; % significant t value for two-sided ttest 29 degrees of freedom reference effect size (significant effect)

% source: http://www.ttable.org/uploads/2/1/7/9/21795380/9754276.png?852

clear bf_log10 res_Bayes bf bayesFactor
for iRoi =1:6
xdat = coupling_strenght_roi(iRoi,:)-median_cs_gray

xref = mean(xdat)/std(xdat)
xrefIroi(iRoi) = xref
[bf_log10(iRoi)]= my_ttest_bayes(xdat, xref);
[res_Bayes{iRoi}, bf(iRoi)] = interpret_Bayes(bf_log10(iRoi))
bayesFactor(iRoi) = 1/(10^bf(iRoi))

end




nobs = 30; % number of observations
% xref = +1.699; % significant t value for one-sided ttest 29 degrees of freedom reference effect size (significant effect)
xref = +2.042; % significant t value for two-sided ttest 29 degrees of freedom reference effect size (significant effect)


figure
bar(mean(coupling_strenght_roi,2)-mean(median_cs_gray))
legend

xticks([1:6])
xticklabels({'rPI','rDAI','rVAI','lPI','lDAI','lVAI'})
ylabel('Coupling strenght')


clear bf_log10 res_Bayes bf bayesFactor

% xref = +2.042; % significant t value for two-sided ttest 29 degrees of freedom reference effect size (significant effect)
% xref = +1.699; % significant t value for one-sided ttest 29 degrees of freedom reference effect size (significant effect)

for iRoi =1:6
xdat = coupling_strenght_roi(iRoi,:)-median_cs_gray

xref = mean(xdat)/std(xdat)
xrefIroi(iRoi) = xref
[bf_log10(iRoi)]= my_ttest_bayes(xdat, xref);
[res_Bayes{iRoi}, bf(iRoi)] = interpret_Bayes(bf_log10(iRoi))
bayesFactor(iRoi) = 1/(10^bf(iRoi))

end

figure
bar(bf)

figure
bar(bf_log10)

xticks([1:6])
xticklabels({'rPI','rDAI','rVAI','lPI','lDAI','lVAI'})
ylabel('BF log 10')


%% dummy test

clear bf_log10 res_Bayes bf bayesFactor
% xdat = coupling_strenght_roi(iRoi,:)-median_cs_gray


xdat = randn(1,100000);
xref = mean(xdat)/std(xdat)
% xrefIroi(iRoi) = xref
[bf_log10]= my_ttest_bayes(xdat, xref);
[res_Bayes, bf] = interpret_Bayes(bf_log10)
bayesFactor = 1/(10^bf)


xdat = randn(1,100);
xdat= xdat +1;
xref = mean(xdat)/std(xdat)
% xrefIroi(iRoi) = xref
[bf_log10]= my_ttest_bayes(xdat, xref);
[res_Bayes, bf] = interpret_Bayes(bf_log10)
bayesFactor = 1/(10^bf)

%%
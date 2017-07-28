function networks_coherence_getClusterTimeseries(subj_idx,cfgMain)

%{
Compute and stores the average BOLD time series across each significant
cluster of the gastric network. Later used for coherence and angle analysis. 

Input: Bold timeseries 
Y:\Subjects\Subject13\Timeseries\MRItimeseries\csfRegressionResiduals_FB_S_13_kw3
and cluster map:
Y:\ClusterResults\kw3\CA0050\Cluster_nR10000_CA0050_kw3_fir2_fspread_015_fOrder_5_tw_15csfr_clusterMap.nii

Output: Cluster timeseries
Y:\Subjects\Subject13\Timeseries\PhasesAnalysis\ClusterTimeseries_csfr_S_13_kw3_CA0050

07/11/2016 Commented the line that cut timeseries since it has to be done
later



%}
%% set output filename
clusterTimeseries_filename = global_filename(subj_idx,cfgMain,'clusterTimeseries_filename');
%% insideBrain
insideBrain = tools_getIndexBrain('inside');
%% Get cluster map
clusterMap = ft_read_mri (strcat(global_filename(0,cfgMain,'clusterOutputFilename'),'_clusterMap.nii'));
clusterMap = clusterMap.anatomy(:);
below_significance = find(clusterMap > 12); % only use first 12 clusters, which correspond to significant clusters
clusterMap (below_significance) = 0;
clusterMap = clusterMap(insideBrain);
%% Load BOLD timeseries
BOLDtimeseriesFilename = global_filename(subj_idx,cfgMain,'filename_csfr_Residuals_FB');
load(BOLDtimeseriesFilename)
% error_csf_z = error_csf_z(cfgMain.beginCut:cfgMain.endCut,:); 
%% get cluster timeseries
nCluster = max(clusterMap);

clusterTimeseries = zeros(450,nCluster);

for iCluster = 1:nCluster
    voxelsInCluster = find(clusterMap==iCluster);
clusterTimeseries(:,iCluster) = nanmean(error_csf_z(:,voxelsInCluster),2);
end

%% save
save(clusterTimeseries_filename,'clusterTimeseries')

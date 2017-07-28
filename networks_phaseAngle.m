%{

This scripts loads the angle of phase locking of each cluster and subject, average them, 
perform the watson-william statistical test, plot them, and store it in nifti format
(Used for Figure 4 of the paper)
IR 29/06/2017

Inputs:
Phase angle per cluster
    Y:\Subjects\Subject13\Timeseries\PhasesAnalysis\ClusterTimeseriesPhaseAngle_csfr_S_13_kw3_CA0050
Cluster map
    Y:\ClusterResults\kw3\CA0050\Cluster_nR10000_CA0050_kw3_fir2_fspread_015_fOrder_5_tw_15csfr_clusterMap.nii


Output:
    Angle x Cluster nifti image
Y:\ClusterResults\Angle\angleXClusterGasnet



%}
%% definitions

subjects = global_subjectList
cfgMain = global_getcfgmain

%% load the angle of each cluster and subject

angleClusters_allsubjects = zeros(30,12);
for iSubj = 1:length(subjects)
    subj_idx = subjects(iSubj)
    clusterTimeseries_phaseAngle_filename = global_filename(subj_idx,cfgMain,'clusterTimeseries_phaseAngle_filename');
    load(clusterTimeseries_phaseAngle_filename)
    angleClusters_allsubjects(iSubj,:) = anglePLVXCluster_aligned
end

group_angleCluster = mean (exp (1i* (angle(angleClusters_allsubjects)))); % average across subjects

% [values indexes] = sort(angle(group_angleCluster),'descend') % sort them according to 

    colorTable = ...
[0 1 0;0 0.9 0.1; 0 0.8 0.2; 0 0.7 0.3; 0 0.6 0.4; 0 0.55 0.45; 0 0.45 0.55 ; 0 0.4 0.6;0 0.3 0.7; 0 0.2 0.8; 0 0.1 0.9; 0 0 1 ]
labels = {'vOcc','dOcc','RCZp','SIr','CCZ','dPrec','PCC','SII l','OccTem','SII r','la dPrec','IPS',};

labels = labels(indexes)
% colorTable = colorTable(indexes,:)

%% Statistical test

magnitudeClusters_allsubjects = abs(mean (exp (1i* (angle(angleClusters_allsubjects)))));
angle_Clusters_allsubjects = angle(mean (exp (1i* (angle(angleClusters_allsubjects)))));
data2test = angle(angleClusters_allsubjects)
indexes = repmat([1:12],30)
indexes = indexes (:,1:12)
[p table]= circ_wwtest(data2test(:),indexes(:))

%% Plot polar plot with angles

figure
for iCluster = 1:12
    i=iCluster
    iCluster = indexes(iCluster)
h1 = polar([angle(group_angleCluster(iCluster));angle(group_angleCluster(iCluster))],[zeros(1,1);abs(group_angleCluster(iCluster))])
    hold on

set(h1,'color',colorTable(i,:),'linewidth',7)

end

lgs = legend(labels)
set(lgs,'position',[0.40 0.5 0.08 0.1])

th = findall(gcf,'Type','text');
% legend('Individual phase differences','Mean phase difference')
for i = 1:length(th),
set(th(i),'FontSize',16)
end
shg
title('mean AnglexCluster across subjects', 'FontSize', 20)


%% Write nifti file with the angle of each cluster

% load clustermap of group level
insideBrain = tools_getIndexBrain('inside');

clusterMap = ft_read_mri (strcat(global_filename(0,cfgMain,'clusterOutputFilename'),'_clusterMap.nii'));
clusterMap = clusterMap.anatomy(:);
below_significance = find(clusterMap > 12); % only first 12 clusters are significant
clusterMap (below_significance) = 0;
clusterMap = clusterMap(insideBrain);

filenameAngle4exportGasnet = strcat(global_path2root,'ClusterResults',filesep,'Angle',filesep,'angleXClusterGasnet')

anglePLVGasnet = zeros(1,71759); % 71759, size of vector with voxels inside the brain
for iCluster = 1:12
voxelsInCluster = find(clusterMap==iCluster);
anglePLVGasnet(voxelsInCluster)= angle(group_angleCluster(iCluster));
end
angle4export_all = zeros(1,153594); % size of whole volume

angle4export_all(insideBrain) = anglePLVGasnet;

tools_writeMri(angle4export_all,filenameAngle4exportGasnet)


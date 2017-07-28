%{
paper_randomization_figure

this scripts loads all the stat files output from paper_randomization_statsCluster
and plot the histogram used in figure 2C, the distribution of summary
statistics across randomization of EGG time shifts
%}
dataDir = [global_path2root 'ClusterResults' filesep 'Randomizations' filesep']

clusterStats  = dir( fullfile( dataDir,'RandomCluster_N*_summary.mat')); %# list all offset nii files
nFiles = size(clusterStats,1);
clusterStats={clusterStats.name};


n_sig_clusters = zeros(1,nFiles);
sum_t_sigCluster = zeros(1,nFiles);
sum_AbsT = zeros(1,nFiles);



for iFile = 1:nFiles
load([dataDir clusterStats{iFile}])
    
n_sig_clusters(iFile)= summaryResults.NSPClusters; 
sum_t_sigCluster(iFile) = summaryResults.sumTsigCluster;
sum_AbsT(iFile) = summaryResults.sumOfAbsT;
end

path2empiricalresults = 'D:\Physiens\ClusterResults\kw3\CA0050\Cluster_nR10000_CA0050_kw3_fir2_fspread_015_fOrder_5_tw_15csfr_summary.mat'
empirical = load(path2empiricalresults)
sum_AbsT(end+1) = empirical.summaryResults.sumOfAbsT;

figure
nhist(sum_AbsT,'maxx',12000);shg
title(['Sum of Abs T across timeshifts randomizations N = ' num2str(nFiles)],'FontSize',20)
h = gca
set(h,'FontSize',20)
xlabel('Sum of Abs T','FontSize',20)

%% New file new figure hist sig cluster

dataDir = 'D:\Physiens\ClusterResults\Randomizations\'

clusterStats  = dir( fullfile( dataDir,'RandomCluster_N*_stats.mat')); %# list all offset nii files
nFiles = size(clusterStats,1);
clusterStats={clusterStats.name};

n_sig_clusters = zeros(1,nFiles);
sum_t_sigCluster = zeros(1,nFiles);

path2empiricalresults = 'D:\Physiens\ClusterResults\kw3\CA0050\Cluster_nR10000_CA0050_kw3_fir2_fspread_015_fOrder_5_tw_15csfr_stats.mat'
empirical = load(path2empiricalresults)

for iFile = 1:nFiles
    iFile
    load([dataDir clusterStats{iFile}])

    
NClusters = length(stat.posclusters);%)fieldnames(stat.posclusters)
for iCluster = 1:NClusters
    indNPClustersSig(iCluster)  = stat.posclusters(1,iCluster).prob <= 0.025 ; % Find significant clusters in the stat structure
end
NPClusters = sum(indNPClustersSig);
sumTsigCluster = 0

if NPClusters >= 1
    
    for iSigCluster = 1:NPClusters
        
        sumTsigCluster  =  sumTsigCluster + stat.posclusters(1,iSigCluster).clusterstat;
    end
end

    
n_sig_clusters(iFile)= NPClusters; 
sum_t_sigCluster(iFile) = sumTsigCluster;
% sum_AbsT(iFile) = summaryResults.sumOfAbsT;
end

figure
nhist(sum_t_sigCluster(find(n_sig_clusters)),'maxx',2500);shg
title(['Sum(t) across randomizations with significant clusters N = ' num2str(sum((n_sig_clusters>=1)))],'FontSize',20)
h = gca
set(h,'FontSize',20)
xlabel('Sum(t)','FontSize',20)

sumTsigClusterEmpirical = round(sum([374.159819445886;285.067444924685;269.545526671763;216.545478953478;196.660269189218;186.615257795487;175.579587285756;167.040572272570;163.343974580993;146.078118639996;125.455949662034]))



find(sum_t_sigCluster>sumTsigClusterEmpirical)
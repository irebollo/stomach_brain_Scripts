function timeseries_statsCluster_Heart(subjects,cfgMain)

%{


Perform group level statistics comparing empirical versus chance PLV
by using the clustering randomization procedure provided by fieldtrip (Maris and Oostenveld
2007)

First it load the data of all subject (empirical and chance PLV) and puts
them into fieldtrip format. Thens it performs the statistics and then it
save into Physiens/ClusterResults a mask containign all significant voxels
(OutputFilename_mask),  each significant cluster (OutputFilename_Ncluster)
the tvalue in significant clusters (OutputFilename_tmap) and the tmap in
the whole brain (OutputFilename_WholeBrainTmap), it also saves a summary of
the statistics (OutputFilename_stats)

inputs:
subjects = subjectects in which the analysis will be performed -> global_subject list
cfgMain must contain fields
    kernelWidth,Timeseries2Regress,frequencySpread ,fOrder,beginCut,endCut
kernelWidth: with of the smoothing kernel from preprocessing, paper  =
3mm
cfgMain.Timeseries2Regress should be 'csf' to load residuals of csf
regression
fOrder : mutiplicative factor of the order of the filter
frequencySpread: spead of the time domain filter in hz * 1000, paper = 0.015 hz = 15,
begin and end cut are the voulmes that are discarded to avoid the filter
ringing artifact
cfgMain.transitionWidth is the transition width of the filter, paper is 15 
offset is with respect to EGG peaking filter, only for control analysis.
offset is in hz x 1000 e.g. and offset of 0.006 hz is a value of 6
clusterAlpha = first level threshold one-sided for determining candidate
clusters

numberofrandomizations  = number of times the labels empirical and
chance PLV will be switched. For this study we use 10000

Commented by IR 28/06/2017

Output: Y:\ClusterResults\kw3\CA0050\Cluster_nR10000_CA0050_kw3_fir2_fspread_015_fOrder_5_tw_15csfr
suffix _mask significant voxels exceeding monte carlo p > 0.025

%}


%% cfgMain parameters used in the script 


numrandomization =cfgMain.numberofrandomizations;
clusterAlpha = cfgMain.clusterAlpha;

%% Output filename
indInside =  tools_getIndexBrain('inside') ;
indOutside = tools_getIndexBrain('outside') ;
% Output filename
  clusterOutputFilename = global_filename(0,cfgMain,'HeartclusterOutputFilename');

%% Load data of all subjects

    empirical = zeros(length(subjects),153594); % Preallocate
    surrogate = empirical; % for calculating t value
    
for iS=1:length(subjects)
    subj_idx = subjects(iS);
    
    % empirical and chance PLV filenames
   

filenamePLV = strcat(global_path2root,'\Subjects\','Subject',sprintf('%.2d',subj_idx),'\Timeseries','\Heart\empPLV_heart_subject_',sprintf('%.2d',subj_idx),'.nii');

% filenamePLVSurrogate = strcat(global_filename(subj_idx,cfgMain,strcat('medianRotationFilename_',cfg.Timeseries2Regress)),'.nii');
filenamePLVSurrogate  = strcat(global_path2root,'\Subjects\','Subject',sprintf('%.2d',subj_idx),'\Timeseries','\Heart\surrPLV_heart_subject_',sprintf('%.2d',subj_idx),'.nii');

    % Load empirical PLV  
    
    PLVGroupEmpirical{iS} = ft_read_mri(filenamePLV); % Put into cell
    PLVGroupEmpirical{iS}.Nsubject = subjects(iS);
    
    % Preparing the FieldTrip structure needed for randomization     
    PLVGroupEmpirical{iS}.coh = PLVGroupEmpirical{iS}.anatomy; 
    PLVGroupEmpirical{iS} = rmfield(PLVGroupEmpirical{iS},'anatomy');
    PLVGroupEmpirical{iS}.inside = indInside;
    PLVGroupEmpirical{iS}.outside = indOutside;
    
    % Load surrogate PLV and prepare structure for surrogate PLV
      
    PLVGroupSurrogate{iS} = ft_read_mri(filenamePLVSurrogate);
    PLVGroupSurrogate{iS}.Nsubject = subjects(iS);
    PLVGroupSurrogate{iS}.coh = PLVGroupSurrogate{iS}.anatomy;
    PLVGroupSurrogate{iS} = rmfield(PLVGroupSurrogate{iS},'anatomy');
    PLVGroupSurrogate{iS}.inside = indInside;
    PLVGroupSurrogate{iS}.outside = indOutside;
    
    empirical(iS,:) = PLVGroupEmpirical{iS}.coh(:);
    surrogate(iS,:) = PLVGroupSurrogate{iS}.coh(:);
    
end


%% Run stats
% run statistics over subjects %
cfgStats=[];
cfgStats.dim         = PLVGroupEmpirical{1}.dim;
cfgStats.method      = 'montecarlo';
cfgStats.statistic   = 'ft_statfun_depsamplesT';
cfgStats.parameter   = 'coh';
cfgStats.correctm    = 'cluster';
cfgStats.numrandomization = numrandomization;
% cfgStats.alpha       = 0.05; % note that this only implies single-sided testing
cfgStats.alpha       = 0.025; % note that this only implies single-sided testing
cfgStats.clusteralpha = clusterAlpha;
cfgStats.tail        = 0;
cfgStats.inside = indInside;
cfgStats.outside = indOutside;


% con, the second condition is the median rotation PLV!
nsubj=numel(PLVGroupEmpirical);
cfgStats.design(1,:) = [1:nsubj 1:nsubj];
cfgStats.design(2,:) = [ones(1,nsubj) ones(1,nsubj)*2];
cfgStats.uvar        = 1; % row of design matrix that contains unit variable (in this case: subjects)
cfgStats.ivar        = 2; % row of design matrix that contains independent variable (the conditions)

stat = ft_sourcestatistics(cfgStats,PLVGroupEmpirical{:}, PLVGroupSurrogate{:});
% Actual call to the statistic function

%% Get indexes of clusters

NClusters = length(stat.posclusters);%)fieldnames(stat.posclusters)
for iCluster = 1:NClusters
    indNPClustersSig(iCluster)  = stat.posclusters(1,iCluster).prob <= cfgStats.alpha ; % Find significant clusters in the stat structure
end

NPClusters = sum(indNPClustersSig);

% write them into HDD
for iCluster = 1:NPClusters 
    data = zeros(size(stat.mask));
    indCluster = find(stat.posclusterslabelmat == iCluster);
    data(indCluster)  = iCluster*10 ;
    tools_writeMri(data,strcat(clusterOutputFilename,'_map_clusterN',num2str(iCluster)))
end

%% Sumary of results


if NPClusters >= 1
sumTsigCluster = stat.posclusters(1,1:NPClusters).clusterstat;
summaryResults.NSPClusters = NPClusters;
summaryResults.sumTsigCluster=sumTsigCluster;

else
  summaryResults.NSPClusters  = 0;
  summaryResults.sumTsigCluster= 0;

end


% Perform a ttest at every voxel to obtain a tmap
empiricalBrain = empirical(:,indInside);
surrogateBrain = surrogate (:,indInside);

[h,p,ci,statsTtest] = ttest(empiricalBrain,surrogateBrain);
sumOfAbsT = sum(abs(statsTtest.tstat));
sumOfPLV = sum(empiricalBrain(:)); % empirical PLV
sumOfPLVChance = sum(surrogateBrain(:)); % chance PLV

summaryResults.sumOfAbsT=sumOfAbsT;
summaryResults.sumOfPLV=sumOfPLV;
summaryResults.sumOfPLVChance=sumOfPLVChance;


tstatCluster = zeros(1,153594);
tstatCluster(stat.inside) =statsTtest.tstat;
tstatCluster = reshape(tstatCluster,53,63,46);
tstatCluster(stat.mask==0) = 0;


tXvoxel = zeros(1,153594);
tXvoxel(indInside) = statsTtest.tstat;
tXvoxel = reshape(tXvoxel,53,63,46);

%% Save and  write results
save(strcat(clusterOutputFilename,'_stats.mat'),'stat')
save(strcat(clusterOutputFilename,'_summary.mat'),'summaryResults')

% data = stat.posclusterslabelmat;
% data(stat.mask)=0;
% tools_writeMri(data,strcat(clusterOutputFilename,'_ClusterMap'))

data = stat.mask;
tools_writeMri(data,strcat(clusterOutputFilename,'_mask'))
tools_writeMri(tstatCluster,strcat(clusterOutputFilename,'_tmap'))
tools_writeMri(tXvoxel,strcat(clusterOutputFilename,'_WholeBrainTmap'))


end
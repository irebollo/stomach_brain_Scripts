function output = paper_MakeTable1(cfgMain)
%{
paper_MakeTable1
The output of this function is used in table 1 of the paper to identify the regions
belonging to the gastric network using automatic anatomical labeling(AAL).
Input:  cfgMain structure used to obtain results, which allows to load the group
level masks. The output is the area number, the number and volume of voxels
in that area, the percentage of the area in the cluster, the tvalue at peak
coordinates and the x y z MNI coordinates of the peak

%}
%% Load data

% use global_filename to load group level stats
statsFilename = [global_filename(0,cfgMain,'clusterOutputFilename') '_stats.mat'];
load(statsFilename)
% this will load the structure 'stat'


% Load atlas fieldtrip lookup
atlas = ft_read_atlas([global_path2root 'scripts4paper' filesep 'files' filesep 'ROI_MNI_V4.nii']);
labelsAAL = load ([global_path2root 'scripts4paper' filesep 'files' filesep 'labels_AAL.mat']);% load labels filename


% (stat = gastric network)
listVoxels =stat.mask(:); %mask of significant voxels
indListSignificantVoxels = find(listVoxels); %% all significant voxels

% preallocate
indAreasSignificantVoxels = zeros(length(indListSignificantVoxels),1); % area number of each significant voxel
labelAreasSignificantVoxels = cell(length(indListSignificantVoxels),1); % area name of each significant voxel


%% iterate through all significant voxels and obtain MNI coordinates, AAL label, index of label,
for iVoxel=1:length(indListSignificantVoxels)
    
    % ind2sub: get the coordinates of a matrix based on a linear index
    [x y z ] = ind2sub([53 63 46], indListSignificantVoxels(iVoxel));
    mniSpace = tools_vox2mni([x y z ]',1);
    [label ,indLabel]= tools_atlas_lookup(atlas,mniSpace' ,'queryrange',1 , 'inputcoord','mni');
    
    indAreasSignificantVoxels(iVoxel) = cell2mat(indLabel);
    labelAreasSignificantVoxels{iVoxel} = label;
end

clear x y z mniSpace label indLabel

% clusterMap: Which gastric network cluster belongs to which cluster
clusterMapFilename = [global_filename(0,cfgMain,'clusterOutputFilename') '_clusterMap.nii'];
clusterMap = ft_read_mri(clusterMapFilename);
clusterMap = clusterMap.anatomy(:); % vector

% only significant clusters
ClusterXVoxel = clusterMap(indListSignificantVoxels)

% stat.stat correspond to the t-value
stats_vector = stat.stat(indListSignificantVoxels);



%% iterate through the 12 significant clusters and make table as a function of the clusters

% preallocate output table
output_table = zeros (1,9);

for iCluster = 1:12
    clear nVoxels indVoxelsInCluster areasInCluster iAreaInCluster output_temp tInPeak percentageAreaInCluster MNIcoordinates nVoxelsIArea indexAreaInAtlas XPeak YPeak ZPeak xMNI yMNI zMNI
    
    % find voxels belonging to the cluster
    indVoxelsInCluster = find(ClusterXVoxel == iCluster)
    
    % make a list of the areas in the cluster, the call to unique is to
    % avoid repetitions
    areasInCluster = unique(indAreasSignificantVoxels(indVoxelsInCluster))
    
    
    %         Loop through areas in iCluster
    for iAreaInCluster = 1:length(areasInCluster)
        
        % number of voxels in that area
        nVoxels(iAreaInCluster,1) = sum(indAreasSignificantVoxels(indVoxelsInCluster) == areasInCluster(iAreaInCluster));
        
        
        % look for the max value of t in that area
        tInPeak (iAreaInCluster,1) =  max(stats_vector(indVoxelsInCluster...
            (find(indAreasSignificantVoxels(indVoxelsInCluster) == areasInCluster(iAreaInCluster))...
            )));  % this search in the stats_vector the voxels in the cluster that
        % belong to the current area
        indListSignificantVoxels(find(stats_vector == tInPeak(iAreaInCluster)))
        
        %in voxel space
        [XPeak, YPeak, ZPeak] = ind2sub([53 63 46], ...
            indListSignificantVoxels(find(stats_vector == tInPeak(iAreaInCluster))));
        
        % mni
        [MNIcoordinates] = tools_vox2mni([XPeak, YPeak, ZPeak]',1);
        xMNI(iAreaInCluster,1)=MNIcoordinates(1);
        yMNI (iAreaInCluster,1)=MNIcoordinates(2);
        zMNI  (iAreaInCluster,1)=MNIcoordinates(3);
        
        
        % get percentage area in cluster
        indexAreaInAtlas = find (atlas.tissue(:) == areasInCluster(iAreaInCluster)); % find all voxels in atlas that belong to that area

        nVoxelsIArea = length  (indexAreaInAtlas) /5.8767; % voxels in atlas(=91*109*91)/153594 (voxels in physiens) to account for voxel size difference between atlas and physiens results
        percentageAreaInCluster(iAreaInCluster) = (nVoxels(iAreaInCluster,1) * 100 )/ nVoxelsIArea ;
        
    end
    
    output_temp = [];
    output_temp(:,2) = areasInCluster;
    output_temp(1:length(areasInCluster),1) = iCluster;
    output_temp(:,3) = nVoxels;
    output_temp(:,4) = nVoxels *3;
    output_temp(:,5) = percentageAreaInCluster;
    output_temp(:,6) = tInPeak;
    output_temp(:,7) = xMNI;
    output_temp(:,8) = yMNI;
    output_temp(:,9) = zMNI;
    
    output_table = [output_table; output_temp] ;
end

% remove from the table unkwown areas (indArea =0)
output_table(find(output_table(:,2)==0),:) = [];
% remove from the table entries with less than 1% of percentageAreaInCluster
output_table(find(output_table(:,5)<1),:) = [];

output.table = output_table
output.labels = cell(length(output.table),1)
for iEntry = 1:length(output.table)
    output.labels{iEntry} = labelsAAL.labels{output.table(iEntry,2)}
end


end

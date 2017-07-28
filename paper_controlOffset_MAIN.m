%{
paper_controlOffsetMAIN

Perform filter offset control
Filters BOLD and EGG timeseries with frequency offsets from the EGGpeak in
steps of 0.001 Hz, calculates empirical and chance PLV and performs group
level statistics.
It then loads the group level stats of all offsets and plots them to be
used in figure 2B

Modifying the parameter cfgMain.offset changes the filenames of all relevant files when calling global_filename
in order to reflect the desired filter offset.

Output is the figure 2B of the paper

IR 04/07/2017
%}


%% 1 perform the gastric network analysis with all filter offsets
cfgMain = global_getcfgmain;
subjects= global_subjectList;

offsetList = [-7 -5 -3 -1 1 3 5 7]
% divided by 1000 i.e -0.007 hz, -0.005 hz, -0.003 hz, etc 

% iterates through offsets
for jOffset = 1:length(offsetList)
    cfgMain.offset = offsetList(jOffset)
    % iterates through subjects
    for iSubj = 1: length(subjects)
    close all
    subj_idx = subjects(iSubj)
        prepro_egg(subj_idx,cfgMain)
        timeseries_preparePhases_Regression(subj_idx,cfgMain) % Filter and hilbert transform residuals of csf regression
        timeseries_mapPLV_Regression(subj_idx,cfgMain) % Obtain PLV per voxel
        timeseries_medianRotation_Regression(subj_idx,cfgMain) % Obtain chance PLV per voxel  
    end
    timeseries_statsCluster_Regression(subjects,cfgMain)
end

%% Plot

offsetlist = [-7 -5 -3 -1 0 1 3 5 7 ] % 0 offset is added so is also plotted
nOffsets = length(offsetlist)


sumSigT = zeros(1,nOffsets);
sumAllT = zeros(1,nOffsets);

plotDir = strcat (global_path2root,'Figures',filesep);
plotFilenameSigT = strcat(plotDir,'offsetControl_sigT');
plotFilenameAbsT = strcat(plotDir,'offsetControl_absT');

% iterate through offsets and load relevant data to plot
for iOffset = 1:nOffsets
    
    currentOffset = offsetlist(iOffset);
    cfgMain.offset = currentOffset


   % load group level stats output from the correpsonding offset
    load(strcat(global_filename(0,cfgMain,'clusterOutputFilename'),'_stats'))
    
   % number of positive clusters (empircal PLV larger than chance PLV)
   nClusters = length(stat.posclusters);
        
   
   %find which clusters are significant
   %initialize
   ClustersSignificant = zeros(nClusters,1);
   allClustersT = zeros(nClusters,1);
   for iCluster = 1:nClusters    
                allClustersT(iCluster) = stat.posclusters(iCluster).clusterstat;
        ClustersSignificant(iCluster) = stat.posclusters(iCluster).prob < 0.026; % i.e cluster stat sig at a 0.025 (one-sided) threshold
   end
    
   % sum of t in significant clusters
    sumSigT(iOffset) = sum(allClustersT(logical(ClustersSignificant)));
    
% Not used in paper; sum t across the whole brain
current_offset_T = ft_read_mri(strcat(global_filename(0,cfgMain,'clusterOutputFilename'),'_WholeBrainTmap.nii'));
sumAllT(iOffset) = sum(abs(current_offset_T.anatomy(:)));

end


labelsOffset = {}
for iOffset = 1:nOffsets
labelsOffset{iOffset} = num2str(offsetlist(iOffset)/1000)
end


% t in significant regions
figure; 
sigTPlot = bar(sumSigT,'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7],'LineWidth',5)
h = gca; 
set(h,'FontSize',20)
title(['KW' 32 num2str(cfgMain.kernelWidth) 32 'CA' 32 num2str(cfgMain.clusterAlpha)...
    32 'sum of t in significant clusters'],'Fontsize',30)
xlabel('offset','Fontsize',25)
ylabel('sum of t','Fontsize',25)
title(['sum of t in significant clusters as a function of filter offset'],'Fontsize',30)
numberOfXTicks = length(labelsOffset);
xData = get(sigTPlot,'XData');
set(gca,'Xtick',linspace(xData(1),xData(end),numberOfXTicks))
set(h,'XTickLabel',labelsOffset)
   
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    set(gcf, 'PaperPositionMode', 'auto');
    
    print ('-dpng', '-painters', eval('plotFilenameSigT'))
    print ('-depsc2', '-painters', eval('plotFilenameSigT'))
    saveas(sigTPlot,strcat(plotFilenameSigT,'.fig'))
    


% abs t whole brain
figure; 
absTplot = bar(sumAllT,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',5)
h = gca; 
set(h,'FontSize',20)
title(['KW' 32 num2str(cfgMain.kernelWidth) 32 'CA' 32 num2str(cfgMain.clusterAlpha)...
    32 'sum of t in significant clusters'],'Fontsize',30)
xlabel('offset','Fontsize',25)
ylabel('sum of abs t','Fontsize',25)
title(['sum of abs t as a function of filter offset'],'Fontsize',30)
numberOfXTicks = length(labelsOffset);
xData = get(absTplot,'XData');
set(gca,'Xtick',linspace(xData(1),xData(end),numberOfXTicks))
set(h,'XTickLabel',labelsOffset)
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    set(gcf, 'PaperPositionMode', 'auto');
    
    print ('-dpng', '-painters', eval('plotFilenameAbsT'))
    print ('-depsc2', '-painters', eval('plotFilenameAbsT'))
    saveas(absTplot,strcat(plotFilenameAbsT,'.fig'))

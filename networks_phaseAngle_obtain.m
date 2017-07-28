function networks_phaseAngle_obtain (subj_idx,cfgMain)
%{
This function 
1 Loads cluster timeseries for each subject and EGG
2 Get angle and value of phase locking for each cluster
3 Substract mean phase across clusters to the angle of each cluster
4 Plots previous steps and store them into HDD
5 Stores the angle and magnitude of PLV for each cluster

cfgMain filter parameters are applied here

Inputs
Cluster timeseries
    Y:\Subjects\Subject13\Timeseries\PhasesAnalysis\ClusterTimeseries_csfr_S_13_kw3_CA0050
EGG phase per volume
    Y:\Subjects\Subject13\Timeseries\EGGtimeseries\PhaseXvolume_S_13_fir2_fspread_015_ord_5_tw_15

Outputs:
PLV angle per cluster
    Y:\Subjects\Subject13\Timeseries\PhasesAnalysis\ClusterTimeseriesPhaseAngle_csfr_S_13_kw3_CA0050

IR 07/11/2016

%}
%% Definitions

% Load timeseries
clusterTimeseries_filename = global_filename(subj_idx,cfgMain,'clusterTimeseries_filename');
load(clusterTimeseries_filename)
EGGPhaseXVolumeFilename = global_filename(subj_idx,cfgMain,'EGGPhaseXVolumeFilename');
load(EGGPhaseXVolumeFilename)

% Define output filename
clusterTimeseries_phaseAngle_filename = global_filename(subj_idx,cfgMain,'clusterTimeseries_phaseAngle_filename');

% EGG peak info
peaksAllsubjects = global_getEGGpeaks;
indPeak = find (peaksAllsubjects(:,1) == subj_idx);
mostPowerfullFrequency = peaksAllsubjects(indPeak,3);

% Sanity plots
plotDir = strcat (global_path2subject(subj_idx),'PreprocessingLog',filesep);
plotFilename = strcat(plotDir,'S_',sprintf('%.2d',subj_idx),'_AngleClusters');
%% Do


% Badnpass Filter BOLD timeseries 
centerFrequency = mostPowerfullFrequency; %
filter_frequency_spread=cfgMain.frequencySpread/1000; % In hz
sr = 0.5; % 1 TR = 2s
filterOrder=(cfgMain.fOrder*fix(sr/(centerFrequency-filter_frequency_spread))-1);%in nsamples
transition_width= cfgMain.transitionWidth/100; % in normalised units
filteredClusters=tools_bpFilter(clusterTimeseries,sr,filterOrder,centerFrequency,filter_frequency_spread,transition_width,cfgMain.filterType);
% Hilbert
phaseClusters = hilbert(filteredClusters(cfgMain.beginCut:cfgMain.endCut,:));

% PLV and angle
anglePLVXCluster = mean (exp (1i* (bsxfun (@minus , angle(phaseClusters), angle (phaseXVolume)'))));

mean_angleAcrossClusters = mean (exp (1i* (angle(anglePLVXCluster))));


% We substract from each cluster angle the mean angle across all clusters

anglePLVXCluster_aligned = zeros(1,12);
for iCluster = 1:12
anglePLVXCluster_aligned (iCluster) = ...
    exp (1i* (angle(anglePLVXCluster(iCluster)) - angle(mean_angleAcrossClusters))) ;
end 

%% Save
save(clusterTimeseries_phaseAngle_filename,'anglePLVXCluster_aligned')

%% Sanity plots
if cfgMain.savePlots == 1
       
    if cfgMain.plotFigures == 0;
        SanityPlot = figure('visible','off');
    else
        SanityPlot = figure('visible','on');
    end
    
    
    colorTable = ...
[0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 0.5 0.5 0.5; 1 1 0; 0.3 0.6 1; 0.5 0.3 1; 0.5 0 0; 1 0.4 0.2; 0 0 0.5];
labels = {'Mean','vOcc','dOcc','RCZp','SIr','CCZ','dPrec','PCC','SIIl','Fus','SIIr','dPrecLeftAnt','IPS',};

subplot(1,2,1)

h2 = polar([angle(mean_angleAcrossClusters);angle(mean_angleAcrossClusters)],[zeros(1,1);abs(mean_angleAcrossClusters)])
set(h2,'color',[0 0 0],'linewidth',14)
hold on
for iCluster = 1:12
    
h1 = polar([angle(anglePLVXCluster(iCluster));angle(anglePLVXCluster(iCluster))],[zeros(1,1);abs(anglePLVXCluster(iCluster))])

set(h1,'color',colorTable(iCluster,:),'linewidth',7)

end
lgs = legend(labels)
set(lgs,'position',[0.1 0.15 0.05 0.05])

th = findall(gcf,'Type','text');
% legend('Individual phase differences','Mean phase difference')
for i = 1:length(th),
set(th(i),'FontSize',16)
end
shg
title(['AnglexCluster Notaligned subjects ' num2str(subj_idx)], 'FontSize', 20)



% figure;
subplot(1,2,2)

for iCluster = 1:12
  
h1 = polar([angle(anglePLVXCluster_aligned(iCluster));angle(anglePLVXCluster_aligned(iCluster))],[zeros(1,1);1])
hold on

set(h1,'color',colorTable(iCluster,:),'linewidth',7)

end
% legend(labels{2:end})
th = findall(gcf,'Type','text');
% legend('Individual phase differences','Mean phase difference')
for i = 1:length(th),
set(th(i),'FontSize',16)
end
shg
title(['AnglexCluster aligned subjects ' num2str(subj_idx)], 'FontSize', 20)


    
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    set(gcf, 'PaperPositionMode', 'auto');
    
    print ('-dpng', '-painters', eval('plotFilename'))
    print ('-depsc2', '-painters', eval('plotFilename'))
    saveas(SanityPlot,strcat(plotFilename,'.fig'))

end
end
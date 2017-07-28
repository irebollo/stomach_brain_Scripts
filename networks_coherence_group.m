%{
 
This Script loads the BOLD-EGG coherence values of all subjects in a loop
and average them at gastric peaking frequency separatly for each subject.
The resulting values are used in table 1 of the paper

Input: Each subject coherence between clusters
Y:\Subjects\Subject13\Timeseries\PhasesAnalysis\ClusterTimeseriesCoherence_csfr_S_13_kw3_CA0050

Output:
Descriptive statistic of group average coherence values

IR 28/06/2017
%}

%% Define
% Load EGG peaks info, average with respect to peak
subjects = global_subjectList
peaksAllsubjects = global_getEGGpeaks;
%% Execute

CoherenceXSubjectsXCluster = zeros(length(subjects),13,93); % 13 =12 clusters + global, 93 = freqs bins
CoherenceXSubjectsXCluster_peakEGG = zeros(length(subjects),13,61); % 13 = 12 clusters + global, 60 = freqs bins

for iSubj = 1:length(subjects)
    subj_idx = subjects(iSubj)
    %load
    clusterTimeseries_coherence_filename = global_filename(subj_idx,cfgMain,'clusterTimeseries_coherence_filename');
    load(clusterTimeseries_coherence_filename)
    CoherenceXSubjectsXCluster(iSubj,:,:) = coherence.cohspctrm(1:13,:);
    
    indPeaks(iSubj) = find(coherence.freq == peaksAllsubjects(iSubj,3)); % find the frequency bin corrdesponding to the EGG peak of the subject
    CoherenceXSubjectsXCluster_peakEGG(iSubj,:,:) = coherence.cohspctrm(1:13,indPeaks(iSubj)-30:indPeaks(iSubj)+30);
end

%% Average them separatly in absoulte frequency or relative frequency with respecto to EGG peak (frequency bin 31)

meanCoherenceXCluster_peak = squeeze(mean(CoherenceXSubjectsXCluster_peakEGG(:,:,31)));
serCoherenceXCluster_peak = squeeze(std (CoherenceXSubjectsXCluster_peakEGG(:,:,31))./sqrt(length(subjects)));

meanSquaredCoherenceXCluster_peak = squeeze(mean(CoherenceXSubjectsXCluster_peakEGG(:,:,31).*CoherenceXSubjectsXCluster_peakEGG(:,:,31)));
serSquaredCoherenceXCluster_peak = squeeze(std (CoherenceXSubjectsXCluster_peakEGG(:,:,31).*CoherenceXSubjectsXCluster_peakEGG(:,:,31))./sqrt(length(subjects)));

minSquaredCoherenceXCluster_peak = squeeze(min(CoherenceXSubjectsXCluster_peakEGG(:,:,31).*CoherenceXSubjectsXCluster_peakEGG(:,:,31)));

maxSquaredCoherenceXCluster_peak = squeeze(max(CoherenceXSubjectsXCluster_peakEGG(:,:,31).*CoherenceXSubjectsXCluster_peakEGG(:,:,31)));

SquaredCoherenceXClusterXSubject_peak = squeeze(CoherenceXSubjectsXCluster_peakEGG(:,:,31).*CoherenceXSubjectsXCluster_peakEGG(:,:,31));
meanSquaredCoherenceXSubject_peak = squeeze(mean(CoherenceXSubjectsXCluster_peakEGG(:,2:end,31).*CoherenceXSubjectsXCluster_peakEGG(:,2:end,31),2));% 1 is global signal

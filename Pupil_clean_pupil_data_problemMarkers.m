function Pupil_clean_pupil_data(subj_idx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Based on GUTSEE_SCRIPT_CLEAN_EYE_DATA
%%%
%%% Cleans eye data from blinks and saccades by interpolating windows
%%% containing artifacts
%%%
%%%
%%%
%%% Author: Nicolai Wolpert
%%% Email: Nicolai.Wolpert@ens.fr
%%% Version: 26/05/2017
%%% Modified by Ignacio Rebollo for Pupil data on Physiens 28/02/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% location of the GUTSEE server
% path2root = strcat(['Y:' filesep 'MEG' filesep]);

% addpath(genpath(strcat([GUTSEEMEG_location 'Scripts'])))

% %%% Add Fieldtrip to the path
% addpath('C:\Program Files\MATLAB\fieldtrip-20170315');
% ft_defaults;

% enter numbers of all subjects
% subjects = global_subjectList;

% nsubjects = length(subjects);

% work_dir = strcat([GUTSEEMEG_location 'Data_work' filesep]);

% reject trials from averaging with saccades larger than this threshold
thresh_sacc = 1.5;

% Eyelink sampling frequency
eyelink_fs = 1000;

% choose padding window in samples for blink saccade artifacts
padding_eyeblinks_seconds = 0.1;
padding_eyeblinks_seconds_after = 0.4;

padding_eyeblinks_samples = padding_eyeblinks_seconds/(1/eyelink_fs);

% specidfy the maximum length of a blink artfiact window to interpolate
maximum_artwindow_interpolation = 1.1121;

% samples of pupil data wher ethe derivative is larger than this value will
% be considered as artifacts
% thres_pupil_derivative = 0.1;
% thres_pupil_derivative = 10;

plotDir = strcat (global_path2subject(subj_idx),'PreprocessingLog',filesep);
%% Clean
% global_filename
% for iSubject=1:length(subjects)
    
    %% Load data
%     subj_idx = subjects(iSubject);
    
    fprintf(['\nSubject ' subj_idx '...\n\n']);
    
    events_eyelink = pupil_prepare_eyeklink_events(['Y:\Subjects\Subject' sprintf('%.2d',subj_idx) '\Eyelink\with MRI\Sujet'  sprintf('%.2d',subj_idx) 'i.asc']);
        [data_pupil,markers_raw] = prepro_egg_loadData(subj_idx,'Pupil',1); % Only loads those channel with label 'EGG'
        events_eyelink_benoit = readEyelink(['Y:\Subjects\Subject' sprintf('%.2d',subj_idx) '\eyelink\with MRI\Sujet'  sprintf('%.2d',subj_idx) 'i.asc']);
    
    
    
    %% Clean Pupil
    %     nBlocks = length(events_original);
    close all;
    data_pupil_clean = cell(1,1);
    
    
    % note samples corresponding to saccades larger 1.5 degree
    
    sacc_artf = [];
    
    
    amplitudes_saccades = [events_eyelink.esacc.saccamplitude];
    
    for isacc=1:length(events_eyelink.esacc)
        if amplitudes_saccades(isacc)>thresh_sacc
            samples_sacc = [events_eyelink.esacc(isacc).samplebegin events_eyelink.esacc(isacc).sampleend];
            sacc_artf = [sacc_artf; samples_sacc];
        end
    end
    

    
    blinks_artdef = []
    for iBlink=1:length(events_eyelink.eblink)
        samples_blink = [events_eyelink.eblink(iBlink).samplebegin events_eyelink.eblink(iBlink).sampleend];
        blinks_artdef = [blinks_artdef; samples_blink];
    end
    
    
    
    artdetectionplot = figure; 
    plot(events_eyelink.dat(4,:))
    for i=1:size(blinks_artdef, 1)
        hold on;
        vline(blinks_artdef(i, 1), 'r')
        hold on; vline(blinks_artdef(i, 2), 'r');      
    end

    for i=1:size(sacc_artf, 1) 
        hold on;
        vline(sacc_artf(i, 1), 'k')
        hold on; vline(sacc_artf(i, 2), 'k');
    end
    title(['s' num2str(subj_idx) 'rBlink kSacc'])
    
    
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    set(gcf, 'PaperPositionMode', 'auto');
    artdetectionplotFilename = strcat(plotDir,'S_',sprintf('%.2d',subj_idx),'_PupilArtDetect');
    print ('-dpng', '-painters', eval('artdetectionplotFilename'))
    print ('-depsc2', '-painters', eval('artdetectionplotFilename'))
    saveas(artdetectionplot,strcat(artdetectionplotFilename,'.fig'))
    
    %% Put data in fieldtrip strucutre
    
    
    timeseries2Coherence = events_eyelink.dat(4,:);
    load(strcat(global_path2root,'\scripts4paper\files\sampleFieldtripStruc.mat'))
    labelsChannelsMAIN = {'Pupil'};
    labelsChannels = labelsChannelsMAIN;
    clusterRegionsComparisons = timeseries2Coherence;
    dataPupilFT.hdr = EGG_downsampled.hdr;
    dataPupilFT.fsample = 1000;
    dataPupilFT.time{1,1}  = linspace(0,length(clusterRegionsComparisons)/dataPupilFT.fsample,length(timeseries2Coherence));
    dataPupilFT.label = labelsChannels;%channelStr;
    dataPupilFT.cfg = EGG_downsampled.cfg;
    dataPupilFT.trial{1,1} = clusterRegionsComparisons;
    dataPupilFT.sampleinfo = [1 length(timeseries2Coherence)]
    
    % rows artifact number
    % column 1 begin, second column end
    
    % replace artifacts with nans
    %     cfg = [];
    %     cfg.reject                      = 'nan';
    %     cfg.artfctdef.zvalue.artifact   = [sacc_artf; blinks_artdef];
    %     data_pupil_clean    = ft_rejectartifact(cfg,dataPupilFT);
%     
%     figure;
%     plot(data_pupil_clean.trial{1,1})
%     
    
    
    artdef_sacc_padded = [];
    for iart=1:length(sacc_artf)
        % add some padding
        artdef_sacc_padded = [artdef_sacc_padded; sacc_artf(iart,1)-padding_eyeblinks_samples sacc_artf(iart,2)+padding_eyeblinks_seconds_after];
    end
    
    artdef_blinks_padded = [];
    for iart=1:length(blinks_artdef)
        % add some padding
        artdef_blinks_padded = [artdef_blinks_padded; blinks_artdef(iart,1)-padding_eyeblinks_samples blinks_artdef(iart,2)+padding_eyeblinks_seconds_after];
    end
    
    
    
% merge blinks and saccades artifacts and sort them by onset

art_def_blinksSacc = [sacc_artf; blinks_artdef]
art_def_blinksSacc=sortrows(art_def_blinksSacc)
    
%     % get rid of overlaps
    iart = 1;
    while iart<length(art_def_blinksSacc(:,1))
        iart
        % scan for overlapping artifacts
        iart2 = iart+1;
        while iart2<length(art_def_blinksSacc(:,1))
            
            iart2
            if ~isempty(intersect([art_def_blinksSacc(iart,1):art_def_blinksSacc(iart,2)], [art_def_blinksSacc(iart+1,1):art_def_blinksSacc(iart+1,2)]))
                
                artbegin = min([art_def_blinksSacc(iart,1) art_def_blinksSacc(iart+1,1)]);
                artend = max([art_def_blinksSacc(iart,2) art_def_blinksSacc(iart+1,2)]);
                
                art_def_blinksSacc(iart,1) = artbegin;
                art_def_blinksSacc(iart,2) = artend;
                art_def_blinksSacc(iart+1, :) = [];
            end
            iart2 = iart2+1;
        end
        iart=iart+1;
    end
    
       
        
    % merge artifact windows seperated by less than 200 ms
    iart = 1;
    while iart<length(art_def_blinksSacc(:,1))
        
        if art_def_blinksSacc(iart+1, 1)-art_def_blinksSacc(iart, 2) < 200
            
            art_def_blinksSacc(iart, 2) = art_def_blinksSacc(iart+1, 2);
            art_def_blinksSacc(iart+1, :) = [];
            
        else
            iart = iart+1;
        end
    end
    % sample of first artifact might be negative due to padding
    if ~isempty(find(art_def_blinksSacc<0))
        art_def_blinksSacc(find(art_def_blinksSacc<0))=1;
    end
    
    % replace artifacts with nans
    cfg = [];
    cfg.reject                      = 'nan';
    cfg.artfctdef.zvalue.artifact   = art_def_blinksSacc;
    data_pupil_clean    = ft_rejectartifact(cfg,dataPupilFT);
 

% interpolate NaNs, but only if certain threshold for artifact wiondow
% length is not exceeded

data_pupil_clean_gaps = data_pupil_clean;
data_pupil_clean.trial{1} = fixgaps(data_pupil_clean.trial{1,1});

cleandataPlot = figure;
plot(data_pupil_clean.trial{1});
title(['s' num2str(subj_idx) 'Pupil clean'])

   set(gcf,'units','normalized','outerposition',[0 0 1 1])
    set(gcf, 'PaperPositionMode', 'auto');
    cleandataPlotFilename = strcat(plotDir,'S_',sprintf('%.2d',subj_idx),'_PupilCleanFB');
    print ('-dpng', '-painters', eval('cleandataPlotFilename'))
    print ('-depsc2', '-painters', eval('cleandataPlotFilename'))
    saveas(cleandataPlot,strcat(cleandataPlotFilename,'.fig'))
    
% ind_windows_long_saccades = find(([sacc_artf(:,2)-sacc_artf(:,1)]/eyelink_fs)>maximum_artwindow_interpolation);
% for iart=1:length(ind_windows_long_saccades)
%     data_pupil_clean.trial{1}(sacc_artf(ind_windows_long_saccades(iart),1):sacc_artf(ind_windows_long_saccades(iart),2)) = nan;
% end
% ind_windows_long_blinks = find(([blinks_artdef(:,2)-blinks_artdef(:,1)]/eyelink_fs)>maximum_artwindow_interpolation);
% for iart=1:length(ind_windows_long_blinks)
%     data_pupil_clean.trial{1}(blinks_artdef(ind_windows_long_blinks(iart),1):blinks_artdef(ind_windows_long_blinks(iart),2)) = nan;
% end
% 
% data_pupil_clean.trial{1} = fixgaps(data_pupil_clean.trial{1});



%% Cut the data at the begining and the end of the MRI acquisition
% 

indx_markers = find([events_eyelink.input.value]==253);
index_sample_begin= find(events_eyelink.dat(1,:) == [events_eyelink.input(indx_markers(5)).timestamp]+1000); % the pupil data is shifted 1s to account for the pupillry response delay -> neural signal of pupil dilation occurs 1s earlier than pupil peak dilation
index_sample_ends= find(events_eyelink.dat(1,:) == [events_eyelink.input(indx_markers(451)).timestamp])+2800; % takes 1.8s after the sample were the volume started + delay. 1.8 isntead of 2s due to resampling to 0.5, if I put 2s i end up with data from next volume (451 instead of 450)
% index_sample_ends= find(events_eyelink.dat(1,:) == [events_eyelink.input(indx_markers(447)).timestamp])+2800; % takes 1.8s after the sample were the volume started + delay. 1.8 isntead of 2s due to resampling to 0.5, if I put 2s i end up with data from next volume (451 instead of 450)


index_sample_begin_allVolumes = zeros(1,450);
for iVolume = 5:454
index_sample_begin_allVolumes(iVolume-4) = find(events_eyelink.dat(1,:) == [events_eyelink.input(indx_markers(iVolume)).timestamp]+1000);
end

figure;
plot(index_sample_begin_allVolumes,'ok')


figure;
plot([events_eyelink.input(indx_markers(1:450)).timestamp]+1000,'ok') ;

%% Take regressors  per volume

nansXvolume = zeros(1,448);
propotion_nansXvolume =zeros(1,448);
for iVolume = 5:452
nansXvolume(iVolume-4) = sum(isnan(data_pupil_clean_gaps.trial{1,1}(index_sample_begin_allVolumes(iVolume):index_sample_begin_allVolumes(iVolume)+1999)));
propotion_nansXvolume(iVolume-4) = sum(isnan(data_pupil_clean_gaps.trial{1,1}(index_sample_begin_allVolumes(iVolume):index_sample_begin_allVolumes(iVolume)+1999)))/2000; 
end

% total proportion of missing data during the scan (more than )
propotion_nan_wholeScan = sum(nansXvolume)/(448*2000);

%%
data_pupil_MRI = data_pupil_clean;
data_pupil_MRI.trial{1,1} = data_pupil_MRI.trial{1,1}(index_sample_begin:index_sample_ends) ;
data_pupil_MRI.time{1,1}=data_pupil_MRI.time{1,1}(index_sample_begin:index_sample_ends);
data_pupil_MRI.sampleinfo=[1 length(data_pupil_MRI.time{1,1})];


disp('Resampling...')
cfg = [];  %initialize configuration structure
cfg.detrend = 'no'; % remove linear trend from the data (done per trial)
cfg.demean = 'yes';
cfg.resamplefs= 0.5; % 4 x top-freq (15 cpm = 0.25 Hz) - Nyquist = 30 cpm  frequency at which the data will be resampled
Pupil_downsampled = ft_resampledata(cfg,data_pupil_MRI); % This procedure also lowpass filter the data at half the new sr
figure; plot(Pupil_downsampled.trial{1,1}); title(['Subject' subj_idx ' Pupil clean ']);


data_pupil_bp = ft_preproc_bandpassfilter (Pupil_downsampled.trial{1,1},0.5,[1/128 0.1]);
data_pupil_bp = ft_preproc_standardize (data_pupil_bp);

pupilRegressorPlot = figure;
subplot(2,1,1); plot(data_pupil_bp); title(['Subject' subj_idx ' Pupil clean ']);
title(['s' num2str(subj_idx) 'Pupil clean BP in scan zscored proportionArtefacts=' num2str(propotion_nan_wholeScan)])

pupil_clean_bp_derivative = diff(data_pupil_bp);
subplot(2,1,2)
plot(pupil_clean_bp_derivative)
% title('derivative')

   set(gcf,'units','normalized','outerposition',[0 0 1 1])
    set(gcf, 'PaperPositionMode', 'auto');
    pupilRegressorPlotFilename = strcat(plotDir,'S_',sprintf('%.2d',subj_idx),'_PupilRegressor');
    print ('-dpng', '-painters', eval('pupilRegressorPlotFilename'))
    print ('-depsc2', '-painters', eval('pupilRegressorPlotFilename'))
    saveas(pupilRegressorPlot,strcat(pupilRegressorPlotFilename,'.fig'))

pupil_data_bp_filename = ['Y:\Subjects\Subject' sprintf('%.2d',subj_idx) '\Timeseries\Pupil\Pupil_bp_MRI_'  sprintf('%.2d',subj_idx)];

save(pupil_data_bp_filename,'data_pupil_bp','pupil_clean_bp_derivative','nansXvolume','propotion_nansXvolume','propotion_nan_wholeScan')

end

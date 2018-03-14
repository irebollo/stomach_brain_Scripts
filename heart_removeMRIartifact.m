
%% parameters
cfgMain = global_getcfgmain;

%%

subjects = global_subjectList


 
for iSubj = 4:length(subjects)
    subj_idx = subjects(iSubj)

% EEG.data = double(EEG.data)

 [ALLEEG EEG CURRENTSET ALLCOM] = eeglab; % start EEGLAB under Matlab 

   EEG.etc.eeglabvers = '14.1.1'; % this tracks which version of EEGLAB is being used, you may ignore it
   
     subjectFolder= global_path2subject(subj_idx);
    brainampDir = strcat(subjectFolder,'brainamp',filesep,'with MRI',filesep);
    files= dir( fullfile( brainampDir,'*.eeg')); %# list all *.vhdr files
    filename = {files.name}';%'# file names
    FilenameBrainAmp = char(strcat(brainampDir,filename));
        
     EEG = pop_fileio(FilenameBrainAmp, 'channels',[5:7] );
%      EEG = pop_fileio(FilenameBrainAmp, 'channels',[5] );
     EEG.setname='ECG';
     EEG = eeg_checkset( EEG );
     EEG = eeg_checkset( EEG );
     EEG = pop_fmrib_fastr(EEG,0,10,30,'R128',0,1,0,0,0,0.03,[1:3],'auto');
%           EEG = pop_fmrib_fastr(EEG,0,10,30,'R128',0,1,0,0,0,0.03,[1],'auto');

     EEG = eeg_checkset( EEG );
     
     
     %% Cut data and keep only with fMRI on
     
     
for i=1:length(EEG.event)
fMRIvolume(i) = strcmp(EEG.event(i).type,'R128')
end

for i=1:length(EEG.event)
samplesStructure(i) = EEG.event(i).latency;
end

fMRIvolume_samples = samplesStructure(fMRIvolume);

dataECG = EEG.data;

dataECG_cutted =dataECG(:,fMRIvolume_samples(1):fMRIvolume_samples(450));

data_ECG_cutted_filename = global_filename(subj_idx,cfgMain,'data_ECG_cutted')
save(data_ECG_cutted_filename,'dataECG_cutted')

end
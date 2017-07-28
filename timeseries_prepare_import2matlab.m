function timeseries_prepare_import2matlab(subj_idx,cfg)

%{
Construct fMRI timeseries (time,voxels) from swaf images located at each subject folder
and saves them as a structure in the HDD.

Input files
hdr files from each volume
Y:\Subjects\Subject13\fMRI\acquisition1\RestingState\s3wafPHYSIENS_Sujet13-0002-00001-000001-01.hdr

Output 
concatenated MRI timeseries
Y:\Subjects\Subject13\Timeseries\MRItimeseries\fMRItimeseries_S13_kw3.mat

Ignacio Rebollo 16/3/2015
commented 28/06/2017

%}

% Get parameters from cfg structure

kernelWidth = cfg.kernelWidth;

% Import

dataDir = strcat(global_path2subject(subj_idx),'fMRI',filesep,'acquisition1',filesep,'RestingState',filesep);

output_filename = global_filename(subj_idx,cfg,'BOLDTimeseriesFilename')

BOLDtimeseries = []; % "raw" fmri data
BOLDtimeseries.fsample  = 1/2;  %0.5 hz, ~ TR=2s


% identify the files corresponding to the deisred spatial kernel smoothing
% obtained from the output from mri preprocessing 

if kernelWidth == 6
    filename= dir( fullfile( dataDir,'s6waf*.hdr')); %# list all *.hdr files of preprocessed images, this are the 450 volumes
elseif kernelWidth == 3
    filename= dir( fullfile( dataDir,'s3waf*.hdr')); %# list all *.hdr files of preprocessed images, this are the 450 volumes
elseif kernelWidth == 4
    filename= dir( fullfile( dataDir,'s4waf*.hdr')); %# list all *.hdr files of preprocessed images, this are the 450 volumes
elseif kernelWidth == 8
    filename= dir( fullfile( dataDir,'s8waf*.hdr')); %# list all *.hdr files of preprocessed images, this are the 450 volumes
elseif kernelWidth == 0
    filename= dir( fullfile( dataDir,'waf*.hdr')); %# list all *.hdr files of preprocessed images, this are the 450 volumes
elseif kernelWidth == 5
    filename= dir( fullfile( dataDir,'s5waf*.hdr')); %# list all *.hdr files of preprocessed images, this are the 450 volumes
end

% list all files
filename = {filename(~[filename.isdir]).name}';

for i=1:length(filename)
  filename{i} = fullfile(dataDir, filename{i});% Only interested in the name of the files
end

% data has to be stored from a 3D matrix to one big vector tmp.anatomy(:,:,:)
datVector = zeros(length(filename),153594); % data in vector format(time,x*y*z) 
% 153594 correpond to the number of voxels of physiens EPI sequences (3mm)

for i=1:length(filename) %iterates through all volumes, load them, and put them in matrix
  disp(i);
  tmp = ft_read_mri(filename{i}); % calls fieldtrip to read the mri intensity images of each volume
  datVector(i,:) = tmp.anatomy (:); % concatenate them
end

BOLDtimeseries.time  = [0:2:(length(filename)*2)-1]; % create time axis
BOLDtimeseries.trialVector = datVector;
BOLDtimeseries.transform=tmp.transform(:,:);%transformation matrix; 

save(output_filename,'BOLDtimeseries')
end
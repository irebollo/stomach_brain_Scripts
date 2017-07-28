function [output] = prepro_egg_path2data(subject,fileType)
%{
This funcion get the full path of different kind of files, e.g. brainamp markers, for a particular
subject. This function is called by prepro_egg_loadData

input:
subject = string, e.g. '01' '02' '03'
filtype= 'brainamp', 'brainampMarkers' 'fmri' , 'brainampWOIRM', 'brainampMarkersWOIRM'


example call:
pathtBrainAmp = getFilesPath('08','brainamp')

Output: string containing path to data

Commented IR 27/06/2017

%}
%%

subjectFolder= global_path2subject(subject);

% Check which kind of file to retrieve
brainamp=strcmp(fileType,'brainamp'); 
brainampMarkers=strcmp(fileType,'brainampMarkers'); 
brainampWOIRM=strcmp(fileType,'brainampWORIM'); 
brainampMarkersWOIRM=strcmp(fileType,'brainampMarkersWOIRM'); 
fmri=strcmp(fileType,'fmri'); 


if brainamp==1
    brainampDir = strcat(subjectFolder,'brainamp',filesep,'with MRI',filesep);
    files= dir( fullfile( brainampDir,'*.vhdr')); %# list all *.vhdr files
    filename = {files.name}';%'# file names
    output = char(strcat(brainampDir,filename));
end

if brainampMarkers==1
    brainampDir = strcat(subjectFolder,'brainamp',filesep,'with MRI',filesep);
    files= dir( fullfile( brainampDir,'*.vhdr')); %# list all *.vhdr header files
    filename = {files.name}';%'# file names
    output = char(strcat(brainampDir,filename));
end   

if brainampWOIRM==1
    brainampDir = strcat(subjectFolder,'brainamp',filesep,'without MRI',filesep);
    files= dir( fullfile( brainampDir,'*.vhdr')); %# list all *.vhdr files
    filename = {files.name}';%'# file names
    output = char(strcat(brainampDir,filename));
end 

if brainampMarkersWOIRM==1
    brainampDir = strcat(subjectFolder,'brainamp',filesep,'without MRI',filesep);
    files= dir( fullfile( brainampDir,'*.vhdr')); %# list all *.vhdr header files
    filename = {files.name}';%'# file names
    output = char(strcat(brainampDir,filename));
end

if fmri==1
    fmriDir=strcat(subjectFolder,'fMRI',filesep,'acquisition1',filesep,'RestingState');
    files= dir( fullfile( fmriDir,'swaf*.img')); %# list all swaf files
    filenames = {files.name}';%'# file names
    output = filenames;
end

end
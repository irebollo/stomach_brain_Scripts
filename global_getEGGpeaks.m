function EGGpeaks =  global_getEGGpeaks

%{
Load the file containing the channel number and exact frequency for all subjects during the fMRI acquisition,
which was manually estimated in the first call of prepro_egg and stored in 
the log file together with the EGG phase timeseries in each subject timeseries folder.
The file to be loaded  was created by the function prepro_egg_saveChannelInfo
Commented IR 27/06/2017

%}

load([global_path2root filesep 'scripts4paper' filesep 'files' filesep 'EGG_peaks_info'])
end
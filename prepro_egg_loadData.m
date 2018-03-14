function [dataInFieldtrip, markersInFieldtrip] = prepro_egg_loadData(subj_idx,dataType,irm)
%{

Function to retrieve subject EGG and IRM volume markers recorded inside or
outside the scanner. This function is called by prepro_egg
 
Input:
subj_idx= eg '8','12'
dataType= 'EGG' for electrogastrogram or 'EKG' for electrocardiogram
irm = if 1 will load data used inside the scanner, if 0 will load data from
outside the scanner

Output: Data and markers in fieldtrip format

IR 27/06/2017 Commented

%}

if irm==1 % data when the irm scanner was on
    pathToBrainAmp = prepro_egg_path2data(subj_idx,'brainamp');
    pathToMarkers = prepro_egg_path2data(subj_idx,'brainampMarkers');
    
    %Markers
    cfg = []; %structure de configuration
    cfg.dataset = pathToMarkers; % nom du fichier d'intérêt
    cfg.trialfun = 'ft_trialfun_general'; % fonction définissant les essais
    cfg.trialdef.eventtype = 'Response'; % type d'évènement
    cfg.trialdef.eventvalue = 'R128'; % valeur d'évènement
    cfg.trialdef.prestim = 1; % secondes écoulées avant l'évènement
    cfg.trialdef.poststim = 1; % secondes écoulées après l'évènement
    markersInFieldtrip = ft_definetrial(cfg);
    
elseif irm==0
    pathToBrainAmp = prepro_egg_path2data(subj_idx,'brainampWORIM');
    markersInFieldtrip = [];
end

EGG=strcmp(dataType,'EGG');
EKG=strcmp(dataType,'EKG');
Pupil= strcmp(dataType,'Pupil');


if EGG==1
        cfg = [];    %configuration structure
    cfg.dataset = pathToBrainAmp;
    cfg.channel = {'EGG_1', 'EGG_2', 'EGG_3', 'EGG_4'};
    dataInFieldtrip = ft_preprocessing(cfg); %structure with data (as one long continuous segment)
elseif EKG==1
    %Data
    cfg = [];    %configuration structure
    cfg.dataset = pathToBrainAmp;
    cfg.channel = {'ECG_6', 'ECG_7', 'ECG_8'};
    dataInFieldtrip = ft_preprocessing(cfg); %structure with data (as one long continuous segment)
elseif Pupil==1
    %Data
    cfg = [];    %configuration structure
    cfg.dataset = pathToBrainAmp;
    cfg.channel = {'D'};
    dataInFieldtrip = ft_preprocessing(cfg); %structure with data (as one long continuous segment)
  
end
end
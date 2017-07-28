function cfgMain = global_getcfgmain

%{

Function defining a matlab structure containing the final parameters used through the preprocessing and analysis pipeline. 
This structure is a mandatory input to most PHYSIENS functions, 
either to know which parameter apply or which filename retrieve/store.
In order to retrieve the path to a file, the different parts of the cfgMain
structure are concatenated when calling global_filename

IR 28/06/2017
%}

%Prepro
cfgMain.kernelWidth = 3; % HWHM of the smoothing kernel applied during fMRI preprocessing with SPM
cfgMain.Timeseries2Regress = 'csfr'; % timeseries to be regressed: csfr cerebro-spinal-fluid 

%Filter 
cfgMain.saveFiltered = 0; %XXXX ? save filtered what ?
cfgMain.filterType = 'fir2';
cfgMain.frequencySpread = 15; %in Hz divided by 1000 e.g. 15 = 0.015 Hz
cfgMain.beginCut = 16; % first and last 15 volumes from the total 450 are discarded systematically after bandpass filtering
cfgMain.endCut = 435;
cfgMain.fOrder = 5; % filter order multiplicative factor i.e. 5 = five times number of samples require to sample the slowest frequency
cfgMain.transitionWidth = 15; % in arbitrary units /100 e.g. 15 = 0.15. Bandpass filter parameter that controls the slope of the filter transition width used in tools_bpFilter


% Stats
cfgMain.numberofrandomizations = 10000; % number of randomizations of the clustering procedure (swtch labels empirical chance)
cfgMain.automaticChannelSelection = 0; % the first time is set to 0 then after calling prepro_egg_saveChannelInfo(cfgMain) is turn to 1
cfgMain.clusterAlpha = 0.005; % exact one-sided p value that will be used for the cluster forming threshold XXXX check


% Randomization control
cfgMain.randomized=0; % if 0 (default). if 1 uses a random EGG time-shift and perofrm analaysis for EGG control XXXX SCRIPT NAME

% Offset control
cfgMain.offset = 0; % filter offset from peak EGG in Hz /1000 e.g. 7 = 0.007

% Other
cfgMain.plotFigures= 1; % if 1 plots will be displayed if 0 they won't
cfgMain.savePlots = 1; % if 1 scripts store plots in disk if 0 not
cfgMain.scrubbed =0; % NOT IN USE
cfgMain.hrfconvolved = 0 % NOT IN USE. If one convolves the EGG with the HRF functions

end
%{
paper_quantifyAutonomicOverlap

This scripts quantify the overlap between voxels of  the gastric network and the
electrodermal activity related regions provided by Beisnner 2013.

Input: Y:\MasksFromOtherPapers\AutonomicMasks\EDA_clust_bin_3mm.hdr
Output: percentage of overlap in the command line
IR 05/07/2017
%}

% load gastric network
gasnet = global_getGastricNetwork;
ind_gasnet = find(gasnet);

% load EDA
EDA3mm = ft_read_mri([ global_path2root 'MasksFromOtherPapers' filesep 'AutonomicMasks' filesep 'EDA_clust_bin_3mm.hdr'])
EDA3mm.anatomy = round(EDA3mm.anatomy)

% calculate percentage overlap
(sum(EDA3mm.anatomy(ind_gasnet))*100)/length(ind_gasnet)
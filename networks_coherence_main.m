% Simple script that calls the relevant functions to perform the coherence analysis between
% gastric network and the EGG at the group level
% IR 28/06/2017

for iSubj = 1:length(subjects)
    subj_idx = subjects(iSubj)
   networks_coherence_getClusterTimeseries(subj_idx,cfgMain) 
   
 networks_coherence_getFullBandEGG(subj_idx,cfgMain) 

   
 networks_coherence_getCoherence(subj_idx,cfgMain) 

end
   
networks_coherence_group(cfgMain)

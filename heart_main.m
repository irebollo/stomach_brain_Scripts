heart_main

subjects = [8 9 13 15 16 17 19 20 22 25 26 29 31 32 34 35 36]

cfgMain = global_getcfgmain
for iSubj = 1:length(subjects)
    subj_idx = subjects(iSubj)
 heart_cohHRV_BOLD(subj_idx,cfgMain)   
end


%%

allS = zeros(length(subjects),153594);
for iSubj = 1:length(subjects)
    subj_idx = subjects(iSubj)
           HRV_BOLD_coherence_map_filename = global_filename(subj_idx,cfgMain,'HRV_BOLD_coherence_map');

temp = ft_read_mri([HRV_BOLD_coherence_map_filename '.nii']);
temp.flat =temp.anatomy(:);
allS(iSubj,:)=temp.flat;
end

emptyBrain = reshape(mean(allS),53,63,46);

tools_writeMri(emptyBrain,'meanSubjects')

temp.anatomy

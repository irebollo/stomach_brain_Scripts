
subjects = global_subjectList
cfgMain = global_getcfgmain


insideBrain = tools_getIndexBrain('inside');
gasnet = global_getGastricNetwork;
gasnet = gasnet(:);
gasnet = gasnet(insideBrain);

PLV_allsubjects = zeros(length(subjects),length(insideBrain));
chancePLV_allsubjects = zeros(length(subjects),length(insideBrain));
PPC_allsubjects = zeros(length(subjects),length(insideBrain));

for iSubj = 1:length(subjects)

    subj_idx = subjects(iSubj)
        
    PLV_filename = global_filename(subj_idx,cfgMain,'PLVXVoxelFilename_csfr');
    chancePLV_filename = global_filename(subj_idx,cfgMain,'medianRotationFilename_csfr');
    PPC_filename = global_filename(subj_idx,cfgMain,'PPC_filename');
    
    PLV = ft_read_mri([PLV_filename '.nii']);
    PLV = PLV.anatomy;
    PLV = PLV(insideBrain);
    
    chancePLV = ft_read_mri([chancePLV_filename '.nii']);
    chancePLV = chancePLV.anatomy;
    chancePLV = chancePLV(insideBrain);
    
    PPC = ft_read_mri([PPC_filename '.nii']);
    PPC = PPC.anatomy;
    PPC = PPC(insideBrain);
    
    
    PLV_allsubjects(iSubj,:) =PLV;
    chancePLV_allsubjects(iSubj,:) =chancePLV;
    PPC_allsubjects(iSubj,:) =PPC;
    
end


mean_PLV = mean(PLV_allsubjects);
mean_chancePLV = mean(chancePLV_allsubjects);
mean_PPC = mean(PPC_allsubjects);


mean_PLV_gasnet = mean(PLV_allsubjects(:,gasnet));
mean_chancePLV_gasnet = mean(chancePLV_allsubjects(:,gasnet));
mean_PPC_gasnet = mean(PPC_allsubjects(:,gasnet));

figure
plot(mean_PLV,mean_PPC,'ok')

figure
plot(mean_PLV.*mean_PLV,mean_PPC,'ok')

figure
plot(mean_PLV-mean_chancePLV,mean_PPC,'ok')

[r p] = corrcoef(mean_PLV,mean_PPC)


figure
plot(mean_PLV_gasnet,mean_PPC_gasnet,'ok')
[r p] = corrcoef(mean_PLV_gasnet,mean_PPC_gasnet)

%% Correlation each subject

corrallsubjects = zeros(1,30)
for iSubj =1:30
[r p] = corrcoef(PLV_allsubjects(iSubj,:),PPC_allsubjects(iSubj,:))
corrallsubjects(iSubj) = r(3) 
end

mean(corrallsubjects)

fisherz(corrallsubjects)
[h p ci stats]=ttest(fisherz(corrallsubjects))

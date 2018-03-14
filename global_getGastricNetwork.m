function [gasnet] = global_getGastricNetwork
% retrieves the final mask of the gastric network obtained from PHYSIENS 30
% subjects 
% 05/10/2016

gasnetRaw = ft_read_mri ([global_path2root filesep 'ClusterResults' filesep 'kw3' filesep...
    'CA0050',filesep,'Cluster_nR10000_CA0050_kw3_fir2_fspread_015_fOrder_5_tw_15csfr_mask_2LVLP0_0250_oneSIDED.nii']);
gasnet = logical(gasnetRaw.anatomy);

end
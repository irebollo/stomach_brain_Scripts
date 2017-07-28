function    [output] = tools_getIndexBrain(insideOutside)

% Get index of voxel outside/inside brain in vector format based on SPM
% apriori mask

% insideOutside= string saying if the index required is 'inside'or 'outside'

% IR commented 28/06/2017

mask = ft_read_mri (strcat(global_path2root,filesep,'scripts4paper',filesep,'files',filesep,'SPM_mask_THRESHOLDED_physiens.hdr'));

mask.anatomyVector=logical(mask.anatomy(:));

inside = strcmp('inside',insideOutside);
outside = strcmp('outside',insideOutside);

if inside == 1
    output = find (mask.anatomyVector == 1); % find values where there is a 1, therefore brain
elseif outside == 1
    output = find (mask.anatomyVector == 0); % find values where there is a zero, therefore no brain brain
end

end
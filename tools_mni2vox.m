function VOX= tools_mni2vox(MNI,transform)
% Transformation of MNI coordiantes to voxelspace

% MNI has to be a column vector
%function VOX=mni2vox(MNIcoordinates,transform)
% if transform = 1, use the transformation from voxels ogf physiens 2015
% study
% if transform = 2 use the one from AAL atlas
% IR 28/06/2017

if transform == 1
M = ...
    [ -3 0 0 81;...
    0 3 0 -115;...
    0 0 3 -53; ...
    0 0 0 1];
elseif transform == 2
M = ...
    [1,0,0,-91;...
    0,1,0,-126;...
    0,0,1,-72;...
    0,0,0,1];
end


%MNI-space to Voxel-space

 T=M(1:3,4);

 M=M(1:3,1:3);

 for i=1:3

MNI(i,:)=MNI(i,:)-T(i);

 end

 VOX=round(inv(M)*MNI);

return

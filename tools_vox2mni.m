function MNI = tools_vox2mni(VOX,transform)

% VOXel space to MNI space
% VOX has to be a column vector
% if transform = 1, use the transformation from voxels physiens study
% if transform = 2 use the one from AAL atlas
% IR 29/06/2017


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

 T=M(1:3,4);

 M=M(1:3,1:3);
 
 MNI=M*VOX;

 for i=1:3

MNI(i,:)=MNI(i,:)+T(i);

 end

return
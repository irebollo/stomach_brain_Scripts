function tools_writeMri(data,filename)
% data, filename
% Write 3d or 4d data in nifti file using fieldtrip
% parameters are specified for PHYSIENS experiment (dimensions)
% IR 26/03/2015

mri=[];

mri.transform = ...
    [ -3 0 0 81;...
    0 3 0 -115;...
    0 0 3 -53; ...
    0 0 0 1];
% physiens voxels transform

mri.coh = data;
mri.dim=[53,63,46];

cfg.parameter = 'coh';
cfg.filename =  filename;
cfg.filetype    = 'nifti';
cfg.coordsys      = 'spm';
cfg.scaling = 'no';
cfg.datatype = 'double';


ft_volumewrite(cfg,mri)

end
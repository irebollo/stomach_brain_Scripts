function MRI_group_preproc(subj_idx,smoothingSize)
% Toplevel batch for preprocessing with SPM8 by Anne-Dominique Lodeho
%
% If an input parameter 'confirm' is provided then questions will be
% asked interactively to perform or not each step of the preprocessings.
%
% To use this batch, please edit the beginning of this file and run it.
% The different parameters to modify are:
%  o rootdir:     Absolute path of the study 
%  o subjects:    Cell array of subjects name (that matches folders name)
%  o sessions:    Cell array of session folders to be used as template for
%                 all sessions. If a char array is given then all subsequent 
%                 folders are automatically extracted
%  o anatdir:     Relative path of the anatomical scan
%  o anatwc:      Wildcard for anatomical scan (regular expression)
%  o funcdirs:    Cell array of strings, each string being a
%                 relative path containing functional scans of a session.
%  o funcwc       Wildcard for functional scans (regular expression)
%  o TR:          Repetition time (for slice timing correction)
%  o ref_slice:   Reference slice (for slice timing correction)
%  o slice_order: Acquisition slice order: ['sequential_ascending'  | 'sequential_descending' |
%                                           'interleaved_ascending' | 'interleaved_descending' |
%                                           'interleaved_middletop']
%  o voxelsize:   Voxel size of the normalized data
%  o smoothing:   FWHM of kernel for smoothing functional scans
%  o logfile:     Filename of the text logfile
%
% At last, you need to specify which particular tasks to perform, through
% the TODO structure with binary fields (1|0):
%  o slice_timing: Slice Timing correction
%  o realign:      Movement correction
%  o normalize:    Combined segmentation and spatial normalisation
%  o coregister:   Anatomical and functional scan coregistration
%  o apply_norm:   Functional scans reslicing in template space
%  o smooth:       Functional scans smoothing
%  o run:          Execute batch job

%- Parameters for preprocessing a group of subjects
%----------------------------------------------------------------------
% rootdir = 'D:\Physiens Data\'; % Chemin du repertoire de la manip
rootdir =global_path2root;

% subjects = {'Subjects\subject15'}
subjects = cellstr(strcat('Subjects\','Subject', sprintf('%.2d',subj_idx)));

%sessions = {'fMRI/acquisition1/run1','fMRI/acquisition1/run2','fMRI/acquisition1/run3','fMRI/acquisition1/run4'};
% sessions = {'fMRI\acquisition1\RestingStateNoSTC'};
sessions = {'fMRI\acquisition1\RestingState'};


params.anatdir = 't1mri\acquisition1';
params.anatwc  = '^sPHYSIENS.*\.img$';%% les images anatomiques originales commencent par anat et sont en .nii c'était avant maintenant (au CHSA sur images GE anonymisées, ça commence par sAW)
%params.anatwc  = '^.*\.nii$'; 

[params.funcdirs{1:length(subjects)}] = deal(sessions);
params.funcwc  = '^fPHYSIENS.*\.img$';%% les images fonctionelles originales commencent par un f !
%params.funcwc  = '^.*\.nii$';

params.TR = 2;
params.ref_slice = 1;% otherwise the first slice can be used
%params.ref_slice = 40; % when slice_order=interleaved it is better to use the middle slice as reference. 
params.slice_order = 'sequential_ascending'; %'sequential_ascending', 'sequential_descending', 'interleaved_ascending', 'interleaved_descending', 'interleaved_middletop'
params.voxelsize = [3 3 3];
params.smoothing = smoothingSize;

params.logfile = 'preproc.log';

%- List of tasks to be performed
%----------------------------------------------------------------------
if nargin
    todo.slice_timing = 1;
    todo.realign      = 1;
    todo.normalize    = 1;
    todo.coregister   = 1;
    todo.apply_norm   = 1;
    todo.smooth       = 1;
    todo.run          = 1;
else
    todo.slice_timing = ask('yes','Apply slice timing');
    todo.realign      = ask('yes','Realign functional images');
    todo.normalize    = ask('yes','Compute normalization for anat');
    todo.coregister   = ask('yes','Coregister anat onto target functional image');
    todo.apply_norm   = ask('yes','Apply normalization to functional images');
    todo.smooth       = ask('yes','Smooth functional images');
    todo.run          = ask('yes','Run Batch');
end

%- Loop over subjects
%----------------------------------------------------------------------
p = params;
t0 = clock;
for n=1:length(subjects)
    t = clock;
    p.rootdir = fullfile(rootdir,subjects{n});
    p.funcdirs = params.funcdirs{n};
    MRI_single_preproc(p,todo);
    fprintf('Elapsed time for subject %s is %.2f seconds.\n',subjects{n},etime(clock,t));
end
fprintf('Overall elapsed time is %.2f seconds.\n',etime(clock,t0));

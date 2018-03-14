function global_createSubjectFolder(subj_idx)

rootDir= strcat(global_path2root,'\Subjects\');
% cd (rootDir)


% mkdir(strcat(rootDir,'Subject',sprintf('%.2d',subj_idx),'\Timeseries','\Heart'))

% mkdir(strcat(rootDir,'Subject',sprintf('%.2d',subj_idx),'\Timeseries','\ControlEGG_othesubject'))

mkdir(strcat(rootDir,'Subject',sprintf('%.2d',subj_idx),'\Timeseries','\Pupil'))


% cd 'C:\PHYSIENS\Physiens\SPM_Tools\mfiles'






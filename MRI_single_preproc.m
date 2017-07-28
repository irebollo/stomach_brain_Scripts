function MRI_single_preproc(parameters,todo)
% Create an SPM8 job for single subject's preprocessing by Anne-Dominique Lodeho
%
% PARAMETERS = structure with fields:
%   o rootdir:     Absolute path of the considered subject 
%   o anatdir:     Relative path of the anatomical scan
%   o funcdirs:    Cell array of strings, each string being a
%                  relative path containing functional scans of a session
%   o anatwc:      Wildcard for anatomical scan (regular expression)
%   o funcwc:      Wildcard for functional scans (regular expression)
%   o TR:          Repetition time (for slice timing correction)
%   o ref_slice:   Reference slice (for slice timing correction)
%   o slice_order: Acquisition slice order (for slice timing correction)
%   o voxelsize:   Size of voxels for writing functional scans (e.g.
%                  [2 2 2])
%   o smoothing:   Size of kernel for smoothing functional scans
%   o logfile:     Filename of the text logfile
%
% TODO = structure with binary fields (1 | 0):
%   o slice_timing: Slice Timing correction
%   o realign:      Movement correction
%   o normalize:    Combined segmentation and spatial normalisation
%   o coregister:   Anatomical and functional scan coregistration
%   o apply_norm:   Functional scans reslicing in template space
%   o smooth:       Functional scans smoothing
%   o run:          Execute batch job

%-Parameters checking
%----------------------------------------------------------------------
error(nargchk(2,2,nargin));

fields = fieldnames(parameters);
valid_fields = {'rootdir', 'anatdir', 'funcdirs', 'anatwc', 'funcwc', ...
    'ref_slice', 'slice_order', 'TR', 'voxelsize', 'smoothing', 'logfile'};
fielddiff = setdiff(fields,valid_fields);
if ~isempty(fielddiff)
    error('Unknown parameters: %s',sprintf('%s  ',fielddiff{:}));
end
fielddiff = setdiff(valid_fields,fields);
if ~isempty(fielddiff)
    error('Mandatory parameters: %s',sprintf('%s  ',fielddiff{:}));
end

fields = fieldnames(todo);
valid_fields = {'slice_timing', 'realign', 'normalize', ...
    'coregister', 'apply_norm', 'smooth', 'run'};
if ~isempty(setxor(fields,valid_fields))
    error('Invalid todo list.');
end

%- Logfile
%----------------------------------------------------------------------
rootdir = parameters.rootdir;

logfile = fullfile(rootdir,parameters.logfile);

logmsg(logfile,sprintf('Preprocessing data in folder %s...',rootdir));

%- Scanning for anatomical scan
%----------------------------------------------------------------------
logmsg(logfile,'Scanning for anatomical scan...');
spm_select('clearvfiles'); 
parameters.anatdir = spm_select('CPath',parameters.anatdir,rootdir);
anat = spm_select('List', parameters.anatdir, parameters.anatwc);
if isempty(anat)
    error('Cannot find anatomical scan "%s" in folder "%s"',parameters.anatwc,parameters.anatdir);
elseif size(anat,1) > 1
    error('Several files match anatomical scan "%s" in folder "%s"',parameters.anatwc,parameters.anatdir);
end
anat = fullfile(parameters.anatdir,deblank(anat(1,:)));
logmsg(logfile,sprintf('  found 1 file named %s.',anat));
manat = addprefixtofilenames(cellstr(anat),'m');

%- Scanning for functional scans
%----------------------------------------------------------------------
logmsg(logfile,'Scanning for functional scans...');
if ischar(parameters.funcdirs)
    logmsg(logfile,'Scanning for session folders...');
    [dummy, d] = spm_select('List',fullfile(rootdir,parameters.funcdirs),'.*');
    d(ismember(d,{'.' '..'}),:) = '';
    parameters.funcdirs = cellstr([repmat(spm_select('CPath','',parameters.funcdirs),size(d,1),1), d]);
    logmsg(logfile,sprintf('  found %d session folders',length(parameters.funcdirs)));
end
for n=1:length(parameters.funcdirs)
    parameters.funcdirs{n} = spm_select('CPath',parameters.funcdirs{n},rootdir);
    f = spm_select('List',parameters.funcdirs{n},parameters.funcwc);
    ff{n} = [repmat(spm_select('CPath','',parameters.funcdirs{n}),size(f,1),1), f];
end
aff = addprefixtofilenames(ff,'a'); % enlever 'a' auqnd on ne fait pas le slice timing
waff = addprefixtofilenames(aff,'w'); % put here the width of smoothing
logmsg(logfile,sprintf('  found %d files in %d session(s) ',sum(cellfun('size',ff,1)),length(ff)));
for n=1:length(parameters.funcdirs)
    logmsg(logfile,sprintf('    with %d files in session %s',size(ff{n},1),parameters.funcdirs{n}));
end

%- Jobs definition
%----------------------------------------------------------------------
jobs = {};
nbjobs = 0;

%- Slice Timing
%----------------------------------------------------------------------
if todo.slice_timing
    logmsg(logfile,sprintf('Slice timing on %d files starting with "%s"...',sum(cellfun('size',ff,1)),ff{1}(1,:)));
    nbjobs = nbjobs + 1;
    for n=1:length(parameters.funcdirs)
        jobs{nbjobs}.temporal{1}.st.scans{n} = cellstr(ff{n});
    end
    V = spm_vol(ff{1}(1,:));
    nbslices = V.dim(3);
    jobs{nbjobs}.temporal{1}.st.nslices = nbslices;
    jobs{nbjobs}.temporal{1}.st.tr = parameters.TR;
    jobs{nbjobs}.temporal{1}.st.ta = parameters.TR - parameters.TR / nbslices;
    switch parameters.slice_order
        case 'sequential_ascending'
            slice_order = 1:1:nbslices;
        case 'sequential_descending'
            slice_order = nbslices:-1:1;
        case 'interleaved_ascending'
            %slice_order = [1:2:nbslices 2:2:nbslices];
            if mod(nbslices,2) %odd slice number
                slice_order = [1:2:nbslices 2:2:nbslices];
                parameters.ref_slice = 1;
            else
                slice_order = [2:2:nbslices 1:2:nbslices];
                parameters.ref_slice = 2;
            end
        case 'interleaved_descending'
            slice_order = [nbslices:-2:1, nbslices-1:-2:1];
        case 'interleaved_middletop'
            for k = 1:nbslices
                slice_order(k) = round((nbslices-k)/2 + (rem((nbslices-k),2) * (nbslices - 1)/2)) + 1;
            end
        otherwise,
            error('Slice Timing specification failed.');
    end
    jobs{nbjobs}.temporal{1}.st.so = slice_order;
    jobs{nbjobs}.temporal{1}.st.refslice = parameters.ref_slice;
end

%- Realign
%----------------------------------------------------------------------
if todo.realign
    logmsg(logfile,sprintf('Realigning %d functional files onto "%s"...',sum(cellfun('size',aff,1)),aff{1}(1,:)));
    nbjobs = nbjobs + 1;
    for n=1:length(parameters.funcdirs)
        jobs{nbjobs}.spatial{1}.realign{1}.estwrite.data{n} = cellstr(aff{n});
    end
    jobs{nbjobs}.spatial{1}.realign{1}.estwrite.eoptions.quality = 1;
    jobs{nbjobs}.spatial{1}.realign{1}.estwrite.eoptions.interp = 5;
    jobs{nbjobs}.spatial{1}.realign{1}.estwrite.roptions.which = [0 1]; %- mean image only
%     jobs{nbjobs}.spatial{1}.realign{1}.estwrite.roptions.msk = 'D:\Physiens\scripts\files\Mask_NoEyeballs_Normalized_Voxels3mm.hdr';
end

%- Normalize
%----------------------------------------------------------------------
if todo.normalize
    logmsg(logfile,sprintf('Normalizing "%s" onto T1.nii',anat));
    nbjobs = nbjobs + 1;

    logmsg(logfile,sprintf('  Unified Segmentation on "%s"',anat));
    jobs{nbjobs}.spatial{1}.preproc.data = cellstr(anat);
    [pth,nm] = fileparts(anat);
    anat_matfile = fullfile(pth,[nm '_seg_sn.mat']);
    logmsg(logfile,sprintf('  Writing Normalized Bias Corrected "%s"',char(manat)));
    jobs{nbjobs}.spatial{2}.normalise{1}.write.subj.matname = cellstr(anat_matfile);
    jobs{nbjobs}.spatial{2}.normalise{1}.write.subj.resample = manat;
    jobs{nbjobs}.spatial{2}.normalise{1}.write.roptions.vox = [1 1 1];
end

%- Coregister
%----------------------------------------------------------------------
if todo.coregister
    logmsg(logfile,sprintf('Coregistering "%s" onto "%s"',aff{1}(1,:),anat));
    nbjobs = nbjobs + 1;
    jobs{nbjobs}.spatial{1}.coreg{1}.estimate.ref = manat;
    jobs{nbjobs}.spatial{1}.coreg{1}.estimate.source = cellstr(aff{1}(1,:));
    jobs{nbjobs}.spatial{1}.coreg{1}.estimate.other = cellstr(strvcat(aff));
    jobs{nbjobs}.spatial{1}.coreg{1}.estimate.other(1) = [];
end

%- Apply Normalize
%----------------------------------------------------------------------
if todo.apply_norm
    if ~todo.normalize
        logmsg(logfile,'Scanning for normalisation parameters...');
        anat_matfile = spm_select('List', parameters.anatdir, '^.*_sn.mat$');
        anat_matfile(~cellfun('isempty',regexp(cellstr(anat_matfile),'^.*inv_sn.mat$')),:) = [];
        if isempty(anat)
            error('Cannot find normalisation parameters "*_sn.mat" in folder "%s"',parameters.anatdir);
        elseif size(anat,1) > 1
            error('Several files match normalisation parameters "*_sn.mat" "%s" in folder "%s"',parameters.anatdir);
        end
        anat_matfile = fullfile(parameters.anatdir,deblank(anat_matfile(1,:)));
    end
    logmsg(logfile,sprintf('Apply normalization "%s" to %d files starting with "%s"...',anat_matfile,sum(cellfun('size',aff,1)),aff{1}(1,:)));
    nbjobs = nbjobs + 1;
    jobs{nbjobs}.spatial{1}.normalise{1}.write.subj(1).matname = cellstr(anat_matfile);
    jobs{nbjobs}.spatial{1}.normalise{1}.write.subj(1).resample = cellstr(strvcat(aff));
    jobs{nbjobs}.spatial{1}.normalise{1}.write.roptions.vox = parameters.voxelsize;
%     jobs{nbjobs}.spatial{1}.normalise{1}.write.subj(1).roptions.mask =  'D:\Physiens\scripts\files\Mask_NoEyeballs_Normalized_Voxels3mm.hdr';% Im adding an explicit mask of voxels inside brain
%     jobs{nbjobs}.spatial{1}.normalise{1}.write.prefix='w1';
    %jobs{nbjobs}.spatial{1}.normalise{1}.write.subj(1).roptions.wrap = [0 1 0];
    %jobs{nbjobs}.spatial{1}.normalise{1}.write.subj(1).roptions.preserve = 0;
    %jobs{nbjobs}.spatial{1}.normalise{1}.write.subj(1).roptions.bb = NaN;
end

%- Smooth
%----------------------------------------------------------------------
if todo.smooth
    logmsg(logfile,sprintf('Smoothing %d files ("%s"...) with fwhm = %d mm', ...
        sum(cellfun('size',waff,1)),waff{1}(1,:),parameters.smoothing));
    nbjobs = nbjobs + 1;
    jobs{nbjobs}.spatial{1}.smooth.data = cellstr(strvcat(waff));
    jobs{nbjobs}.spatial{1}.smooth.fwhm = [parameters.smoothing,parameters.smoothing,parameters.smoothing];
    jobs{nbjobs}.spatial{1}.smooth.prefix = strcat('s',num2str(parameters.smoothing));

end

%- Save and Run job
%----------------------------------------------------------------------
logmsg(logfile,sprintf('Job batch file saved in %s.',fullfile(rootdir,'jobs_preproc.mat')));
save(fullfile(rootdir,'jobs_preproc.mat'),'jobs');
if todo.run
    spm_jobman('run',jobs);
else
    spm_jobman('interactive',jobs);
    spm('show');
end

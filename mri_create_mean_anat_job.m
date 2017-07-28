%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4252 $)
%-----------------------------------------------------------------------
%%
matlabbatch{1}.spm.util.imcalc.input = {
                                        'D:\PHYSIENS\Subjects\Subject08\t1mri\acquisition1\wmsPHYSIENS_Sujet08-0009-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject09\t1mri\acquisition1\wmsPHYSIENS_Sujet09-0007-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject10\t1mri\acquisition1\wmsPHYSIENS_Sujet10-0009-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject13\t1mri\acquisition1\wmsPHYSIENS_Sujet13-0009-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject15\t1mri\acquisition1\wmsPHYSIENS_Sujet15-0009-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject16\t1mri\acquisition1\wmsPHYSIENS_Sujet16-0009-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject17\t1mri\acquisition1\wmsPHYSIENS_Sujet17-0007-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject18\t1mri\acquisition1\wmsPHYSIENS_Sujet18-0007-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject19\t1mri\acquisition1\wmsPHYSIENS_Sujet19-0009-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject20\t1mri\acquisition1\wmsPHYSIENS_Sujet20-0009-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject21\t1mri\acquisition1\wmsPHYSIENS_Sujet21-0009-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject22\t1mri\acquisition1\wmsPHYSIENS_Sujet22-0009-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject23\t1mri\acquisition1\wmsPHYSIENS_Sujet023-0007-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject24\t1mri\acquisition1\wmsPHYSIENS_Sujet24-0009-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject25\t1mri\acquisition1\wmsPHYSIENS_Sujet25-0007-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject26\t1mri\acquisition1\wmsPHYSIENS_Sujet26-0009-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject29\t1mri\acquisition1\wmsPHYSIENS_Sujet29-0009-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject31\t1mri\acquisition1\wmsPHYSIENS_Sujet31-0011-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject32\t1mri\acquisition1\wmsPHYSIENS_Sujet32-0009-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject33\t1mri\acquisition1\wmsPHYSIENS_Sujet33-0009-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject34\t1mri\acquisition1\wmsPHYSIENS_Sujet34-0009-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject35\t1mri\acquisition1\wmsPHYSIENS_sujet35-0007-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject36\t1mri\acquisition1\wmsPHYSIENS_Sujet36-0007-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject37\t1mri\acquisition1\wmsPHYSIENS_Sujet37-0009-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject38\t1mri\acquisition1\wmsPHYSIENS_Sujet38-0009-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject39\t1mri\acquisition1\wmsPHYSIENS_Sujet39-0007-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject40\t1mri\acquisition1\wmsPHYSIENS_Sujet40-0009-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject41\t1mri\acquisition1\wmsPHYSIENS_Sujet41-0009-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject43\t1mri\acquisition1\wmsPHYSIENS_Sujet43-0007-00001-000224-01.img'
                                        'D:\PHYSIENS\Subjects\Subject44\t1mri\acquisition1\wmsPHYSIENS_Sujet44-0009-00001-000224-01.img'
                                        };
%%
matlabbatch{1}.spm.util.imcalc.output = 'meanAnat.img';
matlabbatch{1}.spm.util.imcalc.outdir = {'D:\Physiens\Subjects\SecondLevel\'};
matlabbatch{1}.spm.util.imcalc.expression = '(sum(X))/30';
matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

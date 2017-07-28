 function outpoints = mni2tal(inpoints)
 
 % MNI2TAL - MNI to Talairach coordinates (best guess)
 %
 % outpoints = mni2tal(inpoints)
 %
 % inpoints - 3xN matrix, ie inpoints = [X Y Z]'
 %
 % outpoints - the Talairach coordinates (3xN matrix)
 %
 % See also, TAL2MNI, MNI2TAL_MATRIX &
 % http://www.mrc-cbu.cam.ac.uk/Imaging/Common/mnispace.shtml
 %
 
 % $Revision: 1.3 $ $Date: 2004/12/10 18:15:20 $
 
 % Licence:  GNU GPL, no express or implied warranties
 % Matthew Brett 10/8/99, matthew.brett@mrc-cbu.cam.ac.uk
 % modified 02/2003, Darren.Weber_at_radiology.ucsf.edu
 %                   - removed dependence on spm_matrix and
 %                     abstracted the matrix in case it changes
 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 [s1,s2] = size(inpoints);
 if s1 ~= 3,
     error('input must be a 3xN matrix')
 end
 
 % Transformation matrices, different zooms above/below AC
 M2T = mni2tal_matrix;
 
 inpoints = [inpoints; ones(1, size(inpoints, 2))];
 
 tmp = inpoints(3,:) < 0;  % 1 if below AC
 
 inpoints(:,  tmp) = (M2T.rotn * M2T.downZ) * inpoints(:,  tmp);
 inpoints(:, ~tmp) = (M2T.rotn * M2T.upZ  ) * inpoints(:, ~tmp);
 
 outpoints = inpoints(1:3, :);
 
return
 function outpoints = tal2mni(inpoints)
 
 % TAL2MNI - Talairach to MNI coordinates
 %
 % outpoints = tal2mni(inpoints)
 %
 % inpoints  - 3xN matrix of coordinates
 %
 % outpoints - the coordinate matrix with MNI points
 %
 % See also, MNI2TAL & the best guess discussion at
 % http://www.mrc-cbu.cam.ac.uk/Imaging/Common/mnispace.shtml
 %
 
 % $Revision: 1.3 $ $Date: 2005/01/21 06:22:14 $
 
 % Licence:  GNU GPL, no express or implied warranties
 % Matthew Brett 2/2/01, matthew.brett@mrc-cbu.cam.ac.uk
 % modified 02/2003, Darren.Weber_at_radiology.ucsf.edu
 %                   - swapped inv() for slash equivalent
 %                   - removed dependence on spm_matrix
 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 [s1,s2] = size(inpoints);
 if s1 ~= 3,
     error('input must be a 3xN matrix')
 end
 
 % Transformation matrices, different zooms above/below AC
 M2T = mni2tal_matrix;
 
 inpoints = [inpoints; ones(1, size(inpoints, 2))];
 
 tmp = inpoints(3,:) < 0;  % 1 if below AC
 
 inpoints(:,  tmp) = (M2T.rotn * M2T.downZ) \ inpoints(:,  tmp);
 inpoints(:, ~tmp) = (M2T.rotn * M2T.upZ  ) \ inpoints(:, ~tmp);
 
 outpoints = inpoints(1:3, :);
 
return



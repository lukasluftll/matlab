function a = ht2affine3d(ht)
% HT2AFFINE3D Convert 4x4 homogeneous transformation to affine3d object.
% 
%   See also AFFINE3D. 

% Copyright 2016 Alexander Schaefer

a = affine3d(ht');

end
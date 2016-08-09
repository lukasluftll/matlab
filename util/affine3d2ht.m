function ht = affine3d2ht(a)
% AFFINE3D2HT Convert affine3d object to 4x4 homogeneous transformation.
% 
%   See also AFFINE3D. 

% Copyright 2016 Alexander Schaefer

ht = a.T';

end

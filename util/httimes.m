function c = httimes(a, b)
% HTTIMES Homogeneous transformation multiply.
%   C = HTTIMES(A, B) multiplies N homogeneous transformation matrices by
%   N homogeneous transformation matrices or N transformation matrices
%   by N homogeneous coordinate vectors.
%
%   A is a 4x4xN page-wise concatenation of homogeneous transformation
%   matrices.
%
%   Matrix-by-matrix multiplication
%   -------------------------------
%   B is a 4x4xN page-wise concatenation of homogeneous transformation
%   matrices.
%
%   C is a 4x4xN page-wise concatenation of homogeneous transformation
%   matrices, where the n-th page contains the product of 
%   A(:,:,n) * B(:,:,n).
%
%   Matrix-by-vector multiplication
%   -------------------------------
%   B is a 4xN matrix. The rows of B contain homogenous coordinate vectors.
%
%   C is a 4xN matrix. The n-th row of C contains the product of 
%   A(:,:,n) * B(:,n).
%
%   Example:
%      httimes([ones(3), rand(3,1); 0,0,0,1], [rand(3,1); 1])
%
%   See also TIMES, MTIMES, ISHRT.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check number of input arguments.
narginchk(2, 2)

% Check contents of A.
if ~isnumeric(a)
    error('A must contain numeric values only.')
end

% Check size of A.
if ndims(a) > 3 || size(a, 1) ~= 4 || size(a, 2) ~= 4
    error('A must be a 4x4xN matrix.')
end

% Check contents of B.
if ~isnumeric(b)
    error('B must contain numeric values only.')
end

% Check size of B.
if ndims(b) > 3 || size(b, 1) ~= 4 || (ndims(b)==3 && size(b, 2)~=4)
    error('B must be a 4x4xN matrix or a 4xN matrix.')
end

%% Perform multiplication.
if ndims(b) == 3   % Matrix-matrix multiply.
    c = zeros(size(a));
    for i = 1 : size(a, 3)
        c(:,:,i) = a(:,:,i) * b(:,:,i);
    end
else   % Matrix-vector multiply.
    c = zeros(size(b));
    for i = 1 : size(b, 2)
        c(:,i) = a(:,:,i) * b(:,i);
    end
    c = reshape(c, 4, []);
end

end

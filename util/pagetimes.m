function c = pagetimes(a, b)
% PAGETIMES Page-wise matrix-matrix multiply.
%   C = PAGETIMES(A, B) multiplies K matrices by K matrices or K matrices
%   by K vectors.
%
%   A is a IxJxK page-wise concatenation of IxJ matrices.
%
%   Matrix-by-matrix multiplication
%   -------------------------------
%   B is a IxJxK page-wise concatenation IxJ matrices.
%
%   C is a IxJxK page-wise concatenation of matrices, where the k-th page 
%   contains the product of A(:,:,k) * B(:,:,k).
%
%   Matrix-by-vector multiplication
%   -------------------------------
%   B is a JxK matrix. The columns of B contain J-element vectors.
%
%   C is a IxK matrix. The k-th column of C contains the product of 
%   A(:,:,k) * B(:,k).
%
%   Example:
%      pagetimes([eye(3), rand(3,1); 0,0,0,1], rand(4,1))
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

% Check dimensionality of A.
if ndims(a) > 3
    error('A must be a 3D matrix.')
end

% Check contents of B.
if ~isnumeric(b)
    error('B must contain numeric values only.')
end

% Check dimensionality of B.
if ndims(b) > 3 
    error('B must be a 2D or a 3D matrix.')
end

% Make sure the sizes of A and B match.
msg = 'Sizes of A and B do not match.';
if ismatrix(b) 
    if size(a,2)~=size(b,1) || size(a,3) ~= size(b,2)
        error(msg)
    end
elseif ~all(size(a)==size(b))
    error(msg)
end

%% Perform multiplication.
if ismatrix(b)   % Matrix-vector multiply.
    c = zeros(size(a,1), size(b,2));
    for i = 1 : size(b, 2)
        c(:,i) = a(:,:,i) * b(:,i);
    end
else   % Matrix-matrix multiply.
    c = zeros(size(a));
    for i = 1 : size(a, 3)
        c(:,:,i) = a(:,:,i) * b(:,:,i);
    end
end

end

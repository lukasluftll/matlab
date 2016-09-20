function c = pagetimes(a, b)
% PAGETIMES Page-wise matrix multiply.
%   C = PAGETIMES(A,B) multiplies the pages of two matrices A and B.
%
%   A is a IxJxN matrix. 
%   B is a JxKxN matrix. 
%   C is a IxKxN matrix whose n-th page contains A(:,:,n)*B(:,:,n).
%
%   Example:
%      pagetimes(magic(3,3,3), rand(3,3,3))
%
%   See also TIMES, MTIMES.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check number of input arguments.
narginchk(2, 2)

% Check the size of the input matrices.
if size(a,3)~=size(b,3) || size(a,2)~=size(b,1)
    error('Sizes of A and B do not match.')
end

%% Perform page-wise multiplication.
c = zeros(size(a,1), size(b,2), size(a,3));
parfor i = 1 : size(a, 3)
    c(:,:,i) = a(:,:,i) * b(:,:,i);
end

end

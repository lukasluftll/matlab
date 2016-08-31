function x = getol(a, b, tol)
% GETOL Greater than or equal with tolerance.
%   X = GETOL(A, B) does element by element comparisons between A and B and
%   returns a logical matrix X of the same size. X(i) is true if 
%   A(i) > B(i) or if A(i) is within tolerance of B(i). 
% 
%   Two values A(i) and B(i) are within tolerance if:
%      abs(A(i)-B(i)) <= TOL * max(abs([A(:);B(:)]))
%
%   The tolerance is 1e-6 for single-precision inputs and 1e-12 for 
%   double-precision inputs.
%
%   X = GETOL(A, B, TOL) specifies an absolute scalar tolerance TOL.
%
%   Example:
%      a = [1 11 111];
%      b = a + eps(a);
%      getol(a, b)
%
%   See also LETOL, ISMEMBERTOL, EPS.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check number of input arguments.
narginchk(2, 3)

% Check the type of the inputs.
if ~isnumeric(a) || ~isnumeric(b) || (nargin==3 && ~isnumeric(tol))
    error('All input arguments must be numeric values.')
end

% If the tolerance is not given, define it.
if nargin < 3
    if isa(a, 'single') || isa(b, 'single')
        tol = 1e-6;
    else
        tol = 1e-12;
    end
end    

%% Compare input.
% Do not use MATLAB's built-in function ISMEMBERTOL because it leads to
% memory leaks!
x = a > b | abs(a-b) <= tol*max(abs([a(:);b(:)]));

end

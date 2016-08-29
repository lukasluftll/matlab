function x = getol(a, b, tol)
% GETOL Less than or equal with tolerance.
%   X = GETOL(A, B) does element by element comparisons between A and B and
%   returns a logical matrix X of the same size. X(i) is true if 
%   A(i) > B(i) or if A(i) is within tolerance of B(i). The tolerance is
%   1e-6 for single-precision inputs and 1e-12 for double-precision inputs.
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

switch nargin
    case 2
        x = a > b | ismembertol(a, b);
    case 3
        x = a > b | ismembertol(a, b, tol);
    otherwise
        narginchk(2, 3)
end

end

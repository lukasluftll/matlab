function sizechk(varargin)
% SIZECHK Check equal matrix sizes.
%   SIZECHK(X1,X2,X3,...) checks whether all input matrices X1, X2, X3,
%   etc. have the same size.
%
%   Example:
%      a = magic(4);
%      b = magic(5);
%      sizechk(a, b)
%
%   See also SIZE, GVCHK.

% Copyright 2016 Alexander Schaefer

%% Validate input and output.
nargoutchk(0, 0)
narginchk(0, +Inf)

%% Check sizes pairwise.
% Get the name of the input argument.
if nargin > 0
    ref = varargin{1};
    refname = upper(inputname(1));
end

% Compare the size of the first matrix to all other matrices.
for i = 2 : nargin
    % Get the name of the other matrix.
    comp = varargin{i};
    compname = upper(inputname(i));
    
    % Generate a string that identifies the mismatching matrices.
    inputnames = blanks(0);
    if ~isempty(refname) || ~isempty(compname)
        inputnames = ['of ', refname, ' and ', compname, ' '];
    end
    
    % Check the number of dimensions of both matrices.
    if ndims(ref) ~= ndims(comp)
        error(['Numbers of dimensions ', inputnames, 'do not match.'])
    end
    
    % Check the size of both matrices.
    if any(size(varargin{1}) ~= size(varargin{i}))
        error(['Matrix sizes ', inputnames, 'do not match.'])
    end
end

end

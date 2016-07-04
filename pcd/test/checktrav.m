function checktrav(i, iExpected, t, tExpected)
% CHECKTRAV Assert function for results provided by TRAV.
%
%   See also TRAV.

% Check number of input arguments.
narginchk(1, 4);

% Define the exceptions.
iSizeExcpt = MException('trav:iSize', ...
    'Incorrect index matrix size.');
tSizeExcpt = MException('trav:tSize', ...
    'Incorrect line parameter vector size.');
iIncorrectExcpt = MException('trav:iIncorrect', ...
    'Returned indices incorrect.');
tIncorrectExcpt = MException('trav:tIncorrect', ...
    'Returned line parameters incorrect.');
iZeroExcpt = MException('trav:iZero', ...
    'Returned indices contain zero values.');

% Check results.
assert(all(all(i)), iZeroExcpt.identifier, iZeroExcpt.message);

if nargin > 1
    assert(all(size(i) == size(iExpected)), iSizeExcpt.identifier, ...
        iSizeExcpt.message);
    assert(all(size(t) == size(tExpected)), tSizeExcpt.identifier, ...
        tSizeExcpt.message);
    assert(all(all(i == iExpected | isnan(i) & isnan(iExpected))), ...
        iIncorrectExcpt.identifier, iIncorrectExcpt.message);
    assert(all(abs(t-tExpected) < 1e-9 | (isnan(t) & isnan(tExpected))),...
        tIncorrectExcpt.identifier, tIncorrectExcpt.message);
end

end

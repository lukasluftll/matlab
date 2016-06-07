function checktrav(i, iExpected, t, tExpected)
% CHECKTRAV Assert function for results provided by TRAV.
%
%   See also TRAV.

% Define the exceptions.
iIncorrectExcpt = MException('trav:iIncorrect', ...
    'Returned indices incorrect.');
tIncorrectExcpt = MException('trav:tIncorrect', ...
    'Returned line parameters incorrect.');

% Check results.
assert(all(all(i - iExpected < eps)), iIncorrectExcpt.identifier, ...
    iIncorrectExcpt.message);
assert(all(t - tExpected < eps | (isnan(t) & isnan(tExpected))), ...
    tIncorrectExcpt.identifier, tIncorrectExcpt.message);

end

function f = datefile(ext)
% DATEFILE Generate file name based on date and time.
%   DATEFILE returns a string that represents the current date and time in
%   the format '2016-11-30-19-07-59'.
%
%   DATEFILE(EXT) appends the given extension to the date-time string.
%   The extension may or may not include a leading dot.
%
%   Example:
%      datefile('.txt')
%
%   See also DATESTR.

% Copyright 2016 Alexander Schaefer

%% Validate input.
if nargin < 1
    ext = blanks(0);
else
    validateattributes(ext, {'char'}, {'row', 'nonempty'}, '', 'EXT')
    ext = [repmat('.', 1, ext(1)~='.'), ext];
end

%% Generate string.
f = [datestr(now, 'yyyy-mm-dd-HH-MM-SS'), ext];

end

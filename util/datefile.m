function f = datefile(ext)

if nargin < 1
    ext = blanks(0);
else
    validateattributes(ext, {'char'}, {'row', 'nonempty'}, '', 'EXT')
    ext = [repmat('.', 1, ext(1)~='.'), ext];
end

f = [datestr(now, 'yyyy-mm-dd-HH-MM-SS'), ext];

end

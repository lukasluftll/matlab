function sizechk(varargin)

nargoutchk(0, 1)
narginchk(0, +Inf)

if nargin > 0
    ref = varargin{1};
    refname = inputname(1);
end

for i = 2 : nargin
    comp = varargin{i};
    compname = inputname(i);
    
    inputnames = blanks(0);
    if isempty(refname) || isempty(compname)
        inputnames = ['of ',refname,' and ',compname,' '];
    end
       
    if ndims(ref) ~= ndims(comp)
        error(['Numbers of dimensions ', inputnames, 'do not match.'])
    end
    
    if any(size(varargin{1}) ~= size(varargin{i}))
        error(['Sizes ', inputnames, 'do not match.'])
    end
end

end

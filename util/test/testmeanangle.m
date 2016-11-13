% Unit test for MEANANGLE.

% Copyright 2016 Alexander Schaefer

% SHARED VARIABLES SECTION.

%% Empty set of angles.
theta = [];

m = meanangle(theta);
assert(isempty(m));

m = meanangle(theta, 'vectorsum');
assert(isempty(m));

%% Equal angles.
theta = ones(5,30);
msize = [1,30];
mval = 1;

m = meanangle(theta);
assert(all(size(m) == msize));
assert(all(m == mval));

m = meanangle(theta, 'vectorsum');
assert(all(size(m) == msize));
assert(all(ismembertol(m, mval)));

%% Symmetrically distributed angles.
theta = [-0.3 : 0.2 : +0.3; -0.6 : 0.4 : +0.6];
dim = 2;
msize = [2,1];
mval = 0;

m = meanangle(theta, dim);
assert(all(size(m) == msize));
assert(all(m == mval));

m = meanangle(theta, 2, 'vectorsum');
assert(all(size(m) == msize));
assert(all(abs(m-mval) < eps));

%% Randomly distributed angles.
theta = rand(1,100);

m = meanangle(theta);
assert(0 < m && m < 1);

m = meanangle(theta, 'vectorsum');
assert(0 < m && m < 1);

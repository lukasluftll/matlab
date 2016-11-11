% Unit test for MEANANGLE.

% Copyright 2016 Alexander Schaefer

% SHARED VARIABLES SECTION.

%% Empty set of angles.
m = meanangle([]);
assert(isempty(m));

m = meanangle([], 'vectorsum');
assert(isempty(m));

%% Equal angles applied.
m = meanangle(ones(2,30));
assert(m == 1);

m = meanangle(ones(2,30), 'vectorsum');
assert(ismembertol(m, 1));

%% Symmetrically distributed angles.
m = meanangle(-1 : 0.01 : +1);
assert(m == 0);

m = meanangle(-0.3 : 0.2 : +0.3);
assert(m == 0);

m = meanangle(-1 : 0.01 : +1, 'vectorsum');
assert(m < eps);

m = meanangle(-0.3 : 0.2 : +0.3, 'vectorsum');
assert(m < eps);

%% Randomly distributed angles.
m = meanangle(rand(1, 100));
assert(0 < m && m < 1);

m = meanangle(rand(1, 100), 'vectorsum');
assert(0 < m && m < 1);

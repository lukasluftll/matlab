% Unit test for MEANANGLE.

% Copyright 2016 Alexander Schaefer

% SHARED VARIABLES SECTION.

%% Empty set of angles.
m = meanangle([]);
assert(isempty(m));

%% Equal angles.
m = meanangle(ones(2,30));
assert(m == 1);

%% Symmetrically distributed angles.
m = meanangle(-1 : 0.01 : +1);
assert(m == 0 || m == 2*pi);

m = meanangle(-0.3 : 0.2 : +0.3);
assert(m == 0 || m == 2*pi);

%% Randomly distributed angles.
m = meanangle(rand(1, 100));
assert(0 < m && m < 1);

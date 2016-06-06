% Unit test for TRAV.

% Copyright 2016 Alexander Schaefer

% SHARED VARIABLES SECTION.

%% Ray in x direction
origin = [0, 0, 0];
ray = [1, 0, 0];
vol = [-5, -5, -5, 5, 5, 5];
res = 1;
i = trav(origin, ray, vol, res);

result = [5, 6, 6;
    4, 6, 6;
    3, 6, 6;
    2, 6, 6;
    1, 6, 6];
assert(all(all(i == result)));

%% Ray in y direction

%% Ray in z direction

%% Ray in x-y direction

%% Ray in x-z direction

%% Ray in y-z direction

%% Ray outside grid

%% Ray along edge of grid

%% Ray through corner of grid

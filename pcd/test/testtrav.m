% Unit test for TRAV.

% Copyright 2016 Alexander Schaefer

% SHARED VARIABLES SECTION.

%% Origin ray in x direction
origin = [0, 0, 0];
ray = [1, 0, 0];
vol = [-5, -5, -5, 5, 5, 5];
res = 1;
[i, t] = trav(origin, ray, vol, res);

iExpected = [6, 6, 6;
    7, 6, 6;
    8, 6, 6;
    9, 6, 6;
    10, 6, 6];
tExpected = [NaN; 1; 2; 3; 4; 5];
checktrav(i, iExpected, t, tExpected);

%% Origin ray in y direction
origin = [0, 0, 0];
ray = [0, 1, 0];
vol = [-3, -4, -5, 3, 4, 5];
res = 1;
[i, t] = trav(origin, ray, vol, res);

iExpected = [4, 5, 6;
    4, 6, 6;
    4, 7, 6;
    4, 8, 6];
tExpected = [NaN; 1; 2; 3; 4];
checktrav(i, iExpected, t, tExpected);

%% Origin ray in z direction
origin = [0, 0, 0];
ray = [0, 0, 1];
vol = [-1, -2, -3, 1, 2, 3];
res = 1;
[i, t] = trav(origin, ray, vol, res);

iExpected = [2, 3, 4;
    2, 3, 5,
    2, 3, 6];
tExpected = [NaN; 1; 2; 3];

%% Origin ray in x-y direction
origin = [0, 0, 0];
ray = [1, 1, 0];
vol = [-5, -5, -5, 5, 5, 5];
res = 1;
[i, t] = trav(origin, ray, vol, res);

iExpected = [6, 6, 6;
    7, 7, 6;
    8, 8, 6;
    9, 9, 6;
    10, 10, 6];
tExpected = [NaN; 1; 2; 3; 4; 5] * sqrt(2);
checktrav(i, iExpected, t, tExpected);

%% Ray in x-z direction
origin = [0, 0, 0];
ray = [1, 0, 1];
vol = [-5, -5, -5, 5, 5, 5];
res = 1;
[i, t] = trav(origin, ray, vol, res);

iExpected = [6, 6, 6;
    7, 6, 7;
    8, 6, 8;
    9, 6, 9;
    10, 6, 10];
tExpected = [NaN; 1; 2; 3; 4; 5] * sqrt(2);
checktrav(i, iExpected, t, tExpected);

%% Ray in y-z direction
origin = [0, 0, 0];
ray = [0, 1, 1];
vol = [-5, -5, -5, 5, 5, 5];
res = 1;
[i, t] = trav(origin, ray, vol, res);

iExpected = [6, 6, 6;
    6, 7, 7;
    6, 8, 8;
    6, 9, 9;
    6, 10, 10];
tExpected = [NaN; 1; 2; 3; 4; 5] * sqrt(2);
checktrav(i, iExpected, t, tExpected);

%% Ray in x-y-z direction
origin = [0, 0, 0];
ray = [1, 1, 1];
vol = [10, 10, 10, 15, 15, 18];
res = 1;
[i, t] = trav(origin, ray, vol, res);

iExpected = [1, 1, 1;
    2, 2, 2;
    3, 3, 3;
    4, 4, 4;
    5, 5, 5;
    6, 6, 6];
tExpected = [10; 11; 12; 13; 14; 15; 16];
checktrav(i, iExpected, t, tExpected);

%% Ray outside grid
origin = [1.5, 1.5, 1.5];
ray = [1, 0, 0];
vol = [2, 1, 1, 4, 5, 5];
res = 1;
[i, t] = trav(origin, ray, vol, res);

iExpected = [1, 1, 1;
    2, 1, 1;
    3, 1, 1];
tExpected = [0.5; 1.5; 2.5; 3.5];
checktrav(i, iExpected, t, tExpected);

%% Ray along edge of grid
origin = [1, 1, 1];
ray = [0, 0, -1];
vol = [-1, -1, -1, 1, 1, 1];
res = 1;
[i, t] = trav(origin, ray, vol, res);

assert(all(~size(i)));
assert(all(~size(t)));

%% Ray through corner of grid
origin = [-1, 1, 0];
ray = [1, -1, 0];
vol = [0, 0, 0, 10, 10, 10];
res = 1;
[i, t] = trav(origin, ray, vol, res);

iExpected = [1, 1, 1];
tExpected = [1; 1] * sqrt(2);
checktrav(i, iExpected, t, tExpected);

%% Resolution > 1
origin = [0, 47, 48];
ray = [1, 0, 0];
vol = [45, 45, 45; 55, 55, 55];
res = 5;
[i, t] = trav(origin, ray, vol, res);

iExpected = [1, 1, 1;
    2, 1, 1;
    3, 1, 1];
tExpected = [9; 10; 11];
checktrav(i, iExpected, t, tExpected);

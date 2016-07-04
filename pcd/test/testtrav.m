% Unit test for TRAV.

% Copyright 2016 Alexander Schaefer

% SHARED VARIABLES SECTION.

%% Origin ray in x direction
origin = [0, 0, 0];
ray = [10, 0, 0];
xgv = -5 : 5; 
ygv = -5 : 5; 
zgv = -5 : 5;
[i, t] = trav(origin, ray, xgv, ygv, zgv);

iExpected = [6, 6, 6;
    7, 6, 6;
    8, 6, 6;
    9, 6, 6;
    10, 6, 6];
tExpected = [0; 0.1; 0.2; 0.3; 0.4; 0.5];
checktrav(i, iExpected, t, tExpected);

%% Origin ray in y direction
origin = [0, 0, 0];
ray = [0, 10, 0];
xgv = -3 : 3;
ygv = -4 : 4;
zgv = -5 : 5;
[i, t] = trav(origin, ray, xgv, ygv, zgv);

iExpected = [4, 5, 6;
    4, 6, 6;
    4, 7, 6;
    4, 8, 6];
tExpected = [0; 0.1; 0.2; 0.3; 0.4];
checktrav(i, iExpected, t, tExpected);

%% Origin ray in z direction
origin = [0, 0, 0];
ray = [0, 0, 10];
xgv = -1 : 1;
ygv = -2 : 2;
zgv = -3 : 3;
[i, t] = trav(origin, ray, xgv, ygv, zgv);

iExpected = [2, 3, 4;
    2, 3, 5;
    2, 3, 6];
tExpected = [0; 0.1; 0.2; 0.3];
checktrav(i, iExpected, t, tExpected);

%% Origin ray in x-y direction
origin = [0, 0, 0];
ray = [10, 10, 0];
xgv = -5 : 5;
ygv = -5 : 5;
zgv = -5 : 5;
[i, t] = trav(origin, ray, xgv, ygv, zgv);

iExpected = [6, 6, 6;
    7, 7, 6;
    8, 8, 6;
    9, 9, 6;
    10, 10, 6];
tExpected = [0; 0.1; 0.2; 0.3; 0.4; 0.5];
checktrav(i, iExpected, t, tExpected);

%% Ray in x-z direction
origin = [0, 0, 0];
ray = [10, 0, 10];
xgv = -5 : 5;
ygv = -5 : 5;
zgv = -5 : 5;
[i, t] = trav(origin, ray, xgv, ygv, zgv);

iExpected = [6, 6, 6;
    7, 6, 7;
    8, 6, 8;
    9, 6, 9;
    10, 6, 10];
tExpected = [0; 0.1; 0.2; 0.3; 0.4; 0.5];
checktrav(i, iExpected, t, tExpected);

%% Ray in y-z direction
origin = [0, 0, 0];
ray = [0, 10, 10];
xgv = -5 : 5;
ygv = -5 : 5;
zgv = -5 : 5;
[i, t] = trav(origin, ray, xgv, ygv, zgv);

iExpected = [6, 6, 6;
    6, 7, 7;
    6, 8, 8;
    6, 9, 9;
    6, 10, 10];
tExpected = [0; 0.1; 0.2; 0.3; 0.4; 0.5];
checktrav(i, iExpected, t, tExpected);

%% Ray in x-y-z direction
origin = [0, 0, 0];
ray = [100, 100, 100];
xgv = 10 : 15;
ygv = 10 : 15;
zgv = 10 : 18;
[i, t] = trav(origin, ray, xgv, ygv, zgv);

iExpected = [1, 1, 1;
    2, 2, 2;
    3, 3, 3;
    4, 4, 4;
    5, 5, 5];
tExpected = [0.1; 0.11; 0.12; 0.13; 0.14; 0.15];
checktrav(i, iExpected, t, tExpected);

%% Origin ray ending in grid
origin = [0, 0, 0];
ray = [3, 3, 3];
xgv = -5 : 5;
ygv = -5 : 5;
zgv = -5 : 5;
[i, t] = trav(origin, ray, xgv, ygv, zgv);

iExpected = [6, 6, 6;
    7, 7, 7;
    8, 8, 8;
    9, 9, 9];
tExpected = [0; 1/3; 2/3; 1; 1];
checktrav(i, iExpected, t, tExpected);

origin = [0, 0, 0];
ray = [3.5, 3.5, 3.5];
xgv = -5 : 5;
ygv = -5 : 5;
zgv = -5 : 5;
[i, t] = trav(origin, ray, xgv, ygv, zgv);

iExpected = [6, 6, 6;
    7, 7, 7;
    8, 8, 8;
    9, 9, 9];
tExpected = [0; 1/3.5; 2/3.5; 3/3.5; 1];
checktrav(i, iExpected, t, tExpected);

%% Ray entering grid from lower limits
origin = [1.5, 1.5, 1.5];
ray = [10, 0, 0];
xgv = 2 : 4;
ygv = 1 : 5;
zgv = 1 : 5;
[i, t] = trav(origin, ray, xgv, ygv, zgv);

iExpected = [1, 1, 1;
    2, 1, 1];
tExpected = [0.05; 0.15; 0.25];
checktrav(i, iExpected, t, tExpected);

%% Ray entering grid from upper limits
origin = [5.3, 0.3, 0.3];
ray = [-10, 0, 0];
xgv = -1 : 5;
ygv = -5 : 5;
zgv = -5 : 5;
[i, t] = trav(origin, ray, xgv, ygv, zgv);

iExpected = [6, 6, 6;
    5, 6, 6;
    4, 6, 6;
    3, 6, 6;
    2, 6, 6;
    1, 6, 6];
tExpected = [0.03; 0.13; 0.23; 0.33; 0.43; 0.53; 0.63];
checktrav(i, iExpected, t, tExpected);

%% Ray entering grid from upper x-limit
origin = [100, 47, 48];
ray = [-100, 0, 0];
xgv = 45 : 5 : 55;
ygv = 45 : 5 : 55;
zgv = 45 : 5 : 55;
[i, t] = trav(origin, ray, xgv, ygv, zgv);

iExpected = [2, 1, 1;
    1, 1, 1];
tExpected = [0.45; 0.50; 0.55];
checktrav(i, iExpected, t, tExpected);

%% Ray along edge of grid
origin = [1, 1, 1];
ray = [0, 0, -10];
xgv = -1 : 1;
ygv = -1 : 1;
zgv = -1 : 1;
[i, t] = trav(origin, ray, xgv, ygv, zgv);

assert(all(~size(i)));
assert(all(~size(t)));

%% Ray through corner inside grid
origin = [-1, 1, 0];
ray = [10, -10, 0];
xgv = 0 : 10;
ygv = 0 : 10;
zgv = 0 : 10;
[i, t] = trav(origin, ray, xgv, ygv, zgv);

iExpected = [1, 1, 1];
tExpected = [0.1; 0.1];
checktrav(i, iExpected, t, tExpected);

%% Ray through corner outside grid
origin = [-1, 1, 0];
ray = [10, -10, 0];
xgv = -5 : 0;
ygv = -5 : 0;
zgv = -5 : 0;
[i, t] = trav(origin, ray, xgv, ygv, zgv);

assert(all(~size(i)));
assert(all(~size(t)));

%% Ray length not equal 1
origin = [-0.9, -0.5, -0.5];
ray = [0.01, 0, 0];
xgv = -1 : 1;
ygv = -1 : 1;
zgv = -1 : 1;
[i, t] = trav(origin, ray, xgv, ygv, zgv);

iExpected = [1, 1, 1];
tExpected = [0; 1];
checktrav(i, iExpected, t, tExpected);

%% Zero-index error
origin = [0, 0, 0];
ray = [-8887, -3785, 2588];
xgv = -8.5380 : 0.1 : 22.1352;
ygv = -3.6631 : 0.1 : 18.7318;
zgv = -5.9373 : 0.1 : 5.0340;
[i, t] = trav(origin, ray, xgv, ygv, zgv);

checktrav(i);

%% Infinite loop error
origin = [0, 0, 0];
ray = [9962, 0, -872];
xgv = -8.5380 : 0.1 : 22.1352;
ygv = -3.6631 : 0.1 : 18.7318;
zgv = -5.9373 : 0.1 : 5.0340;
[i, t] = trav(origin, ray, xgv, ygv, zgv);

checktrav(i);

%% gpuArray input arguments
origin = gpuArray([0, 0, 0]);
ray = gpuArray([3, 3, 3]);
xgv = gpuArray(-5 : 5);
ygv = gpuArray(-5 : 5);
zgv = gpuArray(-5 : 5);
[i, t] = trav(origin, ray, xgv, ygv, zgv);

iExpected = NaN(numel(xgv) + numel(ygv) + numel(zgv), 3);
iExpected(1:4,:) = [6, 6, 6;
    7, 7, 7;
    8, 8, 8;
    9, 9, 9];
tExpected = NaN(numel(xgv) + numel(ygv) + numel(zgv), 1);
tExpected(1:5) = [0; 1/3; 2/3; 1; 1];
checktrav(i, iExpected, t, tExpected);

%% Mixed input arguments
origin = [0, 0, 0];
ray = [3.5, 3.5, 3.5];
xgv = -5 : 5;
ygv = -5 : 5;
zgv = gpuArray(-5 : 5);
[i, t] = trav(origin, ray, xgv, ygv, zgv);

iExpected = NaN(numel(xgv) + numel(ygv) + numel(zgv), 3);
iExpected(1:4,:) = [6, 6, 6;
    7, 7, 7;
    8, 8, 8;
    9, 9, 9];
tExpected = NaN(numel(xgv) + numel(ygv) + numel(zgv), 1);
tExpected(1:5) = [0; 1/3.5; 2/3.5; 3/3.5; 1];
checktrav(i, iExpected, t, tExpected);

origin = gpuArray([1.5, 1.5, 1.5]);
ray = [10, 0, 0];
xgv = 2 : 4;
ygv = 1 : 5;
zgv = 1 : 5;
[i, t] = trav(origin, ray, xgv, ygv, zgv);

iExpected = NaN(numel(xgv) + numel(ygv) + numel(zgv), 3);
iExpected(1:2,:) = [1, 1, 1;
    2, 1, 1];
tExpected = NaN(numel(xgv) + numel(ygv) + numel(zgv), 1);
tExpected(1:3) = [0.05; 0.15; 0.25];
checktrav(i, iExpected, t, tExpected);

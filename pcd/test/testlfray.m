% Unit test for LFRAY.

% Copyright 2016 Alexander Schaefer, Lukas Luft

% SHARED VARIABLES SECTION.
prior = 0.5;
map = gridmap([[0.4,0.7,0.8,0.1]; ones(4)], 1, [0,0,0], prior);
stdsp = trvec2tform([0,0,0]);
stdrlim = [1,3];
az = 0;
ele = 0;
pnr =.4;
normalizer = (1-pnr)/(.7+.8);

%% NaN return
r = NaN;
ls = laserscan(sdtsp, az, ele, r, stdrlim);
[p,l] = lfray(ls, map,pnr);

pExpected = pnr;
lExpected = 0;
assert(eqtol(p, pExpected));
assert(eqtol(l, lExpected));

%% Real return
r = 2.2;
ls = laserscan(stdsp, az, ele, r, stdrlim);
[p,l] = lfray(ls, map,pnr);

pExpected = 1;
lExpected = normalizer*.8;
assert(eqtol(p, pExpected));
assert(eqtol(l, lExpected));

%% Radius range exceeds extent of map
% No return.
rlim =  [2,10];
r = NaN;
ls = laserscan(stdsp, az, ele, r, rlim);
[p,l] = lfray(ls, map,pnr);

pExpected = pnr;
lExpected = 0;
assert(eqtol(p, pExpected));
assert(eqtol(l, lExpected));

% Ray returns.
rlim =  [2,10];
normalizer = (1-pnr)/(.8+.1+6*prior);
r = 7.5;
ls = laserscan(stdsp, az, ele, r, rlim);
[p,l] = lfray(ls, map,pnr);

pExpected = 1;
lExpected = log(normalizer*prior);

assert(eqtol(p, pExpected));
assert(eqtol(l, lExpected));

%% Ray from outside into map
sp = trvec2tform([-2,0,0]);
r = 2.5;
ls = laserscan(sp, az, ele, r, stdrlim);
[p,l] = lfray(ls, map,pnr);
normalizer= normalizer = (1-pnr)/(prior+0.4);

pExpected = 1;
lExpected = log(normalizer*.4);
assert(eqtol(p, pExpected));
assert(eqtol(l, lExpected));

%% Ray from outside through map
sp = trvec2tform([1.5,-1,0]);
rlim = [0,20];
r = 15;
ls = laserscan(sp, as, ele, r, rlim);
[p,l] = lfray(ls, map,pnr);
normalizer=(1-pnr)/(20*prior);

pExpected = 1;
lExpected = log(normalizer*prior);
assert(eqtol(p, pExpected));
assert(eqtol(l, lExpected));

%% Ray does not touch map
% No return.
r = NaN;
ls = laserscan(stdsp, az, ele, r, stdrlim);
[p,l] = lfray(ls, map,pnr);

pExpected = pnr;
lExpected = 0;
assert(eqtol(p, pExpected));
assert(eqtol(l, lExpected));

% Ray returns.
sp = trvec2tfrom([100,100,100]);
r = 2.5;
ls = laserscan(sp, az, ele, r, stdrlim);
[p,l] = lfray(ls, map,pnr);
normalizer=(1-pnr)/(2*prior);


pExpected = 1;
lExpected = log(normalizer*prior);
assert(eqtol(p, pExpected));
assert(eqtol(l, lExpected));

%% Zero ray
r = 0;
rlim = [0,10];
normalizer = (1-pnr)/(0.4+0.7+.8+.1+6*prior);
ls = laserscan(stdsp, as, ele, r, rlim);
[p,l] = lfray(ls, map, pnr);

pExpected = 1;
lExpected = log(.4*normalizer);
assert(eqtol(p, pExpected));
assert(eqtol(l, lExpected));

%% Ray through map corner.
sp = treul2tform([-1,2,0], [-pi/4,0,0]);
r = 2;
ls = laserscan(sp, as, ele, r, stdrlim);
normalizer=(1-pnr)/(prior*sqrt(2)+(2-sqrt(2))*.4);

[p,l] = lfray(ls, map, pnr);

pExpected = 1;
lExpected = log(.4*normalizer);
assert(eqtol(p, pExpected));
assert(eqtol(l, lExpected));

%% Ray along upper edge of map
sp = trvec2tform([0,0,1]);
r = NaN;
ls = laserscan(sp, az, ele, r, stdrlim);
[p,l] = lfray(ls, map, pnr);

pExpected = pnr;
lExpected = 0;
assert(eqtol(p, pExpected));
assert(eqtol(l, lExpected));

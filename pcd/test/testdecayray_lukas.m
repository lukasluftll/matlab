% Test for decayray

% define grid, prior and map
gridx = 0:4;
gridy = 0:4;
gridz = 0:1;
prior = .5;
data = [.4,.7,.8,.1;1,1,1,1;1,1,1,1;1,1,1,1]';
map =   voxelmap(data,gridx,gridx,gridz,prior);
RLIM =  [1,3];
% define Laserbeam properties
SP =        eye(4);
AZIMUTH =   0;
ELEVATION = 0;
RADIUS =    3;


%% 
RADIUS =NaN;
ls =     laserscan(SP, AZIMUTH, ELEVATION, RADIUS, RLIM);
[p l] = decayray(ls,map);
pExpected = 1-exp(-.4)+exp(-.4-.7-.8);
lExpected = 0;
assert(p==pExpected);
assert(l==lExpected);

%% 
RADIUS =2.2;
ls =     laserscan(SP, AZIMUTH, ELEVATION, RADIUS, RLIM);
[p l] = decayray(ls,map);
pExpected = 1;
lExpected = -.4-.7-.8*.2+log(.8);
assert(lExpected-eps<l||l<lExpected+eps);
assert(pExpected==p);

%% 
RLIM =  [2,10];
RADIUS =NaN;
ls =     laserscan(SP, AZIMUTH, ELEVATION, RADIUS, RLIM);
[p l] = decayray(ls,map)
pExpected = 1-exp(-.4-.7)+exp(-.4-.7-.8-.5-5*.5)
lExpected = 0;
assert(pExpected-eps<l||l<pExpected+eps);
assert(l==lExpected+eps);

%% 
RLIM =  [2,10];
RADIUS =7.5;
ls =     laserscan(SP, AZIMUTH, ELEVATION, RADIUS, RLIM);
[p l] = decayray(ls,map)
pExpected = -.4-.7-.8-.1-5.5*.5
lExpected = 0;
assert(pExpected-eps<l||l<pExpected+eps);
assert(l==lExpected+eps);
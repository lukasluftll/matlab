% Test for decayray

% define grid, prior and map
gridx = 0:4;
gridy = 0:4;
gridz = 0:1;
prior = 0;
data = .4*ones(5)
data = [.4,.7,.8,.1;1,1,1,1;1,1,1,1;1,1,1,1]';
map =   voxelmap(data,gridx,gridx,gridz,prior);
RLIM =  [1,3];
% define Laserbeam properties
SP =        eye(4);
AZIMUTH =   0;
ELEVATION = 0;
RADIUS =    3;


%% 
% define Laserbeam
RADIUS =.2
ls =     laserscan(SP, AZIMUTH, ELEVATION, RADIUS, RLIM);
% calculate ray probabilities
[p l] = refray(ls,map)
pExpected = log((1-.4)*.3*.8)
%
%assert(l==lExpected);
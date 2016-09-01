% Run all lidar map experiments.

% Upload parameters.
lidarparams

% Run all models for all datasets.
dataset = 'campus'; %#ok<*NASGU>
buildpcmap
model = 'decay';
buildlidarmap, evallidarmap
model = 'ref';
buildlidarmap, evallidarmap
model = 'lf';
buildlidarmap, evallidarmap

dataset = 'schauinsland';
buildpcmap
model = 'decay';
buildlidarmap, evallidarmap
model = 'ref';
buildlidarmap, evallidarmap
model = 'lf';
buildlidarmap, evallidarmap

dataset = 'mooswald';
buildpcmap
model = 'decay';
buildlidarmap, evallidarmap
model = 'ref';
buildlidarmap, evallidarmap
model = 'lf';
buildlidarmap, evallidarmap

% Run all lidar map experiments.

% Upload parameters.
lidarparams

% Run all models for all datasets.
dataset = 'campus'; %#ok<*NASGU>
buildpcmap
model = 'decay';
buildlidarmap, kllidarmap%, klilidarmap
model = 'ref';
buildlidarmap, kllidarmap%, klilidarmap
model = 'lf';
buildlidarmap, kllidarmap%, klilidarmap

dataset = 'schauinsland';
buildpcmap
model = 'decay';
buildlidarmap, kllidarmap%, klilidarmap
model = 'ref';
buildlidarmap, kllidarmap%, klilidarmap
model = 'lf';
buildlidarmap, kllidarmap%, klilidarmap

dataset = 'mooswald';
buildpcmap
model = 'decay';
buildlidarmap, kllidarmap%, klilidarmap
model = 'ref';
buildlidarmap, kllidarmap%, klilidarmap
model = 'lf';
buildlidarmap, kllidarmap%, klilidarmap

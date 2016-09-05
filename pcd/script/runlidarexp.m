% Run all lidar map experiments.

% Run all models for all datasets.
dataset = 'campus'; %#ok<*NASGU>
buildpcmap
model = 'decay';
buildlidarmap, kllidarmap, klilidarmap, problidarmap
model = 'ref';
buildlidarmap, kllidarmap, klilidarmap, problidarmap
model = 'lf';
buildlidarmap, kllidarmap, klilidarmap, problidarmap

dataset = 'schauinsland';
buildpcmap
model = 'decay';
buildlidarmap, kllidarmap, klilidarmap, problidarmap
model = 'ref';
buildlidarmap, kllidarmap, klilidarmap, problidarmap
model = 'lf';
buildlidarmap, kllidarmap, klilidarmap, problidarmap

dataset = 'mooswald';
buildpcmap
model = 'decay';
buildlidarmap, kllidarmap, klilidarmap, problidarmap
model = 'ref';
buildlidarmap, kllidarmap, klilidarmap, problidarmap
model = 'lf';
buildlidarmap, kllidarmap, klilidarmap, problidarmap

% Run all lidar map experiments.

% Declare global variables.
global dataset model

% Upload parameters.
lidarparams

% Define datasets to process.
datasets = {'campus', 'schauinsland', 'mooswald'};

% Process all datasets using all models.
for i = 1 : numel(datasets)
    dataset = datasets{i};
    model = 'decay'; %#ok<*NASGU>
    buildpcmap, buildlidarmap, evallidarmap
    model = 'ref';
    buildlidarmap, evallidarmap
    model = 'lf';
    buildlidarmap, evallidarmap
end

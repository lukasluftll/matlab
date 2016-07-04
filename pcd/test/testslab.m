% Unit test for SLAB.

% Copyright 2016 Alexander Schaefer

% SHARED VARIABLES SECTION.
% Define exceptions.
hitIncorrectExcpt = MException('slab:hitIncorrect', ...
    'Intersection detection incorrect.');
hitFalsePositiveExcpt = MException('slab:hitFalsePositive', ...
    'Intersection detected although ray and box do not intersect.');
hitFalseNegativeExcpt = MException('slab:hitFalseNegative', ...
    'No intersection detected although ray and box intersect.');
paramsNotNanExcpt = MException('slab:paramsNotNan', ...
    'Line parameters for intersection are not NaN.');
paramsNanExcpt = MException('slab:paramsNan', ...
    'Line parameters for intersection are not set to real, but to NaN.');
paramsIncorrectExcpt = MException('slab:paramsIncorrect', ...
    'Line parameters for intersection are incorrect.');

%% Ray through corner inside box
support = [-3, 3, 0];
ray = [1, -1, 0];
box = [0, 0, 0, 10, 10, 10];
[hit, t] = slab(support, ray, box);

assert(hit, hitFalseNegativeExcpt.identifier, ...
    hitFalseNegativeExcpt.message);
assert(all(t == [3, 3]), paramsIncorrectExcpt.identifier, ...
    paramsIncorrectExcpt.message);

%% Ray through corner outside box
support = [1, -1, -1];
ray = [1, 1, 1];
box = [-1, -1, -1, 1, 1, 1];
[hit, t] = slab(support, ray, box);

assert(~hit, hitFalsePositiveExcpt.identifier, ...
    hitFalsePositiveExcpt.message);
assert(all(isnan(t)), paramsNotNanExcpt.identifier, ...
    paramsNotNanExcpt.message);

%% Ray through edge inside box
support = [-22, -44, 5];
ray = [0, 0, -1];
box = [-22, -44, -60, 11, 5, -50];
[hit, t] = slab(support, ray, box);

assert(hit, hitFalseNegativeExcpt.identifier, ...
    hitFalseNegativeExcpt.message);
assert(all(t == [55, 65]), paramsIncorrectExcpt.identifier, ...
    paramsIncorrectExcpt.message);

%% Ray through edge outside box
support = [2, 0, 1];
ray = [-1, 0, 0];
box = [0, 0, 0, 1, 1, 1];
[hit, t] = slab(support, ray, box);

assert(~hit, hitFalsePositiveExcpt.identifier, ...
    hitFalsePositiveExcpt.message);
assert(all(isnan(t)), paramsIncorrectExcpt.identifier, ...
    paramsNotNanExcpt.message);

%% Ray on surface inside box
support = [-5, -5, -1];
ray = [1, 1, 0];
box = [-1, -1, -1, 1, 1, 1];
[hit, t] = slab(support, ray, box);

assert(hit, hitFalseNegativeExcpt.identifier, ...
    hitFalseNegativeExcpt.message);
assert(all(t == [4, 6]), paramsIncorrectExcpt.identifier, ...
    paramsIncorrectExcpt.message);

%% Ray on surface outside box
support = [5, 2, 2];
ray = [-1, 0, 0];
box = [-2, -2, -2, 2, 2, 2];
[hit, t] = slab(support, ray, box);

assert(~hit, hitFalseNegativeExcpt.identifier, ...
    hitFalseNegativeExcpt.message);
assert(all(isnan(t)), paramsNotNanExcpt.identifier, ...
    paramsNotNanExcpt.message);

%% Intersection of ray with box
support = [-3, -3, -3];
ray = [0.4, 0.2, 1.1];
box = [-3, -3, -3, 3, 3, 3];
[hit, t] = slab(support, ray, box);

assert(hit, hitFalseNegativeExcpt.identifier, ...
    hitFalseNegativeExcpt.message);
assert(~any(isnan(t)), paramsIncorrectExcpt.identifier, ...
    paramsIncorrectExcpt.message);

%% No intersection of ray with box
support = [23, 59, 1];
ray = [3.42, 4.61, 0.09];
box = [-1, -1, -1, 1, 1, 1];
[hit, t] = slab(support, ray, box);

assert(~hit, hitFalsePositiveExcpt.identifier, ...
    hitFalsePositiveExcpt.message);
assert(all(isnan(t)), paramsNotNanExcpt.identifier, ...
    paramsNotNanExcpt.message);

%% Ray support inside box.
support = [0, 0, 0];
ray = [1, 1, 1];
box = [-5, -5, -5, 5, 12, 13];
[hit, t] = slab(support, ray, box);

assert(hit, hitFalseNegativeExcpt.identifier, ...
    hitFalseNegativeExcpt.message);
assert(all(t == [-5, 5]), paramsIncorrectExcpt.identifier, ...
    paramsIncorrectExcpt.message);

%% Multiple rays and boxes
support = [-3, 3, 0; ...
    1, -1, -1; ...
    -22, -44, 5; ...
    2, 0, 1; ...
    -5, -5, -1; ...
    5, 2, 2; ...
    23, 59, 1; ...
    0, 0, 0];
ray = [1, -1, 0; ...
    1, 1, 1; ...
    0, 0, -1; ...
    -1, 0, 0; ...
    1, 1, 0; ...
    -1, 0, 0; ...
    3.42, 4.61, 0.09; ...
    1, 1, 1];
box = [0, 0, 0, 10, 10, 10; ...
    -1, -1, -1, 1, 1, 1; ...
    -22, -44, -60, 11, 5, -50; ...
    0, 0, 0, 1, 1, 1; ...
    -1, -1, -1, 1, 1, 1; ...
    -2, -2, -2, 2, 2, 2; ...
    -1, -1, -1, 1, 1, 1; ...
    -5, -5, -5, 5, 12, 13];
[hit, t] = slab(support, ray, box);

assert(...
    all(hit == logical([1; 0; 1; 0; 1; 0; 0; 1])), ...
    hitIncorrectExcpt.identifier, ...
    hitIncorrectExcpt.message);
assert(all(all(isequaln(t, [3, 3; ...
    NaN, NaN; ...
    55, 65; ...
    NaN, NaN; ...
    4, 6; ...
    NaN, NaN; ...
    NaN, NaN; ...
    -5, 5]))), ...
    paramsIncorrectExcpt.identifier, paramsIncorrectExcpt.message);

%% gpuArray input arguments
support = [-3, 3, 0; ...
    1, -1, -1; ...
    -22, -44, 5; ...
    2, 0, 1; ...
    -5, -5, -1; ...
    5, 2, 2; ...
    23, 59, 1; ...
    0, 0, 0];
ray = [1, -1, 0; ...
    1, 1, 1; ...
    0, 0, -1; ...
    -1, 0, 0; ...
    1, 1, 0; ...
    -1, 0, 0; ...
    3.42, 4.61, 0.09; ...
    1, 1, 1];
box = [0, 0, 0, 10, 10, 10; ...
    -1, -1, -1, 1, 1, 1; ...
    -22, -44, -60, 11, 5, -50; ...
    0, 0, 0, 1, 1, 1; ...
    -1, -1, -1, 1, 1, 1; ...
    -2, -2, -2, 2, 2, 2; ...
    -1, -1, -1, 1, 1, 1; ...
    -5, -5, -5, 5, 12, 13];
[hit, t] = slab(gpuArray(support), gpuArray(ray), gpuArray(box));

assert(...
    all(hit == logical([1; 0; 1; 0; 1; 0; 0; 1])), ...
    hitIncorrectExcpt.identifier, ...
    hitIncorrectExcpt.message);
assert(all(all(isequaln(t, [3, 3; ...
    NaN, NaN; ...
    55, 65; ...
    NaN, NaN; ...
    4, 6; ...
    NaN, NaN; ...
    NaN, NaN; ...
    -5, 5]))), ...
    paramsIncorrectExcpt.identifier, paramsIncorrectExcpt.message);

%% Mixed-type input arguments
support = [-3, 3, 0; ...
    1, -1, -1; ...
    -22, -44, 5; ...
    2, 0, 1; ...
    -5, -5, -1; ...
    5, 2, 2; ...
    23, 59, 1; ...
    0, 0, 0];
ray = [1, -1, 0; ...
    1, 1, 1; ...
    0, 0, -1; ...
    -1, 0, 0; ...
    1, 1, 0; ...
    -1, 0, 0; ...
    3.42, 4.61, 0.09; ...
    1, 1, 1];
box = [0, 0, 0, 10, 10, 10; ...
    -1, -1, -1, 1, 1, 1; ...
    -22, -44, -60, 11, 5, -50; ...
    0, 0, 0, 1, 1, 1; ...
    -1, -1, -1, 1, 1, 1; ...
    -2, -2, -2, 2, 2, 2; ...
    -1, -1, -1, 1, 1, 1; ...
    -5, -5, -5, 5, 12, 13];
[hit, t] = slab(gpuArray(support), ray, box);

assert(...
    all(hit == logical([1; 0; 1; 0; 1; 0; 0; 1])), ...
    hitIncorrectExcpt.identifier, ...
    hitIncorrectExcpt.message);
assert(all(all(isequaln(t, [3, 3; ...
    NaN, NaN; ...
    55, 65; ...
    NaN, NaN; ...
    4, 6; ...
    NaN, NaN; ...
    NaN, NaN; ...
    -5, 5]))), ...
    paramsIncorrectExcpt.identifier, paramsIncorrectExcpt.message);

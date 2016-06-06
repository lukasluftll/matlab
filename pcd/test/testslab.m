% Unit test for SLAB.

% Copyright 2016 Alexander Schaefer

% SHARED VARIABLES SECTION.
% Define exceptions.
hitFalsePositiveExcpt = MException('slab:hitFalsePositive', ...
    'Intersection detected although ray and box do not intersect.');
hitFalseNegativeExcpt = MException('slab:hitFalseNegative', ...
    'No intersection detected although ray and box intersect.');
paramsNotNanExcpt = MException('slab:paramsNotNan', ...
    'Line parameters for intersection are not NaN.');
paramsIncorrectExcpt = MException('slab:paramsIncorrect', ...
    'Line parameters for intersection are incorrect.');

%% Ray through corner inside box
support = [-3, 3, 0];
ray = [1, -1, 0];
box = [0, 0, 0, 10, 10, 10];
[hit, t] = slab(support, ray, box);

assert(hit == true, hitFalseNegativeExcpt.identifier, ...
    hitFalseNegativeExcpt.message);
assert(all(t == [3; 3]), paramsIncorrectExcpt.identifier, ...
    paramsIncorrectExcpt.message);

%% Ray through corner outside box
support = [1, -1, -1];
ray = [1, 1, 1];
box = [-1, -1, -1, 1, 1, 1];
[hit, t] = slab(support, ray, box);

assert(hit == false, hitFalsePositiveExcpt.identifier, ...
    hitFalsePositiveExcpt.message);
assert(all(isnan(t)), paramsNotNanExcpt.identifier, ...
    paramsNotNanExcpt.message);

%% Ray through edge inside box
support = [-22, -44, 5];
ray = [0, 0, -1];
box = [-22, -44, -60, 11, 5, -50];
[hit, t] = slab(support, ray, box);

assert(hit == true, hitFalseNegativeExcpt.identifier, ...
    hitFalseNegativeExcpt.message);
assert(all(t == [55; 65]), paramsIncorrectExcpt.identifier, ...
    paramsIncorrectExcpt.message);

%% Ray through edge outside box
support = [2, 0, 1];
ray = [-1, 0, 0];
box = [0, 0, 0, 1, 1, 1];
[hit, t] = slab(support, ray, box);

assert(hit == false, hitFalsePositiveExcpt.identifier, ...
    hitFalsePositiveExcpt.message);
assert(all(isnan(t)), paramsIncorrectExcpt.identifier, ...
    paramsNotNanExcpt.message);

%% Ray on surface inside box
support = [-5, -5, -1];
ray = [1, 1, 0];
box = [-1, -1, -1, 1, 1, 1];
[hit, t] = slab(support, ray, box);

assert(hit == true, hitFalseNegativeExcpt.identifier, ...
    hitFalseNegativeExcpt.message);
assert(all(t == [4; 6]), paramsIncorrectExcpt.identifier, ...
    paramsIncorrectExcpt.message);

%% Ray on surface outside box
support = [5, 2, 2];
ray = [-1, 0, 0];
box = [-2, -2, -2, 2, 2, 2];
[hit, t] = slab(support, ray, box);

assert(hit == false, hitFalseNegativeExcpt.identifier, ...
    hitFalseNegativeExcpt.message);
assert(all(isnan(t)), paramsNotNanExcpt.identifier, ...
    paramsNotNanExcpt.message);

% Intersection

%% Multiple rays and boxes
support = [2, 0, 1; 2, 0, 1];
ray = [-1, 0, 0; -1, 0, 0];
box = [0, 0, 0, 1, 1, 1; 0, 0, 0, 1, 1, 1];
[hit, t] = slab(support, ray, box);

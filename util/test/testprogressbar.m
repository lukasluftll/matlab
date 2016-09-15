% Unit test for PROGRESSBAR.

% Copyright 2016 Alexander Schaefer

% SHARED VARIABLES SECTION.
n = 1000;
nWarn = 300;

%% Single-threaded execution
disp('Single-threaded execution')
progressbar(n)
for i = 1 : n
    progressbar
end

%% Single-threaded execution with warnings
disp('Single-threaded execution with warnings')
progressbar(n)
for i = 1 : n
    progressbar
    
    if rem(i,nWarn) == 0
        warning('warning text')
    end
end

%% Execution in parfor
disp('Multi-threaded execution')
progressbar(n)
parfor i = 1 : n
    progressbar
end

%% Execution in parfor with warnings
disp('Multi-threaded execution with warnings')
progressbar(n)
parfor i = 1 : n
    progressbar
    
    if rem(i,nWarn) == 0
        warning('warning text')
    end
end

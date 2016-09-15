% Unit test for PROGRESSBAR.

% Copyright 2016 Alexander Schaefer

% SHARED VARIABLES SECTION.

%% Single-threaded execution
disp('Single-threaded execution')
n = 1000;
progressbar(n)
for i = 1 : n
    progressbar
end

%% Single-threaded execution with warnings
disp('Single-threaded execution with warnings')
n = 1000;
progressbar(n)
for i = 1 : n
    progressbar
    
    if rem(i,300) == 0
        warning('warning text')
    end
end

%% Execution in parfor
disp('Multi-threaded execution')
n = 1000;
progressbar(n)
parfor i = 1 : n
    progressbar
end

%% Execution in parfor with warnings
disp('Multi-threaded execution with warnings')
n = 1000;
progressbar(n)
parfor i = 1 : n
    progressbar
    
    if rem(i,300) == 0
        warning('warning text')
    end
end

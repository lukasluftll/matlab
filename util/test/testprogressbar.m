% Unit test for PROGRESSBAR.

% Copyright 2016 Alexander Schaefer

% SHARED VARIABLES SECTION.
% Duration of single test case in [s].
t = 3;

%% Single-threaded execution
disp('Single-threaded execution')
n = 1000;
parprogress(n)
for i = 1 : n
    pause(t/n)
    parprogress
end

%% Single-threaded execution with warnings
disp('Single-threaded execution with warnings')
n = 1000;
parprogress(n)
for i = 1 : n
    pause(t/n)
    parprogress
    
    if rem(i/300) == 0
        warning('warning text')
    end
end

%% Execution in parfor
disp('Single-threaded execution')
n = 1000;
parprogress(n)
parfor i = 1 : n
    pause(t/n)
    parprogress
end

%% Execution in parfor with warnings
disp('Single-threaded execution with warnings')
n = 1000;
parprogress(n)
parfor i = 1 : n
    pause(t/n)
    parprogress
    
    if rem(i/300) == 0
        warning('warning text')
    end
end

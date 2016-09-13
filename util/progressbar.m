function progressbar(n)
% PROGRESSBAR Progress bar in command window.
%   PROGRESSBAR(N) displays an empty progress bar in the command window.
%   N is a positive integer that denotes the number of computations to 
%   come.
%
%   PROGRESSBAR increments the progress bar. After N calls to PROGRESSBAR,
%   the bar shows full completion.
%
%   PROGRESSBAR(0) must be called to close the progress bar.
%
%   PROGRESSBAR works both in a single-threaded and in a multi-threaded
%   setting.
%
%   Example:
%      n = 100;
%      progressbar(n);
%      parfor i = 1 : n
%         pause(0.01); 
%         progressbar;
%      end
%      progressbar(0);
%
%   See also WAITBAR.

% Copyright 2016 Alexander Schaefer
% PROGRESSBAR is inspired by TEXTPROGRESSBAR by Paul Proteus
% (https://www.mathworks.com/matlabcentral/fileexchange/28067-
% text-progress-bar) and by progressbar by Jeremy Scheff
% (https://de.mathworks.com/matlabcentral/fileexchange/32101-
% progress-monitor--progress-bar--that-works-with-parfor).

%% Validate input.
% Check number of input arguments.
narginchk(0, 1)

% Check type and content of input argument, if any.
if nargin > 0
    validateattributes(n, {'numeric'}, ...
        {'integer', 'nonnegative', 'numel',1}, '', 'N')
end

%% Define parameters.
% Width of the progress bar.
barwidth = 50;

% Name of the file in which the progress is stored.
filename = '.progressbar.txt';

%% Create, advance, or close progress bar.
if nargin > 0
    if n > 0
        pgbopen(n);
    else
        pgbclose;
    end
else
    pgbadvance;
end

%% Helper functions.
    function pgbopen(n)
        % PGBOPEN Open progress bar.
        
        % Store number of computations to come in file.
        pgbsave(n);
        
        % Display computation progress.
        pgbprint(0);
    end

    function pgbadvance
        % PGBADVANCE Advance progress bar by one step.
        
        % Check if the hidden file exists.
        if ~exist(filename, 'file')
            error('Auxiliary file not found. Run PROGRESSBAR(N) first.')
        end
        
        % Append a line to the auxiliary file.
        fid = fopen(filename, 'a');
        fprintf(fid, '1\n');
        fclose(fid);
        
        % Count the lines in the auxiliary file to get the progress.
        fid = fopen(filename, 'r');
        progress = fscanf(fid, '%d');
        fclose(fid);
        progress = numel(progress(2:end)) / progress(1);
        
        % If the file contains invalid numbers, throw error.
        if ~isfinite(progress)
            error('Auxiliary file corrupted. Run PROGRESSBAR(N).')
        end

        % Display progress.
        if progress < 1
            pgbprint(progress)
        else
            pgbclose;
        end
    end

    function pgbclose
        % PGBCLOSE Close progress bar.
        
        % Reset the value in the file.
        pgbsave(0);
        
        % Display progress.
        pgbprint(1);
    end

    function pgbprint(progress)
        % PGBPRINT Display progress.
        
        % Create percentage string.
        percentStr = sprintf('%3.0f%%', round(progress*100));
        
        % Display progress bar.
        display([repmat(char(8), 1, (barwidth+9) * (progress>0)), ...
            char(10), percentStr, ...
            '[', repmat('=', 1, round(progress*barwidth)), '>', ...
            repmat(' ', 1, barwidth-round(progress*barwidth)), ']'])
    end

    function pgbsave(n)
        % PGBSAVE Save number of computations to come to file.
        
        fid = fopen(filename, 'w');
        if fid < 0
            error('Failed to create auxiliary file in folder %s.', pwd)
        end
        fprintf(fid, '%d\n', n);
        fclose(fid);
    end
end

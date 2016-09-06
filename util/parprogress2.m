function parprogress2(n)
% PARPROGRESS Progress bar in command line.
%   parprogress works by creating a file called parprogress.txt in
%   your working directory, and then keeping track of the parfor loop's
%   progress within that file. This workaround is necessary because parfor
%   workers cannot communicate with one another so there is no simple way
%   to know which iterations have finished and which haven't.
%
%   parprogress(N) initializes the progress monitor for a set of N
%   upcoming calculations.
%
%   parprogress updates the progress inside your parfor loop and
%   displays an updated progress bar.
%
%   parprogress(0) deletes parprogress.txt and finalizes progress
%   bar.
%
%   To suppress output from any of these functions, just ask for a return
%   variable from the function calls, like PERCENT = parprogress which
%   returns the percentage of completion.
%
%   Example:
%
%      N = 100;
%      parprogress(N);
%      parfor i=1:N
%         pause(rand); % Replace with real code
%         parprogress;
%      end
%      parprogress(0);
%
%   See also PARFOR.

% By Jeremy Scheff - jdscheff@gmail.com - http://www.jeremyscheff.com/

%% Validate input.
% Check number of input arguments.
narginchk(0, 1)

% Check type and content of input argument, if any.
if nargin > 0
    validateattributes(n, {'numeric'}, ...
        {'integer', 'nonnegative', 'numel',1}, '', 'N')
end

%%
width = 50;
filename = '.parprogress.txt';

if nargin > 0
    if n > 0
        pgbopen(n);
    else
        pgbclose;
    end
else
    pgbadvance;
end

%%
    function pgbopen(n)
        % PGBINIT Initialize progress bar.
        
        % Store number of calculations to come in hidden file.
        fid = fopen(filename, 'w');
        if fid < 0
            error('Failed to create auxiliary file in folder %s.', pwd)
        end
        fprintf(fid, '%d\n', n);
        fclose(fid);
        
        % Display computation progress.
        disp(['  0%[>', repmat(' ', 1, width), ']']);
    end

    function pgbclose
        % PGBCLOSE Close progress bar.
        
        % Delete the hidden file.
        pgbopen(0);
        
        % Display full progress bar.
        disp([repmat(char(8), 1, (width+9)), char(10), ...
            '100%[', repmat('=', 1, width+1), ']']);
    end

    function pgbadvance
        % PGBUPDATE Update progress bar.
        
        % Check if the hidden file exists.
        if ~exist(filename, 'file')
            error('Auxiliary file not found. Run PARPROGRESS(N) first.');
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
        
        if ~isfinite(progress)
            error('Run PARPROGRESS(N) first.')
        end
        
        if progress >= 1
            pgbclose;
        else
            percent = sprintf('%3.0f%%', round(progress*100));
            disp([repmat(char(8), 1, (width+9)), char(10), ...
                percent, '[', repmat('=', 1, round(progress*width)), ...
                '>', repmat(' ', 1, width - round(progress*width)), ']']);
        end
    end
end

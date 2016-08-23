function progressbar(x)
% PROGRESSBAR Display progress bar in command window.
%   PROGRESSBAR(X) displays a progress bar in the command window.
%   
%   PROGRESSBAR(X) with X being a string prints a message and displays a
%   progress bar showing zero percent task completion.
%
%   PROGRESSBAR(X) with X being a number in [0;1] updates the last progress
%   bar to show 100*X percent task completion.
%
%   Example:
%      progressbar('Generating output ...')
%      for i = 0.01 : 0.01 : 1
%          progressbar(i)
%          pause(0.1)
%      end
%
%   See also WAITBAR.

% Copyright Alexander Schaefer
%
% PROGRESSBAR is inspired by TEXTPROGRESSBAR by Paul Proteus:
% https://www.mathworks.com/matlabcentral/fileexchange/28067-
% text-progress-bar

%% Validate input.
% Check number of input arguments.
narginchk(1, 1)

% Check type of input argument.
if ~(isnumeric(x) || ischar(x))
    error('X must be a number in [0;1] or a message string.')
end

% For numeric input, verify the value.
if isnumeric(x) && (x<0 || x>1)
    error('X must be in [0;1].')
end

%% Preprocessing.
% Set message and progress variables.
prg = 0;
prg(isnumeric(x)) = x(isnumeric(x));
msg(ischar(x),:) = x(ischar(x),:);

% Activate command line log.
diaryfile = 'pgbdiary.tmp';
diary(diaryfile)

% Set the width of the progress bar.
width = 70;

%% Clear last progress bar.
% Clear the last progress bar if there is no message to print.
if isempty(msg)
    % Try to open the command line log.
    fid = fopen(diaryfile);
    if fid ~= -1
        % Read the last line from the command line window.
        readline = fgets(fid);
        while ischar(readline)
            lastline = readline;
            readline = fgets(fid); 
        end
        fclose(fid);
        
        % Check if the last line is a progress bar.
        pattern = ['\[[\# ]{', num2str(width-2), '}\] [\d\s]{3}\%'];
        if ~isempty(regexpi(lastline, pattern))
            % Delete the last line.
            fprintf(repmat('\b', 1, width+6));
        end
    end
end

%% Print message and progress bar.
% Print the message.
if ~isempty(msg)
    fprintf([msg, '\n'])
end

% Print the progress bar.
nHash = round(prg * (width-2));
hashStr = repmat('#', 1, nHash);
spaceStr = repmat(' ', 1, width-2-nHash);
fprintf(['[', hashStr, spaceStr, '] %3u%%\n'], round(100*prg));

end

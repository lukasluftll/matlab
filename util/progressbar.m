function progressbar(x)
% PROGRESSBAR Display progress bar in command window.
%   PROGRESSBAR(X) creates, updates, or closes a progress bar displayed in
%   the command window.
%   
%   If PROGRESSBAR(X) is called for the first time, X must be a message
%   string that tells the user what he is waiting for.
%
%   To update the progress bar, X is a number in interval [0;1] that 
%   represents the progress of the computation.
%
%   To close the progress bar, X is either not provided or another message 
%   string.
%
%   Example:
%      progressbar('Generating output ...')
%      for i = 0 : 0.01 : 1
%          progressbar(i)
%          pause(0.1)
%      end
%      progressbar
%
%   See also WAITBAR.

% Copyright Alexander Schaefer
%
% PROGRESSBAR is inspired by TEXTPROGRESSBAR by Paul Proteus:
% https://www.mathworks.com/matlabcentral/fileexchange/28067-
% text-progress-bar

%% Validate input.
% Check number of input arguments.
nargin(0, 1)

% Depending on whether the progress bar has been initialized, check
% validity of input argument.
if isempty(init) % Not yet initialized.
    if ~ischar(x)
        error('X must be a string when calling PROGRESSBAR first.')
    end
else % Already initialized.
    if isnumeric(x)
        if x < 0 || x > 1
            error('X must be element of [0;1].')
        end
    elseif ~ischar(x)
        
            

% Check type of input argument.
if nargin > 0
    if (isnumeric(x) && (x<0 || x>1)) || ~ischar(x)
        error('X must be a number in [0;1] or a string.')
    end
end

%% Set parameters.
% Persistently stores whether progressbar has been initialized.
persistent init;

% Width of progress bar.
width = 75;


if isempty(strCR) && ~ischar(x),
    % Progress bar must be initialized with a string
    error('The text progress must be initialized with a string');
elseif isempty(strCR) && ischar(x),
    % Progress bar - initialization
    fprintf('%s',x);
    strCR = -1;
elseif ~isempty(strCR) && ischar(x),
    % Progress bar  - termination
    strCR = [];  
    fprintf([x '\n']);
elseif isnumeric(x)
    % Progress bar - normal progress
    x = floor(x);
    percentageOut = [num2str(x) '%%'];
    percentageOut = [percentageOut repmat(' ',1,strPercentageLength-length(percentageOut)-1)];
    nDots = floor(x/100*strDotsMaximum);
    dotOut = ['[' repmat('.',1,nDots) repmat(' ',1,strDotsMaximum-nDots) ']'];
    strOut = [percentageOut dotOut];
    
    % Print it on the screen
    if strCR == -1,
        % Don't do carriage return during first run
        fprintf(strOut);
    else
        % Do it during all the other runs
        fprintf([strCR strOut]);
    end
    
    % Update carriage return
    strCR = repmat('\b',1,length(strOut)-1);
    
else
    % Any other unexpected input
    error('Unsupported argument type');
end
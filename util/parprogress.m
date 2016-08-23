function percent = parprogress(N)
%parprogress Progress monitor (progress bar) that works with parfor.
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

error(nargchk(0, 1, nargin, 'struct'));

if nargin < 1
    N = -1;
end

percent = 0;
w = 50; % Width of progress bar

if N > 0
    f = fopen('parprogress.txt', 'w');
    if f<0
        error('Do you have write permissions for %s?', pwd);
    end
    fprintf(f, '%d\n', N); % Save N at the top of progress.txt
    fclose(f);
    
    if nargout == 0
        disp(['  0%[>', repmat(' ', 1, w), ']']);
    end
elseif N == 0
    delete('parprogress.txt');
    percent = 100;
    
    if nargout == 0
        disp([repmat(char(8), 1, (w+9)), char(10), '100%[', repmat('=', 1, w+1), ']']);
    end
else
    if ~exist('parprogress.txt', 'file')
        error('parprogress.txt not found. Run parprogress(N) before parprogress to initialize parprogress.txt.');
    end
    
    f = fopen('parprogress.txt', 'a');
    fprintf(f, '1\n');
    fclose(f);
    
    f = fopen('parprogress.txt', 'r');
    progress = fscanf(f, '%d');
    fclose(f);
    percent = (length(progress)-1)/progress(1)*100;
    
    if nargout == 0
        perc = sprintf('%3.0f%%', percent); % 4 characters wide, percentage
        disp([repmat(char(8), 1, (w+9)), char(10), perc, '[', repmat('=', 1, round(percent*w/100)), '>', repmat(' ', 1, w - round(percent*w/100)), ']']);
    end
end

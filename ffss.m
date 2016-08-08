function pred = ffss(pf, prev, sigma)

%% Validate input.
% Check number of input arguments.
narginchk(3, 3)

sp = size(pf.Particles);
if ~all(size(prev) == sp)
    error(['PREV must be a matrix of size ', num2str(sp(1)), ...
        'x', num2str(sp(2)), '.'])
end

if ~all(size(sigma) == [sp(2), sp(2)])
    error(['SIGMA must be a matrix of size ', num2str(sp(2)), ...
        'x', num2str(sp(2)), '.'])
end

if ~isspd(sigma)
    error('SIGMA must be a symmetric positive-definite matrix.')
end

pred = mvnrnd(prev, sigma);

pred(:,4:6) = wrapToPi(pred(:,4:6));

end

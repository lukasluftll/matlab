function m = quatmean(q, w)

nargoutchk(0, 1)
narginchk(1, 2)
validateattributes(q, {'numeric'}, {'real', 'rows', 4}, '', 'Q')
nq = size(q,2);
if nargin < 2
    w(1:nq) = 1/nq;
end
validateattributes(w, {'numeric'}, ...
    {'real', 'positive', 'rows', 1, 'cols', nq}, '', 'W')

M = zeros(4);
for i = 1 : size(q,2)
    M = M + w(i) * q(:,i) * q(:,i)';
end

[m,~] = eigs(M,1);

end

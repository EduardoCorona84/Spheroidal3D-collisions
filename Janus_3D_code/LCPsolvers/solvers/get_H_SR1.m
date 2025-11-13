function [H, opts] = get_H_SR1(s, y, opts)

error('Currently this impementation is unfinished')
% if ~exist('bMask', 'var') || isempty(bMask)
%% TODO Allow this to be used
    bMask = false(n,1);
% end

[r, h0, ~, S, Y] = updateQNMemory(k, s, y, opts, bMask);
U = zeros(n, r);
V = zeros(n, r);
H = eye(n) * h0;
B = eye(n) / h0;
for j = 1:r 
    s = S(:,j); 
    y = Y(:,j);
    % update H (and mem)
    denomH = dot(s - H*y,y);
    sigma = sign(denomH);
    denomH = sqrt(sigma*denomH);
    assert(isreal(denomH));
    u = (s - H*y) / denomH;
    H = H + sigma*(u*u');
    % update H (and mem)
    denomB = dot(y - B*s,s);
    lambda = sign(denomB);
    denomB = sqrt(lambda*denomB);
    assert(isreal(denomB));
    v = (y - B*s) / denomB;
    B = B + lambda*(v*v');
    % Store for proof of concept
    if lambda > 0 
        U(:,j) = v;
    else
        V(:,j) = v;
    end
end
mask = arrayfun(@(i) all(U(:,i) == 0),1:r);
U = U(:,~mask);
V = V(:,mask);
assert(norm(B*H - eye(n))/norm(H)/norm(B) < 1e-8)
% TODO: Make this matrix free...
H = @(x) H*x;
% We may have an indefinite matrix from the SR1 update
if ~isempty(U)  
    % TODO do this better with a shifted power method
    assert(min(eig(H)) > 0)
end

end
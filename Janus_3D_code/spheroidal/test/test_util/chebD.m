function [x, D] = chebD(N)
% Chebyshev grid and first-derivative matrix on [-1,1].
if N == 0
    x = 1;
    D = 0;
    return;
end
k = (0:N).';
x = cos(pi * k / N);
c = [2; ones(N - 1, 1); 2] .* (-1).^k;
X = repmat(x, 1, N + 1);
dX = X - X.';
D = (c * (1 ./ c).') ./ (dX + eye(N + 1));
D = D - diag(sum(D, 2));
end
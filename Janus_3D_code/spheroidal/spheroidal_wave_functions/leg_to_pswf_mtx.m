function [D, lambda, deg] = leg_to_pswf_mtx(m, p, gamma, parity)
%{
parity \in {0, 1}
%}

assert(gamma ~= 0, 'Gamma should be non-zero; zero case is not implemented.');

nu = m + parity;
N = floor((p - nu)/2) + 1;
k = (0:N-1).';
deg = nu + 2*k; % Legendre degrees used in this block (order is m)

%% Construct tridiagonal matrix
% See Falloon et al. (2003).
A = zeros(N,1); B = zeros(N,1); C = zeros(N,1);
A(2:end) = -gamma^2 .* ((nu-m+2*k(2:end)-1).*(nu-m+2*k(2:end))) ...
                    ./ ((2*nu+4*k(2:end)-3).*(2*nu+4*k(2:end)-1));
B(:)     = (deg).*(deg+1) ...
         - 2*gamma^2 .* ((deg).*(deg+1) + m^2 - 1) ...
          ./ ((2*deg-1).*(2*deg+3));
C(1:end-1) = -gamma^2 .* ((nu+m+2*k(1:end-1)+1).*(nu+m+2*k(1:end-1)+2)) ...
                      ./ ((2*nu+4*k(1:end-1)+3).*(2*nu+4*k(1:end-1)+5));

% Form a symmetric tridiagonal matrix
w = ones(N,1);
for j = 2:N
    w(j) = w(j-1) * C(j-1) / A(j);
end

t = sqrt(abs(A(2:end).*C(1:end-1)));
T = diag(B) + diag(t,1) + diag(t,-1);

% Eigendecomposition
[V, L] = eig(T);
[lambda, perm] = sort(diag(L), 'ascend');
V = V(:, perm);

% Undo scaling
D = V ./ sqrt(w);

% Have to normalize so the spheroidal wave functions are normalized
% Recall that:
%   ||P_l^m||^2 = 2/(2l+1) * (l+m)!/(l-m)!
normalization_constant = 2./(2*deg+1) .* factorial(deg+abs(m)+1) ./ factorial(deg-abs(m)+1);

for j = 1:N
    s2 = sum( (D(:,j).^2) .* normalization_constant);
    D(:,j) = D(:,j) / sqrt(s2);
end

end %%% END MAIN FUNCTION

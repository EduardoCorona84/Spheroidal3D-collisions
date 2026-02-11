function [D, lambda, deg] = leg_to_pswf_mtx(m, p, gamma, parity)
%{
Build the Legendre-to-PSWF coefficient matrix for fixed order |m|
and a fixed parity block.

Inputs:
    m      - order
    p      - maximum Legendre degree used for truncation
    gamma  - spheroidal parameter
    parity - parity block selector:
             0 -> degrees m, m+2, m+4, ...
             1 -> degrees m+1, m+3, m+5, ...

Outputs:
    D      - change-of-basis matrix from PSWF coefficients to
             associated-Legendre coefficients within this parity block
    lambda - angular spheroidal eigenvalues for this block (ascending)
    deg    - Legendre degrees represented in this block
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

% Preserve complex phase in the symmetrized off-diagonal.
t = sqrt(A(2:end).*C(1:end-1));

% For imaginary gamma, the opposite square-root branch matches cprofcn's convention.
% If one wants to remove this, then one must also deal with the convention in Rnm1 and Rnm3.
if real(gamma^2) < 0
    t = -t;
end
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
normalization_constant = 2./(2*deg+1) .* factorial(deg+abs(m)) ./ factorial(deg-abs(m));

% Normalize each PSWF column in the Legendre basis so that
%   sum_k D(k,j)^2 * ||P_{deg_k}^m||^2 = 1.
for j = 1:N
    s2 = sum( (D(:,j).^2) .* normalization_constant);
    D(:,j) = D(:,j) / sqrt(s2);
end

end %%% END MAIN FUNCTION

function [ASWFnm, lambda_nm] = ASWFnm(n, m, v, phi, c, p)
%{
Angular spheroidal wave function S_n^m(v) e^{i m phi}, v in (-1,1).
Computed via Legendre series using coefficients from leg_to_pswf_mtx,
but evaluated using Ynm to avoid manual Legendre summation.

Inputs:
  n, m   - degree/order (|m| <= n)
  v      - evaluation points in (-1,1)
  phi    - azimuth
  c      - spheroidal parameter
  p      - max Legendre degree for truncation (optional)

Outputs:
  ASWFnm    - column vector of S_n^m(v) e^{i m phi}
  lambda_nm - eigenvalue for (n,m) if requested
%}

if nargin < 5
    error('ASWFnm requires at least n, m, v, phi, c.');
end
if nargin < 6 || isempty(p)
    p = n + max(8, ceil(8 * abs(c)));
end

if isscalar(phi)
    phi = repmat(phi, numel(v), 1);
elseif numel(phi) ~= numel(v)
    error('phi must be scalar or the same length as v.');
end

if abs(m) > n
    ASWFnm = zeros(numel(v), 1);
    if nargout > 1, lambda_nm = []; end
    return;
end

u = real(acos(v));

mm = abs(m);

if c == 0
    ASWFnm = Ynm(n, m, u, phi);
    if nargout > 1, lambda_nm = n * (n + 1); end
    return;
end

if p < n
    error('p must be >= n.');
end

parity = mod(n - mm, 2);
nu = mm + parity;
j = (n - nu) / 2 + 1;
if j < 1 || j > (floor((p - nu) / 2) + 1)
    error('The max Legendre degree p is too small for the given (n,m).');
end

[D, lambda] = leg_to_pswf_mtx(mm, p, c, parity);
coeffs = D(:, j);

degs = (nu:2:p).';
num_terms = numel(degs);
ASWFnm = zeros(numel(v), 1);
for k = 1:num_terms
    l = degs(k);
    ASWFnm = ASWFnm + coeffs(k) * Ynm(l, m, u, phi);
end

if nargout > 1
    lambda_nm = lambda(j);
end
end

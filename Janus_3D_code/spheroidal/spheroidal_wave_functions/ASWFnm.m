function ASWFnm = ASWFnm(n, m, v, phi, gamma, p, iopnorm)
%{
Angular spheroidal wave function S_n^m(v) e^{i m phi}, v in (-1,1).
Computed via Legendre series using coefficients from leg_to_pswf_mtx,
but evaluated using Ynm to avoid manual Legendre summation.

Inputs:
    n, m    - degree/order (|m| <= n). If m = [], returns all orders m=-n:n.
    v       - evaluation points in (-1,1)
    phi     - azimuth
    c       - spheroidal parameter
    p       - max Legendre degree for truncation (optional)
    iopnorm - normalization flag (optional):
              0 -> associated-Legendre (Meixner-Schafke) normalization
                \int |S_n^m|^2 dv = \frac{2}{2n+1}\frac{(n+|m|)!}{(n-|m|)!}
              1 -> unit normalization

Outputs:
    ASWFnm - column vector of S_n^m(v) e^{i m phi), or
             a matrix with columns for m=-n:n if m = []
%}

if nargin < 5, error('ASWFnm requires at least n, m, v, phi, c.'); end
if nargin < 6, p = n + max(8, ceil(8 * abs(gamma))); end
if nargin < 7, iopnorm = 0; end

if isscalar(phi)
    phi = repmat(phi, numel(v), 1);
elseif numel(phi) ~= numel(v)
    error('phi must be scalar or the same length as v.');
end

u = real(acos(v));

if isempty(m)
    m_list = -n:n;
    if gamma == 0
        ASWFnm = Ynm(n, [], u, phi);
        for idx = 1:numel(m_list)
            m_val = m_list(idx);
            ASWFnm(:, idx) = ASWFnm(:, idx) / LOCAL_ylm_norm(n, m_val);
            if iopnorm == 1
                ASWFnm(:, idx) = ASWFnm(:, idx) / LOCAL_ms_norm(n, m_val);
            end
        end
        return;
    end

    ASWFnm = zeros(numel(v), numel(m_list));

    if p < n
        error('p must be >= n.');
    end

    % Precompute coefficients per m.
    coeffs_by_m = cell(2 * n + 1, 1);
    for idx = 1:numel(m_list)
        m_val = m_list(idx);
        mm = abs(m_val);
        parity = mod(n - mm, 2);
        nu = mm + parity;
        j = (n - nu) / 2 + 1;
        D = leg_to_pswf_mtx(mm, p, gamma, parity);
        coeffs = D(:, j);
        coeffs_by_m{idx} = coeffs;
    end

    for l = 0:p
        Yn_all = Ynm(l, [], u, phi); % columns for m=-l:l
        for idx = 1:numel(m_list)
            m_val = m_list(idx);
            mm = abs(m_val);
            parity = mod(n - mm, 2);
            nu = mm + parity;
            if l < nu || mod(l - nu, 2) ~= 0
                continue;
            end
            k = (l - nu) / 2 + 1;
            coeffs = coeffs_by_m{idx};
            if k > numel(coeffs)
                continue;
            end
            term = Yn_all(:, m_val + l + 1);
            term = term / LOCAL_ylm_norm(l, m_val);
            ASWFnm(:, idx) = ASWFnm(:, idx) + coeffs(k) * term;
        end
    end
    if iopnorm == 0
        for idx = 1:numel(m_list)
            ASWFnm(:, idx) = LOCAL_ms_norm(n, m_list(idx)) * ASWFnm(:, idx);
        end
    end
    return;
end

[ASWFnm] = LOCAL_single_m(n, m, u, phi, gamma, p, iopnorm);
end %% END MAIN FUNCTION

function ASWFnm = LOCAL_single_m(n, m, u, phi, c_eval, p, iopnorm)
    if abs(m) > n
        ASWFnm = zeros(numel(u), 1);
        return;
    end

    mm = abs(m);

    if c_eval == 0
        ASWFnm = Ynm(n, m, u, phi);
        ASWFnm = ASWFnm / LOCAL_ylm_norm(n, m);
        if iopnorm == 1
            ASWFnm = ASWFnm / LOCAL_ms_norm(n, m);
        end
        return;
    end

    if p < n, error('p must be >= n.'); end

    parity = mod(n - mm, 2);
    nu = mm + parity;
    j = (n - nu) / 2 + 1;
    if j < 1 || j > (floor((p - nu) / 2) + 1)
        error('The max Legendre degree p is too small for the given (n,m).');
    end

    D = leg_to_pswf_mtx(mm, p, c_eval, parity);
    coeffs = D(:, j);

    degs = (nu:2:p).';
    num_terms = numel(degs);
    ASWFnm = zeros(numel(u), 1);
    for k = 1:num_terms
        l = degs(k);
        term = Ynm(l, m, u, phi);
        term = term / LOCAL_ylm_norm(l, m); % Remove the normalization constant
        ASWFnm = ASWFnm + coeffs(k) * term;
    end
    if iopnorm == 0
        ASWFnm = LOCAL_ms_norm(n, m) * ASWFnm; % Match expected L2 normalization
    end
end

function nrm = LOCAL_ylm_norm(l, m)
    % Ynm normalization constant for degree/order (l,m).
    mm = abs(m);
    nrm = sqrt((2*l + 1) / (4*pi) * factorial(l - mm) / factorial(l + mm));
end

function nrm = LOCAL_ms_norm(n, m)
    % Meixner-Schafke constant
    mm = abs(m);
    nrm = sqrt(2 / (2*n + 1) * factorial(n + mm) / factorial(n - mm));
end

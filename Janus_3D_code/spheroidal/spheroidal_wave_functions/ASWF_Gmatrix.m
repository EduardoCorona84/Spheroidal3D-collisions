function G = ASWF_Gmatrix(p, u0, gamma, fig, oblate, vnp, iopnorm)
%{
Construct the SWF basis-change matrix analogous to Gmatrix.

This matrix maps coefficients in the weighted basis
    S_n^m(v)e^{im phi} / w(v)
to coefficients in the unweighted basis
    S_n^m(v)e^{im phi}.

Inputs
    p       - truncation order (integer >= 0)
    u0      - spheroidal shape parameter:
    gamma   - spheroidal parameter (scalar)
    fig     - optional plotting flag (default 0)
    oblate  - optional shape flag (default false)
    vnp     - optional number of Gauss-Legendre nodes (default 50)
    iopnorm - optional ASWFnm normalization flag (default 0; must be 0 or 1)

Output
    G       - (p+1)^2 x (p+1)^2 basis-change matrix

For each m, we define matrices
    A_m(i,j) = <S_i^m, S_j^m>
    B_m(i,j) = <S_i^m, S_j^m / w(v)>
where
    w(v) = sqrt(u0^2 - v^2)  (prolate)
    w(v) = sqrt(u0^2 + v^2)  (oblate)
where <.,.> is approximated with GL quadrature in v and exact 2*pi in phi (due to the integral in phi).
Then A_m * G_m = B_m, so G_m maps weighted-basis coefficients to
unweighted-basis coefficients:
    a_m = G_m * b_m.
Note that G should be block diagonal by the value of m due to the integral in phi.
%}

if nargin < 3
    error('ASWF_Gmatrix requires at least p, u0, gamma.');
end
if nargin < 4, fig = 0; end
if nargin < 5, oblate = false; end
if nargin < 6, vnp = 50; end
if nargin < 7, iopnorm = 0; end

[vp, vw] = g_grid(vnp);
vp = vp(:);
vw = vw(:);

sp = (p + 1)^2;
G = complex(zeros(sp));

% Precompute all degree blocks S_n^m(v,phi=0) once.
% This avoids repeated ASWFnm calls inside the m-loop and lets us slice out
% the needed column for each (n,m).
Sn_all = cell(p + 1, 1);
for n = 0:p
    Sn_all{n + 1} = ASWFnm(n, [], vp, 0, gamma, p, iopnorm);
end

geti = @(n, m) m + n.^2 + n + 1;
% Quadrature weights for the inner products.
qw_unweighted = vw;
if oblate
    qw_weighted = vw ./ sqrt(u0.^2 + vp.^2);
else
    qw_weighted = vw ./ sqrt(u0.^2 - vp.^2);
end

% The basis is block-diagonal in m, so each m can be assembled independently.
for m = -p:p
    nlist = (abs(m):p).';
    Nm = numel(nlist);
    idx = geti(nlist, m);

    % Columns are S_n^m(v) for n = |m|, ..., p at this fixed m.
    ASWF_m_block = complex(zeros(vnp, Nm));
    for t = 1:Nm
        n = nlist(t);
        Sn = Sn_all{n + 1};
        ASWF_m_block(:, t) = Sn(:, m + n + 1);
    end

    % Gram blocks:
    %   Am: unweighted SWF inner products
    %   Bm: inner products against weighted basis functions
    Am = 2 * pi * (ASWF_m_block' * (qw_unweighted .* ASWF_m_block));
    Bm = 2 * pi * (ASWF_m_block' * (qw_weighted .* ASWF_m_block));

    % Solve A_m * G_m = B_m for the per-m change-of-basis matrix.
    G(idx, idx) = Am \ Bm;
end

if fig
    figure;
    % Sparsity pattern of numerically non-negligible entries.
    spy(abs(G) > 1e-10);
end
end

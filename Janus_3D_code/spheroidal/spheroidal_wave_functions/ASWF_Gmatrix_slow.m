function G = ASWF_Gmatrix_slow(p, u0, gamma, fig, oblate, vnp, iopnorm)
%{
Reference (intentionally inefficient) implementation of ASWF_Gmatrix.

This mirrors the full i,j assembly style of Gmatrix.m:
1) Build full projection matrices A and B entry-by-entry.
2) Solve G = A \ B.

The goal is correctness/reference comparison against ASWF_Gmatrix, not speed.
%}

if nargin < 3
    error('ASWF_Gmatrix_slow requires at least p, u0, gamma.');
end
if nargin < 4, fig = 0; end
if nargin < 5, oblate = false; end
if nargin < 6, vnp = 50; end
if nargin < 7, iopnorm = 0; end

if oblate
    error('Oblate case is not implemented.');
end

[vp, vw] = g_grid(vnp);
vp = vp(:);
vw = vw(:);
weight_denom = sqrt(u0.^2 - vp.^2);

sp = (p + 1)^2;
ii = (1:sp).';
nn = floor(sqrt(ii - 1));
mm = ii - nn.^2 - nn - 1;

A = complex(zeros(sp));
B = complex(zeros(sp));

for i = 1:sp
    ni = nn(i);
    mi = mm(i);
    for j = 1:sp
        nj = nn(j);
        mj = mm(j);

        if mi ~= mj
            A(i, j) = 0;
            B(i, j) = 0;
            continue;
        end

        % Deliberately evaluate both basis functions per (i,j) pair,
        % to mimic the naive Gmatrix implementation.
        Si = ASWFnm(ni, mi, vp, 0, gamma, p, iopnorm);
        Sj = ASWFnm(nj, mj, vp, 0, gamma, p, iopnorm);

        integrand_unweighted = conj(Si) .* Sj;
        integrand_weighted = integrand_unweighted ./ weight_denom;

        A(i, j) = 2 * pi * sum(vw .* integrand_unweighted);
        B(i, j) = 2 * pi * sum(vw .* integrand_weighted);
    end
end

G = A \ B;

if fig
    figure;
    spy(abs(G) > 1e-10);
end
end

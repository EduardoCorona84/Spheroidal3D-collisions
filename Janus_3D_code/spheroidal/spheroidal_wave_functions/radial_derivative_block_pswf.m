function dFr = radial_derivative_block_pswf(p, u_region, c, ioprad)
%{
For fixed |m| and u, one cprofcn_mex call returns all n >= |m| values
of d/du R_n^m. We reuse that vector and place it into both +m and -m
columns instead of calling the wrappers per (n,m).

Inputs
    p        - truncation order
    u_region - target radial coordinates
    c        - spheroidal parameter
    ioprad   - cprofcn radial kind flag
                1 -> first kind  (R^(1))
                2 -> third kind  (R^(3) = R^(1) + i R^(2))

Output
dFr      - nt x (p+1)^2 derivative block where columns follow
            geti(n,m) = n^2 + n + m + 1 ordering.
%}

sp = (p + 1)^2;
u_region = u_region(:);
nt = numel(u_region);

% Derivative block in packed geti(n,m) indexing.
dFr = zeros(nt, sp);
if nt == 0
    return;
end

% cprofcn is evaluated at conj(c) for imag(c)<0, then conjugated back.
% This matches the branch convention used by Rnm wrappers.
neg_imag = imag(c) < 0;
cc = c;
if neg_imag
    cc = conj(c);
end

% Avoid repeated cprofcn calls at identical targets (i.e. for self-evaluation).
[u_unique, ~, map] = unique(u_region, 'stable');
nu = numel(u_unique);
dFr_unique = zeros(nu, sp);

% Precompute column indices for each |m| block.
nvals_by_mm = cell(p + 1, 1);
idx_pos_by_mm = cell(p + 1, 1);
idx_neg_by_mm = cell(p + 1, 1);
for mm = 0:p
    % For fixed |m|=mm, valid degrees are n=mm,...,p.
    nvals = (mm:p).';
    nvals_by_mm{mm + 1} = nvals;
    % Stacked columns for +m and -m.
    idx_pos_by_mm{mm + 1} = nvals.^2 + nvals + mm + 1;
    if mm > 0
        idx_neg_by_mm{mm + 1} = nvals.^2 + nvals - mm + 1;
    else
        idx_neg_by_mm{mm + 1} = [];
    end
end

for mm = 0:p
    % cprofcn returns entries of size (p-mm+1) for fixed mm.
    lnum = p - mm + 1;
    nvals = nvals_by_mm{mm + 1};
    if isempty(nvals)
        continue;
    end
    idx_pos = idx_pos_by_mm{mm + 1};
    idx_neg = idx_neg_by_mm{mm + 1};

    % Map degree n to row index in cprofcn output: n-mm+1.
    nm_pairs = nvals - mm + 1;
    for ku = 1:nu
        % cprofcn radial argument is x1 = u-1.
        x1 = u_unique(ku) - 1;
        if ioprad == 1
            [~, ~, r1dc, ir1de] = cprofcn_mex(cc, mm, lnum, ioprad, x1, 0, 0, 0);
            dr = r1dc(:) .* 10.^ir1de(:);
        elseif ioprad == 2
            [~, ~, r1dc, ir1de, ~, ~, r2dc, ir2de] = ...
                cprofcn_mex(cc, mm, lnum, ioprad, x1, 0, 0, 0);
            dr = r1dc(:) .* 10.^ir1de(:) + 1i * r2dc(:) .* 10.^ir2de(:);
        else
            error('Unsupported ioprad input.');
        end

        if neg_imag
            dr = conj(dr);
        end

        % d/du R_n^{+m}(u) = d/du R_n^{-m}(u)
        dval = dr(nm_pairs);
        dFr_unique(ku, idx_pos) = dval;
        if mm > 0
            dFr_unique(ku, idx_neg) = dval;
        end
    end
end

% Restore original u ordering (including duplicates).
dFr = dFr_unique(map, :);
end

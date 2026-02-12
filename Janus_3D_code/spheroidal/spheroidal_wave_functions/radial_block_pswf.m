function Fr = radial_block_pswf(p, u_region, c, ioprad)
%{
For fixed |m| and u, one cprofcn_mex call returns all n >= |m| values.
We reuse that vector and place it into both +m and -m columns instead
of calling the wrappers per (n,m).

Inputs
    p        - truncation order
    u_region - target radial coordinates
    c        - spheroidal parameter
    ioprad   - cprofcn radial kind flag
                1 -> first kind  (R^(1))
                2 -> third kind  (R^(3) = R^(1) + i R^(2))

Output
Fr       - nt x (p+1)^2 block where columns follow
            geti(n,m) = n^2 + n + m + 1 ordering.
%}

sp = (p + 1)^2;
u_region = u_region(:);
nt = numel(u_region);

% n=0,m=0 starts as 1; remaining columns are overwritten below.
Fr = ones(nt, sp);
if nt == 0
    return;
end

% cprofcn is evaluated at conj(c) for imag(c)<0, then conjugated back.
% This matches the branch convention used by Rnm1/Rnm3 wrappers.
neg_imag = imag(c) < 0;
cc = c;
if neg_imag
    cc = conj(c);
end

% Precompute column indices for each |m| block.
nvals_by_mm = cell(p + 1, 1);
idx_pos_by_mm = cell(p + 1, 1);
idx_neg_by_mm = cell(p + 1, 1);
for mm = 0:p
    % For fixed |m|=mm, valid degrees are n=mm,...,p.
    nvals = (mm:p).';
    nvals_by_mm{mm + 1} = nvals;
    % Packed columns for +m and -m.
    idx_pos_by_mm{mm + 1} = nvals.^2 + nvals + mm + 1;

    if mm > 0
        idx_neg_by_mm{mm + 1} = nvals.^2 + nvals - mm + 1;
    else % mm = 0
        idx_neg_by_mm{mm + 1} = [];
    end
end

for mm = 0:p
    % cprofcn returns entries of size (p-mm+1) for fixed mm.
    % For example, if p = 8, mm = 3, then lnum = 6 corresponds to
    % entries correspond to n = 3,4,5,6,7,8.
    lnum = p - mm + 1;
    nvals = nvals_by_mm{mm + 1};
    if isempty(nvals)
        continue;
    end
    idx_pos = idx_pos_by_mm{mm + 1};
    idx_neg = idx_neg_by_mm{mm + 1};

    % Map degree n to row index in cprofcn output: n-mm+1.
    nm_pairs = nvals - mm + 1; % Picks out correct (n, m)
    for ku = 1:nt
        % cprofcn radial argument is x1 = u-1.
        x1 = u_region(ku) - 1;
        if ioprad == 1
            [r1c, ir1e] = cprofcn_mex(cc, mm, lnum, ioprad, x1, 0, 0, 0);
            r = r1c(:) .* 10.^ir1e(:);
        elseif ioprad == 2
            [r1c, ir1e, ~, ~, r2c, ir2e] = cprofcn_mex(cc, mm, lnum, ioprad, x1, 0, 0, 0);
            r = r1c(:) .* 10.^ir1e(:) + 1i * r2c(:) .* 10.^ir2e(:);
        else
            error('Unsupported ioprad value.');
        end

        if neg_imag
            r = conj(r);
        end

        % R_n^{+m}(u) = R_n^{-m}(u)
        rval = r(nm_pairs);
        Fr(ku, idx_pos) = rval;
        if mm > 0
            Fr(ku, idx_neg) = rval;
        end
    end
end
end

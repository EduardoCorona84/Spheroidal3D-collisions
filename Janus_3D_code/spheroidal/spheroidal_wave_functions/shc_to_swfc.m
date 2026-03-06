function swfc = shc_to_swfc(shc, p, gamma)
%{
Convert spherical-harmonic coefficients to angular SWF coefficients.

Inputs:
    shc    - sp x C array of shAna/Ynm coefficients, sp = (p + 1)^2
    p      - truncation order
    gamma  - spheroidal parameter (same convention as ASWFnm)

Output:
    swfc   - sp x C array of SWF coefficients
%}

if nargin < 3 || isempty(gamma)
    gamma = 0;
end

sp = (p + 1)^2;
if size(shc, 1) ~= sp
    if isvector(shc) && numel(shc) == sp
        shc = reshape(shc, sp, 1);
    else
        error("shc first dimension must be (p+1)^2.");
    end
end
orig_size = size(shc);
shc_mat = reshape(shc, sp, []);

[ylm_norm, ms_norm] = LOCAL_index_norms(p);
% Convert normalized Y_n^m coefficients to associated-Legendre coefficients.
leg_coeffs = shc_mat .* ylm_norm;

if gamma == 0
    swfc_mat = leg_coeffs;
else
    swfc_mat = LOCAL_legendre_to_swf_unitnorm(leg_coeffs, p, gamma) ./ ms_norm;
end

swfc = reshape(swfc_mat, orig_size);
end

function swfc = LOCAL_legendre_to_swf_unitnorm(leg_coeffs, p, gamma)
    %{
    Map associated-Legendre coefficients g_{n,m} to SWF coefficients s_{n,m}
    with the convention that the spheroidal wave functions are normalized.
    
    For fixed |m|, leg_to_pswf_mtx returns a matrix D such that
      g_block = D * s_block
    within one parity block (n = |m|,|m|+2,... or n = |m|+1,|m|+3,...).
    This routine splits each m-block into even/odd parity, solves D \ g for
    each block, and interleaves the two results back into expected ordering.
    %}
    sp = (p + 1)^2;
    num_cols = size(leg_coeffs, 2);

    % Precompute Legendre -> SWF coefficient matrices for each |m|.
    D_even = cell(p + 1, 1);
    D_odd = cell(p + 1, 1);
    for mm = 0:p
        D_even{mm + 1} = leg_to_pswf_mtx(mm, p, gamma, 0);
        if p >= mm + 1
            D_odd{mm + 1} = leg_to_pswf_mtx(mm, p, gamma, 1);
        else
            D_odd{mm + 1} = [];
        end
    end

    geti = @(n, m) m + n.^2 + n + 1;
    swfc = zeros(sp, num_cols);

    for m = -p:p
        mm = abs(m);
        degrees = (mm:p).';
        idx = geti(degrees, m);
        g = leg_coeffs(idx, :);

        % Split parity for n = m,m+2,... and n = m+1,m+3,...
        ge = g(1:2:end, :);
        go = g(2:2:end, :);

        if ~isempty(ge)
            De = D_even{mm + 1};
            se = De \ ge;
        else
            se = zeros(0, num_cols);
        end
        if ~isempty(go)
            Do = D_odd{mm + 1};
            so = Do \ go;
        else
            so = zeros(0, num_cols);
        end

        s_m = zeros(numel(degrees), num_cols);
        s_m(1:2:end, :) = se;
        s_m(2:2:end, :) = so;
        swfc(idx, :) = s_m;
    end
end

function [ylm_norm, ms_norm] = LOCAL_index_norms(p)
    sp = (p + 1)^2;
    ii = (1:sp).';
    nn = floor(sqrt(ii - 1));
    mm = ii - nn.^2 - nn - 1;
    absm = abs(mm);

    ratio = factorial(nn - absm) ./ factorial(nn + absm);
    ylm_norm = sqrt(((2 * nn + 1) ./ (4 * pi)) .* ratio);
    ms_norm = sqrt((2 ./ (2 * nn + 1)) ./ ratio);
end

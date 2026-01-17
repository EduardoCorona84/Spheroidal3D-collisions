function S = pswf(n, m, c, p, u)
%{
    Implements the prolate spheroidal wave function by using 
    a sum of Associated Legendre functions.

    Inputs
        n, m - order/degree
        c - parameter for PSWF
        p - max degree to use for the Legendre functions
        u - evaluation points

    Outputs
        S - prolate spheroidal wave function evaluated at u > 1

    TODO:
    Allow the cache to store the D matrices.
%}

mm = abs(m);
parity = mod(n - mm, 2);
nu = mm + parity;
j = (n - nu)/2 + 1; % column index for (n, m)

D = leg_to_pswf_mtx(mm, p, c, parity); % columns are coeffs a_n^m
L = legendre_otc(p, u, 0, 0, 0); L = L{1};

% \sum_{n, nu \leq 2n+\nu \leq p} a_n^m P_{2n + \nu}^m(c; u)
S = L(geti(nu:2:p, m)).'*D(:, j);

end
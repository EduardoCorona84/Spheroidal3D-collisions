function [R1, dR1, R3, dR3] = modified_laplace_test_radial_mode(p, n, m, u, c, oblate, mex_opts)
if nargin < 6 || isempty(oblate)
    oblate = false;
end
if nargin < 7
    mex_opts = struct();
end

idx = m + n^2 + n + 1;
[R1_all, dR1_all, R3_all, dR3_all] = ...
    modified_laplace_eval_radial_sphwv(p, u, c, oblate, mex_opts);

R1 = R1_all(:, idx);
dR1 = dR1_all(:, idx);
R3 = R3_all(:, idx);
dR3 = dR3_all(:, idx);
end

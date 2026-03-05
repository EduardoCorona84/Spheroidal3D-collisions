function [R1, dR1, R3, dR3] = modifiedLaplaceEvalRadialSphwv(p, u, c, oblate, mex_opts)
% Returns matrices of size N x sp.
if nargin < 5 || isempty(mex_opts)
    mex_opts = struct('precision_bits', 256, 'output_digits', 17, 'max_memory', 2000);
end

n_vec = 0:p;
[u_unique, iu] = LOCAL_unique_u_grid(u);

if oblate
    if ~isreal(c) || ~(real(c) > 0)
        error('Oblate sphwv radial evaluation requires real scalar c > 0.');
    end
    [R1u, R3u, dR1u, dR3u] = obl_radial_r13_pure_imag_mex(n_vec, 'all', real(c), u_unique, mex_opts);
else
    if abs(real(c)) > 1e-12
        error('Prolate sphwv radial evaluation requires purely imaginary scalar c.');
    end
    if imag(c) < 0
        error('Magnitude of parameter is negative: this should never happen.')
    end

    gamma = abs(imag(c));
    [R1u, R3u, dR1u, dR3u] = pro_radial_r13_pure_imag_mex(n_vec, 'all', gamma, u_unique, mex_opts);
end
R1 = R1u(iu, :);
dR1 = dR1u(iu, :);
R3 = R3u(iu, :);
dR3 = dR3u(iu, :);
end

function [u_unique, iu] = LOCAL_unique_u_grid(u)
    u = u(:);
    if isempty(u)
        u_unique = u;
        iu = zeros(0, 1);
        return;
    end
    [u_unique, ~, iu] = uniquetol(u, 1e-12);
end

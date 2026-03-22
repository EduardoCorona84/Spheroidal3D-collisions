function [R1, dR1, R3, dR3] = modified_laplace_eval_radial_sphwv(p, u, c, oblate, mex_opts)
if nargin < 5 || isempty(mex_opts)
    mex_opts = struct();
end
mex_opts = modified_laplace_get_mex_options(mex_opts);

n_vec = 0:p;
[u_unique, u_to_unique_idx] = LOCAL_unique_u_grid(u);

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
R1 = R1u(u_to_unique_idx, :);
dR1 = dR1u(u_to_unique_idx, :);
R3 = R3u(u_to_unique_idx, :);
dR3 = dR3u(u_to_unique_idx, :);
end %% END MAIN FUNCTION

function [u_unique, u_to_unique_idx] = LOCAL_unique_u_grid(u)
    u = u(:);
    if isempty(u)
        u_unique = u;
        u_to_unique_idx = zeros(0, 1);
        return;
    end
    [u_unique, ~, u_to_unique_idx] = uniquetol(u, 1e-12);
end

function [soln, truesoln, sigma_vec, condK, info] = modified_laplace_bie_problem(p, lambda, u0, target_distance, plt, problem_type, rng_seed)
%{
Theoretical background on the BIEs for Yukawa's equation can be found
in Bryan Quaife's PhD thesis.

Helper for a single-body modified-Laplace exterior boundary value problem,
using interior point charges/point fluxes to generate true solutions.

Supports two exterior problem types:
    (1) Dirichlet, combined layer ansatz:
        (0.5*I + SLP + DLP) * sigma = g
    (2) Neumann, SLP ansatz:
        (-0.5*I + Kp) * sigma = f
        where Kp is the adjoint double-layer (dS/dn).

Boundary data is generated analytically from the interior point charges:
    g = u_true|_S,        for Dirichlet
    f = (dn u_true)|_S,   for Neumann.

Inputs:
    p               - spectral order.
    lambda          - Yukawa parameter.
    u0              - spheroid parameter.
    target_distance - exterior offset along the normal.
    plt             - true/false for plotting.
    problem_type    - optional:
                      "exterior_dirichlet" (default) or "exterior_neumann".
    rng_seed        - optional scalar seed for random interior charge placement.

Outputs:
    soln      - reconstructed potential at exterior targets.
    truesoln  - true point-charge Yukawa potential at targets.
    sigma_vec - recovered boundary density.
    condK     - condition number of on-surface BIE matrix.
    info      - diagnostic struct.
%}

if nargin < 7 || isempty(rng_seed)
    rng_seed = 2;
end
if ~isscalar(u0)
    error("This helper only supports a single spheroid (scalar u0).");
end
if any(target_distance <= 0)
    error("target_distance must be positive for exterior target evaluation.");
end

np = 2*p*(p + 1);

params = SpheroidalParameters;
params.isReal = false;
params.u0 = u0;
params.a = 1 / u0;
params.oblate = false;
params.centers = [0 0 0];
params.thetas = 0;
params.phis = 0;
params.sigma = zeros(np, 1);

Y = params.get_X();
nu = get_norm_vecs(p, params.u0, params.oblate);
Xeval = Y + target_distance * nu;

% Interior point charges/point fluxes are used to generate the data on the surface.
num_point_charges = 4;
charge_radius = 0.25 * params.a * sqrt(max(params.u0^2 - 1, eps));
[ptch, Xptch] = LOCAL_place_interior_point_charges( ...
    params.a, params.u0, params.centers(1, :), num_point_charges, charge_radius, rng_seed);
truesolnSurf = LOCAL_yukawa_pt_charge(ptch, Xptch, Y, lambda);
truefluxSurf = LOCAL_yukawa_pt_charge_flux(ptch, Xptch, Y, nu, lambda);

params_op = copy(params);
params_op.sigma = eye(np);
params_op.get_shc();

switch problem_type
    case "exterior_dirichlet"
        rhs = truesolnSurf;

        SL = spheroidalModifiedSLP(params_op, lambda, []);
        DL = spheroidalModifiedDLP(params_op, lambda, []);
        K = 0.5 * eye(np) + SL + DL;

        sigma_vec = K \ rhs;

        params_sigma = copy(params);
        params_sigma.sigma = sigma_vec;
        params_sigma.get_shc();

        SL_eval = spheroidalModifiedSLP(params_sigma, lambda, Xeval);
        DL_eval = spheroidalModifiedDLP(params_sigma, lambda, Xeval);
        soln = SL_eval + DL_eval;
    case "exterior_neumann"
        rhs = truefluxSurf;

        Kp = spheroidalModifiedSP(params_op, lambda, []);
        K = -0.5 * eye(np) + Kp;

        sigma_vec = K \ rhs;

        params_sigma = copy(params);
        params_sigma.sigma = sigma_vec;
        params_sigma.get_shc();

        soln = spheroidalModifiedSLP(params_sigma, lambda, Xeval);
    otherwise
        error('Invalid problem type.');
end

condK = cond(K);
truesoln = LOCAL_yukawa_pt_charge(ptch, Xptch, Xeval, lambda);

if isrow(truesoln), truesoln = truesoln.'; end
if isrow(soln), soln = soln.'; end

rel_soln_abs = norm(soln - truesoln);
rel_soln_err = rel_soln_abs / norm(truesoln);

info = struct();
info.problem_type = problem_type;
info.ptch = ptch;
info.Xptch = Xptch;
info.truesolnSurf = truesolnSurf;
info.truefluxSurf = truefluxSurf;
info.Xeval = Xeval;
info.rel_soln_abs = rel_soln_abs;
info.rel_soln_err = rel_soln_err;

if plt
    figure;
    point_err = abs(soln - truesoln);
    scatter3(Xeval(:, 1), Xeval(:, 2), Xeval(:, 3), 30, log10(point_err + eps), 'filled');
    hold on;
    scatter3(Xptch(:, 1), Xptch(:, 2), Xptch(:, 3), 150, ptch, 'filled', 'MarkerEdgeColor', 'k');
    hold off;
    axis equal;
    grid on;
    title('log10 pointwise error');
    colorbar;
end
end %% END MAIN FUNCTION

function [ptch, Xptch] = LOCAL_place_interior_point_charges(a, u0, center, n_charges, charge_radius, rng_seed)
    minor_axis = a * sqrt(max(u0^2 - 1, eps));
    center = reshape(center, 1, []);

    % Keep all interior charges in a ball strictly smaller than the minor axis.
    charge_radius = min(charge_radius, 0.5 * minor_axis);

    % Randomly place point charges inside a sphere around the spheroid center.
    rng(rng_seed);
    dirs = randn(n_charges, 3);
    dirs = dirs ./ max(sqrt(sum(dirs.^2, 2)));
    radii = charge_radius * rand(n_charges, 1).^(1/3);
    Xptch = dirs .* radii + repmat(center, n_charges, 1);

    ptch = 2 * rand(n_charges, 1) - 1;
end

function pcp = LOCAL_yukawa_pt_charge(ptch, Xptch, Y, lambda)
    M = numel(ptch);
    np = size(Y, 1);
    pcp = zeros(np, 1);
    for i = 1:M
        r = sqrt(sum((Xptch(i, :) - Y).^2, 2));
        pcp = pcp + ptch(i) * exp(-lambda * r) ./ (4 * pi * r);
    end
end

function flux = LOCAL_yukawa_pt_charge_flux(ptch, Xptch, Y, Ny, lambda)
    M = numel(ptch);
    np = size(Y, 1);
    flux = zeros(np, 1);
    for i = 1:M
        r_vec = Xptch(i, :) - Y;
        r = sqrt(sum(r_vec.^2, 2));
        r = max(r, eps);
        edotn = sum(Ny .* (r_vec ./ r), 2);
        flux = flux + ptch(i) * exp(-lambda * r) .* (1 ./ r.^2 + lambda ./ r) .* edotn / (4 * pi);
    end
end

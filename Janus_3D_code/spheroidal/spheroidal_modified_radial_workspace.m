function radial_ws = spheroidal_modified_radial_workspace(params, lambda, X_trg, mex_opts)
%{
Build shared radial data for modified-Laplace operators.

The modified-Laplace DLP and DP both need the same spheroidal radial
functions. This helper batches all radial requests and returns a
workspace/data structure that both operators can utilize.

The workspace has two layers of indexing:
1. "groups" are gathered based on the quantities that
    determine which radial MEX backend result is needed.
2. "targets" stores, for each spheroid, how its local surface/target
    requests map back for a spheroid family.

This function also maintains a short-lived persistent cache so the same
workspace can be reused across the DLP and DP builds inside one operator
update cycle. Note that the cache is cleared when the geometry/operators
are updated: see SpheroidalMS_UpdateOperators.m.
%}

persistent radial_cache

if nargin >= 1 && ischar(params)
    switch char(params)
        case 'clear'
            % Clear the persistent workspace cache between update cycles.
            radial_cache = [];
            radial_ws = [];
            return;
        otherwise
            error('Unsupported spheroidal_modified_radial_workspace operation.');
    end
end

if nargin < 4
    mex_opts = struct();
end
mex_opts = modified_laplace_get_mex_options(mex_opts);

% Expand scalar geometry parameters into one entry per spheroid so the rest
% of the file can treat everything uniformly.
[p, u0, a, oblate, ns] = LOCAL_extract_params_info(params);


cache_key = LOCAL_cache_key(p, u0, a, oblate, X_trg);
% Reuse a previously built workspace if the geometry/target was already computed for.
if ~isempty(radial_cache)
    match = find(strcmp({radial_cache.key}, cache_key), 1);
    if ~isempty(match)
        radial_ws = radial_cache(match).workspace;
        radial_ws.cache_hit = true;
        radial_ws.cache_key = cache_key;
        return;
    end
end

% Tolerance so that we can collapse "the same" u values.
tol = 9e-12;

% Build the spheroidal c-parameters used by the angular/radial basis.
[c_ang, c_rad] = LOCAL_build_c_params(lambda, a, oblate);

% "groups" holds one entry per distinct radial family, meaning one unique
% combination of p, c_rad, and oblate/prolate type. Any spheroids in the
% same family can share a single radial evaluation block.
groups = struct( ...
    'key', {}, ...
    'c_rad', {}, ...
    'c_ang', {}, ...
    'oblate', {}, ...
    'surface_indices', {}, ...
    'all_u', {}, ...
    'u_unique', {}, ...
    'R1', {}, ...
    'dR1', {}, ...
    'R3', {}, ...
    'dR3', {}, ...
    'unique_u_count', {});

% "targets" holds per-spheroid bookkeeping: which group it belongs to, what
% its target coordinates are in spheroidal form, and which rows of the
% group-level radial arrays correspond to its local surface/target requests.
targets = repmat(struct( ...
        'has_targets', false, ...
        'group_index', 0, ...
        'surface_group_u_index', [], ...
        'target_group_u_indices', [], ...
        'S', [], ...
        'u_x', [], ...
        'v_x', [], ...
        'phi_x', [], ...
        'indices_interior', [], ...
        'indices_surface', [], ...
        'indices_exterior', [], ...
        'nt', 0 ...
    ), 1, ns);
group_keys = cell(1, ns);

% Temporary per-spheroid storage used by the second pass below:
% 1) u_local_unique_by_surface{k} stores the unique u values requested
% by spheroid k.
% 2) local_to_unique_by_surface{k} maps the original local request ordering
% back into that unique list.
%
% These are the only pieces of intermediate state needed after the global
% grouped u-sets have been formed.
u_local_unique_by_surface = cell(1, ns);
local_to_unique_by_surface = cell(1, ns);

% Helper function for what follows.
function group_index = LOCAL_find_or_add_group()
    % Group by radial physics only. If multiple spheroids share the same
    % group key, they will reuse one radial evaluation block.
    group_index = find(strcmp({groups.key}, group_keys{k}), 1);
    if isempty(group_index)
        group_index = numel(groups) + 1;
        groups(group_index).key = group_keys{k};
        groups(group_index).c_rad = c_rad(k);
        groups(group_index).c_ang = c_ang(k);
        groups(group_index).oblate = oblate(k);
        groups(group_index).surface_indices = [];
        groups(group_index).all_u = [];
    end
end

% Identify which shared radial group this spheroid belongs to.
for k = 1:ns
    group_keys{k} = LOCAL_group_key(p, c_rad(k), oblate(k));
    targets(k).group_index = LOCAL_find_or_add_group();

    if isempty(X_trg)
        % Self-evaluation: only the surface radius u0 is needed.
        u_request = u0(k);
        targets(k).has_targets = false;
        targets(k).nt = 0;
    else
        % Explicit target evaluation: collect the per-spheroid target cloud.
        Xt_k = X_trg(:,:,k); % Target points for spheroid k
        nt_k = size(Xt_k, 1);
        if nt_k == 0
            % Empty target block still needs the surface u index for consistency.
            u_request = u0(k);
            targets(k).has_targets = false;
            targets(k).nt = 0;
        else
            S = cart2spheroidal(Xt_k, a(k), oblate(k));
            u_x = S(:, 1);
            v_x = real(S(:, 2));
            phi_x = S(:, 3);

            targets(k).has_targets = true;
            targets(k).S = S;
            targets(k).u_x = u_x;
            targets(k).v_x = v_x;
            targets(k).phi_x = phi_x;
            targets(k).indices_interior = (u_x < u0(k) - tol);
            targets(k).indices_surface = (abs(u_x - u0(k)) <= tol);
            targets(k).indices_exterior = (u_x > u0(k) + tol);
            targets(k).nt = nt_k;

            % Keep the surface radius in the first slot so code can later
            % recover the surface u index with request_group_u_indices(1)
            % and target indices with request_group_u_indices(2:end).
            u_request = [u0(k); u_x(:)];
        end
    end

    % Remove redundant u-values within this spheroid first.
    [u_local_unique, ~, local_to_unique] = uniquetol(u_request(:), tol);
    u_local_unique_by_surface{k} = u_local_unique;
    local_to_unique_by_surface{k} = local_to_unique;

    groups(targets(k).group_index).all_u = [groups(targets(k).group_index).all_u; u_local_unique(:)];
    groups(targets(k).group_index).surface_indices(end+1) = k;
end

% Collapse all u requests from every spheroid in this radial family, then
% evaluate the radial MEX exactly once on that unique set.
for g = 1:numel(groups)
    group = groups(g);
    [u_unique, ~, ~] = uniquetol(group.all_u(:), tol);
    [R1, dR1, R3, dR3] = LOCAL_eval_group_radial(p, u_unique, group.c_rad, group.oblate, mex_opts);
    groups(g).u_unique = u_unique(:);
    groups(g).R1 = R1;
    groups(g).dR1 = dR1;
    groups(g).R3 = R3;
    groups(g).dR3 = dR3;
    groups(g).unique_u_count = numel(u_unique);
    groups(g).all_u = [];
end

% Map this spheroid's local unique u values into the group-level unique
% array, then expand back to the original local request ordering.
for k = 1:ns
    g = targets(k).group_index;
    local_unique_group_u_indices = LOCAL_lookup_group_u_indices(u_local_unique_by_surface{k}, groups(g).u_unique, tol);
    request_group_u_indices = local_unique_group_u_indices(local_to_unique_by_surface{k});

    % By construction, the first local request is always the source surface u0.
    targets(k).surface_group_u_index = request_group_u_indices(1);
    targets(k).target_group_u_indices = request_group_u_indices(2:end);
end

% Package everything needed by the optimized DLP/DP evaluators. The returned
% struct is intentionally explicit because it also doubles as a profiling and
% debugging data structure for tests.
radial_ws = struct( ...
    'p', p, ...
    'ns', ns, ...
    'u0', u0, ...
    'a', a, ...
    'oblate', oblate, ...
    'c_ang', c_ang, ...
    'c_rad', c_rad, ...
    'groups', groups, ...
    'targets', targets, ...
    'cache_key', cache_key, ...
    'cache_hit', false, ...
    'group_count', numel(groups), ...
    'unique_u_counts', reshape([groups.unique_u_count], 1, []) ...
);

% Store the completed workspace in the persistent cache for reuse by the
% other modified-Laplace operator built in the same update cycle.
cache_entry = struct('key', cache_key, 'workspace', radial_ws);
if isempty(radial_cache)
    radial_cache = cache_entry;
else
    radial_cache(end+1) = cache_entry;
end
end %% END MAIN FUNCTION

function [p, u0, a, oblate, ns] = LOCAL_extract_params_info(params)
% Expand scalar geometry data to one entry per spheroid.
    p = params.p;
    u0 = params.u0;
    a = params.a;
    oblate = params.oblate;
    ns = max(size(params.sigma, 3), 1);

    if isscalar(u0)
        u0 = u0 .* ones(1, ns);
    else
        u0 = reshape(u0, 1, []);
    end

    if isscalar(a)
        a = a .* ones(1, ns);
    else
        a = reshape(a, 1, []);
    end

    if isscalar(oblate)
        oblate = oblate .* ones(1, ns);
    else
        oblate = reshape(oblate, 1, []);
    end
end

function [c_ang, c_rad] = LOCAL_build_c_params(lambda, a, oblate)
% Convert the modified-Laplace parameter lambda into the convention used by
% the spheroidal wave-function backends.
    c_ang = zeros(size(a));
    c_rad = zeros(size(a));
    for k = 1:numel(a)
        if oblate(k)
            c_rad(k) = lambda * a(k);
        else
            c_rad(k) = 1j * lambda * a(k);
        end
        c_ang(k) = c_rad(k);
    end
end

function cache_key = LOCAL_cache_key(p, u0, a, oblate, X_trg)
% Build the persistent workspace cache key.
%
% This key intentionally ignores mex_opts because the cache is scoped to a
% single run configuration and explicitly cleared by the caller between
% operator-update cycles.
    if isempty(X_trg)
        target_repr = '[]';
    else
        target_repr = mat2str(X_trg(:).', 17);
    end

    cache_key = sprintf('p=%d|u0=%s|a=%s|oblate=%s|targets=%s', ...
        p, mat2str(u0, 17), mat2str(a, 17), mat2str(double(oblate), 17), target_repr);
end

function group_key = LOCAL_group_key(p, c_rad, oblate)
    group_key = sprintf('p=%d|cre=%.17g|cim=%.17g|ob=%d', ...
        p, real(c_rad), imag(c_rad), oblate);
end

function [R1, dR1, R3, dR3] = LOCAL_eval_group_radial(p, u_unique, c_rad, oblate, mex_opts)
% Evaluate one grouped radial problem on the unique-u set.
    n_vec = 0:p;
    if oblate
        [R1, R3, dR1, dR3] = obl_radial_r13_pure_imag_mex(n_vec, 'all', real(c_rad), u_unique, mex_opts);
    else
        [R1, R3, dR1, dR3] = pro_radial_r13_pure_imag_mex(n_vec, 'all', abs(imag(c_rad)), u_unique, mex_opts);
    end
end

function group_u_indices = LOCAL_lookup_group_u_indices(u_query, u_unique, tol)
% Map a local u list back into a group-level unique-u array.
    group_u_indices = zeros(numel(u_query), 1);
    for i = 1:numel(u_query)
        match = find(abs(u_unique - u_query(i)) <= tol, 1);
        if isempty(match)
            error('Failed to map local u value into radial workspace.');
        end
        group_u_indices(i) = match;
    end
end

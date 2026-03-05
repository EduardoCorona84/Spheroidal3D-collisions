function Y = SSph_MatVec(V, L, params, modlap_opts)
%{
Matvec handler for the mobility solver. The main purpose of this is to decode
the inputs so that we can pass it to the correct functions.

Note that this is modeled after VSh_MatVec_RB2.m.

IMPORTANT NOTE: This does not implement the principal-valued part only;
this implements the entire BIO on-surface (so there is an associated jump 
relation that is being added).

Inputs
params - parameter struct with fields such as: 
    p - (int) - spherical harmonic degree
    n3 - (int) - number of objects 
    kerd - (int) - kernel dimension 
    dense - (bool) - Dense vs FMM far interactions
    out - (bool) exterior or interior problem
    C - (double n_c x 3) - centers in original box [a,b]^3
    rd -  (double n_c x 1) - sphere radii
    flag_pot - (string) - the following cases are supported: 
        'SL_L_3D'     single layer, Laplace
        'dSL_L_3D'    normal derivative of single layer, Laplace
        'DL_L_3D'     double layer, Laplace
        'dDL_L_3D'    normal derivative of double layer, Laplace
        'SL_LMOD_3D'  single layer, modified Laplace (Yukawa)
        'dSL_LMOD_3D' normal derivative of single layer, modified Laplace
        'DL_LMOD_3D'  double layer, modified Laplace
        'dDL_LMOD_3D' normal derivative of double layer, modified Laplace
        'SL_Stk_3D'   single layer, Stokes
        'DL_Stk_3D'   double layer, Stokes
        'dSL_Stk_3D'  normal derivative of single layer, Stokes
        'TSL_Stk_3D'  traction kernel of single layer, Stokes
    tsl_dealiasing - (bool, optional) enable dealiasing in optimized TSL near-evaluation
    tsl_dealiasing_pad - (int, optional) spherical harmonic padding for dealiasing (default 4)
    Xv - (double 3*np*n3 x 3) duplicated list of points; seems to be only
    used for Kernel_Eval/FMM
    X - (double 3*np*n3 x 3) list of points; should be in the global frame
    (i.e. the points are already rotated and in their respective centered
    positions)
V - (double n_c*N_deg x 1) array of density or densities
    should (probably) be formatted like
    [sig_x(p_1) ; sig_y(p_1) ; sig_z(p_1) ; sig_x(p_2) ; ...]
L - (sparse array) Extra matrix for completion flow / nullspace correction
modlap_opts - (struct)
    optional struct to pass in for the MEX calls used in modified Laplace

Outputs
Y - (double) 3*N_trg x 1 array
    also is in interleaved format; see the comment for the input V.
%}

persistent Gmatrix_cache

%% SETUP
p  = params.p;      % Spharm degree p 
n3 = params.n3;     % Number of objects
kerd = params.kerd; % Kernel dimension
dense = params.dense; % Dense vs FMM off diagonal 
C = params.C;       % object centers
out = params.out;   % outside vs inside sphere

equ_radii = params.equ_radii;
polar_radii = params.polar_radii;

pot = params.flag_pot;
is_modified_laplace = contains(pot, 'LMOD');
if is_modified_laplace
    lambda = params.lambda;
    if nargin < 4 || isempty(modlap_opts)
        modlap_opts = params.sphwv_mex_opts;
    end
else
    lambda = [];
end

if is_modified_laplace && ~dense
    error('Modified Laplace currently requires params.dense = true.');
end

[tsl_dealiasing_flag, tsl_dealiasing_pad] = LOCAL_get_tsl_dealiasing_options(params);

if isempty(Gmatrix_cache)
    body_shape_types = params.shape_type;
    Gmatrix_cache = cell(1, n3);
    for i=1:n3
        body_shape = body_shape_types(i);
        switch body_shape
            case 'sphere'
                    continue
            otherwise
                [u0, ~] = calculate_u0_and_a_from_radii(body_shape, equ_radii(i), polar_radii(i));
                Gmatrix_cache{i} = sparse(Gmatrix(params.p, u0, 0, strcmp(body_shape, 'oblate')));
        end
    end
end

% Rotation of bodies
if isfield(params,'MRot')
    MRot = params.MRot;
    rot = true;
else
    MRot = cell(1,n3);  
    rot = false;
end

% Handle neighbors
if isfield(params,'neigh')
    neigh=params.neigh; % Neighbor list (cell(n3,1)) 
else
    distC = get_distances_between_centers(C);  
    neigh = cell(n3,1);

    % Circumscribe a sphere around each body
    for i=1:n3
        max_radius_i = max(equ_radii(i), polar_radii(i));
        neigh{i} = find(distC(i,:)<max_radius_i*params.mdist); 
    end
    params.neigh=neigh;  
end

% Determine whether Kernel_Eval needs target normals
nortrg = true;
switch pot(1:3)
    case {'DL_','SDL','dDL'}
        nortrg = false;
end

np = 2*p*(p+1); 
Nb = kerd*np; % 3*np
N = Nb*n3; % Number of "data points" needed for entire system
X = params.Xp;
Xv = params.X;
Nor = params.nor;
W2 = params.W2;

if kerd == 1 % For this, get the unstacked versions
    Xv = params.Xp;
    Nor = params.Nrp;
    W2 = params.Wg.';
end

%% ACTUAL MATVEC
if strcmp(V, 'Mat')
    % Start with Kernel Eval (correct for far-interactions)
    % Then, replace self-to-self and self-to-near appropriately.
    params_KE = params; % Copy it to handle this for kerd = 1 (maybe not so efficient)
    params_KE.W2 = W2;
    params_KE.nor = Nor;
    if strncmp(pot, 'dDL', 3)
        params_KE.targnor = Nor;
    end
    Y = Kernel_Eval(Xv,Xv,params_KE);

    if kerd == 3
        prm = zeros(1, Nb);
        prm(1:np) = 1:3:Nb;
        prm(np+1:2*np) = 2:3:Nb;
        prm(2*np+1:3*np) = 3:3:Nb;
        iprm = zeros(1, Nb);
        iprm(prm) = 1:Nb;
    else
        iprm = 1:Nb;
    end

    for body_ind=1:n3
        % Indices for source particle
        I_box = (1:Nb)+Nb*(body_ind-1);

        % Number of neighbors for source particle (i.e. number of particles
        % that are needed for near-evaluation).
        num_ngh = length(neigh{body_ind});

        % Build neighbor indices
        neigh{body_ind} = reshape(neigh{body_ind},1,[]);         
        I_nghv = repmat((1:Nb)',1,num_ngh)+Nb*(repmat(neigh{body_ind},Nb,1)-1);  
        I_nghv = I_nghv(:);

        equ_radius = equ_radii(body_ind);
        polar_radius = polar_radii(body_ind);
        body_shape_type = params.shape_type(body_ind);
        source_Gmatrix = [];
        if ~isempty(Gmatrix_cache)
            source_Gmatrix = Gmatrix_cache{body_ind};
        end

        target_pts = Xv(I_nghv(1:kerd:end),:);
        target_normals = Nor(I_nghv(1:kerd:end),:);

        % Near-interaction
        if rot
            target_pts = (target_pts - C(body_ind, :)) * MRot{body_ind};
            target_normals = target_normals * MRot{body_ind};
        else
            target_pts = (target_pts - repmat(C(body_ind,:), size(target_pts,1), 1));
        end

        if strcmp(body_shape_type, 'prolate') || strcmp(body_shape_type, 'oblate')
            [u0, a, oblate] = LOCAL_calculate_u0_a(equ_radius, polar_radius, body_shape_type);
            if num_ngh > 1
                slf = find(neigh{body_ind}==body_ind); 
                indv_off = [1:Nb*(slf-1) (Nb*slf+1):Nb*num_ngh].';
                ind_off  = [1:np*(slf-1) (np*slf+1):np*num_ngh].';

                Ynear = zeros(Nb*num_ngh, Nb);

                if ~isempty(ind_off)
                    target_pts_off = target_pts(ind_off,:);
                    target_normals_off = target_normals(ind_off,:);
                    Ynear(indv_off,:) = LOCAL_spheroid_near_matrix( ...
                        pot, u0, a, oblate, target_pts_off, target_normals_off, ...
                        np, iprm, source_Gmatrix, lambda, modlap_opts, ...
                        tsl_dealiasing_flag, tsl_dealiasing_pad);
                end
                Ynear(Nb*(slf-1)+1:Nb*slf,:) = LOCAL_spheroid_near_matrix( ...
                    pot, u0, a, oblate, [], [], np, iprm, source_Gmatrix, lambda, modlap_opts, ...
                    tsl_dealiasing_flag, tsl_dealiasing_pad);
            else
                Ynear = LOCAL_spheroid_near_matrix( ...
                    pot, u0, a, oblate, [], [], np, iprm, source_Gmatrix, lambda, modlap_opts, ...
                    tsl_dealiasing_flag, tsl_dealiasing_pad);
            end
        elseif strcmp(body_shape_type, 'sphere')
            switch pot(1:3)
                case 'SL_'
                    pMat = 'SMat';
                case 'SDL'
                    pMat = 'SDMat';
                case 'dSL'
                    pMat = 'SpMat';
                case 'TSL'
                    pMat = 'TSMat';
                case 'DL_'
                    pMat = 'DMat';
                case 'dDL'
                    pMat = 'DpMat';
                case 'TDL'
                    pMat = 'TDMat';
                otherwise
                    error('Unsupported pot for spherical near-eval.');
            end

            % Rescale geometry to unit sphere for spectral near-eval
            Xtrg = (1/equ_radius) * target_pts;

            [th, phi, rho] = cart2sph(Xtrg(:,1), Xtrg(:,2), Xtrg(:,3));
            th(th<0) = th(th<0) + 2*pi;
            phi = pi/2 - phi;

            if kerd == 1
                Nrtrg = target_normals;
                if is_modified_laplace
                    Ynear = Sh_Mod_Kernel_Eval_off([p 1], pMat, 0, out, rho, phi, th, Nrtrg, lambda);
                else
                    Ynear = Sh_Kernel_Eval_off([p 1], pMat, 0, out, rho, phi, th, Nrtrg);
                end
            else
                Nrtrg = Nor(I_nghv,:);
                if rot
                    Nrtrg = Nrtrg * MRot{body_ind};
                end

                Ynear = Vsh_Kernel_Eval_off([p 1], pMat, 0, out, rho, phi, th, Nrtrg);

                if rot
                    Ynear = permute(reshape(Ynear,[3 np*num_ngh Nb]),[2 3 1]); 
                    Ynear = reshape(Ynear,[],3)*MRot{body_ind};   
                    Ynear = permute(reshape(Ynear,[np*num_ngh Nb 3]),[3 1 2]); 
                    Ynear = reshape(Ynear,[],Nb); 
                end
            end

            if strncmp(pot,'SL_',3)
                rda = equ_radius; % SLP scales like r
            elseif strncmp(pot,'dDL',3)
                rda = 1/equ_radius; % dDL scales like 1/r
            else
                rda = 1;
            end
            Ynear = rda * Ynear;
        else
            error('Invalid body shape type; should be "sphere" or "prolate" or "oblate".');
        end

        if rot && kerd==3 && (strcmp(body_shape_type,'prolate') || strcmp(body_shape_type,'oblate'))
            Ynear = reshape(Ynear, [], 3) * MRot{body_ind}; % Reformat vector into matrix and rotate back
            Ynear = reshape(Ynear, 3*size(target_pts,1), Nb); % Place back into a vector
        end

        Y(I_nghv, I_box) = Ynear;
    end

    if out
        ct = 0.5;
    else
        ct = -0.5;
    end

    % Add jump relation on the diagonal
    if ~strcmp(pot(1:3),'SL_') && ~strcmp(pot(1:3),'dDL')
        if isfield(params,'a')
            if strcmp(pot(2:3),'SL') 
                if isfield(params,'eta') && strcmp(pot(1:3),'dSL') 
                    Y = params.eta*Y + (params.a+ct)*eye(N); 
                else
                    Y = Y + (params.a+ct)*eye(N); 
                end
            else
                Y = Y + (params.a-ct)*eye(N); 
            end
        end
    end

    % Add nullspace term
    if ~isempty(L)
        Y = Y + L;
    end
elseif ~isempty(V) && isnumeric(V)
    % Actually do the matvec given a numeric input
    V = reshape(V,N,[]); % TODO: Remove this. Not sure why this is needed.
    Y = zeros(N,size(V,2));

    for body_ind=1:n3
        % Indices for source particle
        I_box = (1:Nb)+Nb*(body_ind-1);

        % Number of neighbors for source particle (i.e. number of particles
        % that are needed for near-evaluation).
        num_ngh = length(neigh{body_ind});

        % Build neighbor sphere index
        neigh{body_ind} = reshape(neigh{body_ind},1,[]);         
        I_nghv = repmat((1:Nb)',1,num_ngh)+Nb*(repmat(neigh{body_ind},Nb,1)-1);  
        I_nghv = I_nghv(:);

        equ_radius = equ_radii(body_ind);
        polar_radius = polar_radii(body_ind);
        body_shape_type = params.shape_type(body_ind);
        source_Gmatrix = [];
        if ~isempty(Gmatrix_cache)
            source_Gmatrix = Gmatrix_cache{body_ind};
        end

        target_pts = Xv(I_nghv(1:kerd:end),:);
        target_normals = Nor(I_nghv(1:kerd:end),:);

        %% Near-interaction
        % No matter what, we need to rotate the target spheres and spheroids to be in the local frame
        % of the current body.
        if rot % It seems that the target points are already in the local frame. Is this a good idea?
            target_pts = (target_pts - C(body_ind, :)) * MRot{body_ind};
            target_normals = target_normals * MRot{body_ind};
        else
            % Translate to local frame (no rotation)
            target_pts = (target_pts - repmat(C(body_ind,:), size(target_pts,1), 1));
        end

        if strcmp(body_shape_type, 'prolate') || strcmp(body_shape_type, 'oblate')
            if kerd==1
                % Scalar Laplace near-eval on spheroids
                sigma = V(I_box,:);

                % Build local parameters
                [u0, a, oblate] = LOCAL_calculate_u0_a(equ_radius, polar_radius, body_shape_type);
                params_i = SpheroidalParameters();
                params_i.sigma = sigma; % set p for geometry generation
                params_i.u0 = u0;
                params_i.a = a;
                params_i.oblate = oblate;
                params_i.centers = [0 0 0];
                params_i.thetas = 0;
                params_i.phis = 0;
                params_i.Rmat = eye(3);

                if num_ngh > 1
                    slf = find(neigh{body_ind}==body_ind); 
                    indv_off = [1:Nb*(slf-1) (Nb*slf+1):Nb*num_ngh].';
                    ind_off  = [1:np*(slf-1) (np*slf+1):np*num_ngh].';

                    Ynear = zeros(Nb*num_ngh, size(sigma,2));

                    if ~isempty(ind_off)
                        target_pts_off = target_pts(ind_off,:);
                        target_normals_off = target_normals(ind_off,:);
                        Ynear(indv_off,:) = LOCAL_eval_scalar_spheroidal_potential( ...
                            pot, params_i, target_pts_off, target_normals_off, ...
                            sigma, source_Gmatrix, lambda, modlap_opts);
                    end

                    Ynear(Nb*(slf-1)+1:Nb*slf,:) = LOCAL_eval_scalar_spheroidal_potential( ...
                        pot, params_i, [], [], sigma, source_Gmatrix, lambda, modlap_opts);
                else
                    Ynear = LOCAL_eval_scalar_spheroidal_potential( ...
                        pot, params_i, [], [], sigma, source_Gmatrix, lambda, modlap_opts);
                end
            else
                % Grab density on source particle
                sig_x = V(I_box(1:3:end));
                sig_y = V(I_box(2:3:end));
                sig_z = V(I_box(3:3:end));

                % Build local parameters for optimized L2Stk
                [u0, a, oblate] = LOCAL_calculate_u0_a(equ_radius, polar_radius, body_shape_type);
                params_i = SpheroidalParameters();
                params_i.sigma = sig_x; % set p for geometry generation
                params_i.u0 = u0;
                params_i.a = a;
                params_i.oblate = oblate;
                params_i.centers = [0 0 0];
                params_i.thetas = 0;
                params_i.phis = 0;
                params_i.Rmat = eye(3);

                if num_ngh > 1
                    % Get indices of target sources (self vs off-diagonal within neighbor list)
                    slf = find(neigh{body_ind}==body_ind); 
                    indv_off = [1:Nb*(slf-1) (Nb*slf+1):Nb*num_ngh].';
                    ind_off  = [1:np*(slf-1) (np*slf+1):np*num_ngh].';

                    Ynear = zeros(Nb*num_ngh,1);

                    if ~isempty(ind_off)
                        target_pts_off = reshape(target_pts(ind_off,:), [], 3, 1);
                        target_normals_off = reshape(target_normals(ind_off,:), [], 3, 1);
                        switch pot
                            case 'SL_Stk_3D'
                                Ynear(indv_off) = LOCAL_eval_l2stk_slp(params_i, target_pts_off, sig_x, sig_y, sig_z, source_Gmatrix);
                            case 'TSL_Stk_3D'
                                Ynear(indv_off) = LOCAL_eval_l2stk( ...
                                    params_i, target_pts_off, target_normals_off, sig_x, sig_y, sig_z, source_Gmatrix, ...
                                    tsl_dealiasing_flag, tsl_dealiasing_pad);
                            case 'DL_Stk_3D'
                                error('No reason to use this.');
                            otherwise
                                error('Incorrect potential passed in.');
                        end
                    end

                    switch pot
                        case 'SL_Stk_3D'
                            Ynear(Nb*(slf-1)+1:Nb*slf) = LOCAL_eval_l2stk_slp(params_i, [], sig_x, sig_y, sig_z, source_Gmatrix);
                        case 'TSL_Stk_3D'
                            Ynear(Nb*(slf-1)+1:Nb*slf) = LOCAL_eval_l2stk( ...
                                params_i, [], [], sig_x, sig_y, sig_z, source_Gmatrix, ...
                                tsl_dealiasing_flag, tsl_dealiasing_pad);
                        case 'DL_Stk_3D'
                            error('No reason to use this.');
                        otherwise
                            error('Incorrect potential passed in.');
                    end
                else
                    switch pot
                        case 'SL_Stk_3D'
                            Ynear = LOCAL_eval_l2stk_slp(params_i, [], sig_x, sig_y, sig_z, source_Gmatrix);
                        case 'TSL_Stk_3D'
                            Ynear = LOCAL_eval_l2stk( ...
                                params_i, [], [], sig_x, sig_y, sig_z, source_Gmatrix, ...
                                tsl_dealiasing_flag, tsl_dealiasing_pad);
                        case 'DL_Stk_3D'
                            error('No reason to use this.');
                        otherwise
                            error('Incorrect potential passed in.');
                    end
                end
            end
        elseif strcmp(body_shape_type, 'sphere')
            switch pot(1:3)
                case 'SL_'
                    pMat = 'SMat';
                case 'SDL'
                    pMat = 'SDMat';
                case 'dSL'
                    pMat = 'SpMat';
                case 'TSL'
                    pMat = 'TSMat';
                case 'DL_'
                    pMat = 'DMat';
                case 'dDL'
                    pMat = 'DpMat';
                case 'TDL'
                    pMat = 'TDMat';
                otherwise
                    error('Unsupported pot for spherical near-eval.');
            end

            % Rescale geometry to unit sphere for spectral near-eval
            % target_pts are already translated (and rotated if rot) into the local frame
            Xtrg = (1/equ_radius) * target_pts;

            % Spherical coordinates on unit sphere frame
            [th, phi, rho] = cart2sph(Xtrg(:,1), Xtrg(:,2), Xtrg(:,3));
            th(th<0) = th(th<0) + 2*pi;
            phi = pi/2 - phi;

            Vloc = V(I_box,:);
            if kerd == 1
                Nrtrg = target_normals;
                Vh_loc = shAna(Vloc);

                if num_ngh>1
                    slf = find(neigh{body_ind}==body_ind); 
                    indv_off = [1:Nb*(slf-1) (Nb*slf+1):Nb*num_ngh].';
                    ind_off  = [1:np*(slf-1) (np*slf+1):np*num_ngh].';

                    Ynear = zeros(Nb*num_ngh, size(Vh_loc,2));

                    % Self term (evaluate at unit radius)
                    if is_modified_laplace
                        Ynear(Nb*(slf-1)+1:Nb*slf, :) = Sh_Mod_Kernel_Eval_off(Vh_loc, pMat, 0, out, 1, [], [], [], lambda);
                    else
                        Ynear(Nb*(slf-1)+1:Nb*slf, :) = Sh_Kernel_Eval_off(Vh_loc, pMat, 0, out, 1, [], [], []);
                    end
                    % Off-diagonal neighbor terms at specified spherical coordinates and normals
                    if ~isempty(ind_off)
                        if is_modified_laplace
                            Ynear(indv_off, :) = Sh_Mod_Kernel_Eval_off( ...
                                Vh_loc, pMat, 0, out, rho(ind_off), phi(ind_off), th(ind_off), Nrtrg(ind_off,:), lambda);
                        else
                            Ynear(indv_off, :) = Sh_Kernel_Eval_off( ...
                                Vh_loc, pMat, 0, out, rho(ind_off), phi(ind_off), th(ind_off), Nrtrg(ind_off,:));
                        end
                    end
                else
                    if is_modified_laplace
                        Ynear = Sh_Mod_Kernel_Eval_off(Vh_loc, pMat, 0, out, 1, [], [], [], lambda);
                    else
                        Ynear = Sh_Kernel_Eval_off(Vh_loc, pMat, 0, out, 1, [], [], []);
                    end
                end
            else
                % Target normals in the local frame, duplicated per component
                Nrtrg = Nor(I_nghv,:);
                if rot
                    Nrtrg = Nrtrg * MRot{body_ind};
                end

                % Densities are assumed in the local body frame
                Vh_loc = VshAna([Vloc(1:3:end,:); Vloc(2:3:end,:); Vloc(3:3:end,:)],'VW');

                if num_ngh>1
                    slf = find(neigh{body_ind}==body_ind); 
                    indv_off = [1:Nb*(slf-1) (Nb*slf+1):Nb*num_ngh].';
                    ind_off  = [1:np*(slf-1) (np*slf+1):np*num_ngh].';

                    Ynear = zeros(Nb*num_ngh, size(Vh_loc,2));

                    % Self term (evaluate at unit radius)
                    Ynear(Nb*(slf-1)+1:Nb*slf, :) = Vsh_Kernel_Eval_off(Vh_loc, pMat, 0, out, 1, [], [], []);
                    % Off-diagonal neighbor terms at specified spherical coordinates and normals
                    if ~isempty(ind_off)
                        Ynear(indv_off, :) = Vsh_Kernel_Eval_off(Vh_loc, pMat, 0, out, rho(ind_off), phi(ind_off), th(ind_off), Nrtrg(indv_off,:));
                    end
                else
                    Ynear = Vsh_Kernel_Eval_off(Vh_loc, pMat, 0, out, 1, [], [], []);
                end

                if size(Ynear,2) > 1
                    Ynear = sum(Ynear, 2);
                end
            end

            % Scale from unit sphere back to radius r
            if strncmp(pot,'SL_',3)
                rda = equ_radius; % SLP scales like r
            elseif strncmp(pot,'dDL',3)
                rda = 1/equ_radius; % dDL scales like 1/r
            else
                rda = 1;
            end
            Ynear = rda * Ynear;
        else
            error('Invalid body shape type; should be "sphere" or "prolate" or "oblate".');
        end

        % Post-processing of Ynear
        if rot && kerd==3
            % Rotate back
            Ynear = [Ynear(1:3:end,:), Ynear(2:3:end,:), Ynear(3:3:end,:)];  
            Ynear = Ynear*MRot{body_ind};       

            % Go back to interleaved format
            Ynear = reshape(reshape(Ynear,[],3).',[],1);
        end

        %% Far-interaction
        % This should not change too much from VSh_MatVec_RB2.m.
        far_idx_v = true(N,1);
        far_idx_v(I_nghv) = false; 

        % Local params copy for Kernel_Eval
        parnear = params;

        if dense
            if nortrg
                parnear.nor = Nor(far_idx_v,:);
            else
                % DL/SDL: source normals per source block
                parnear.nor = Nor(I_box,:);
            end
            parnear.W2 = W2(I_box);
            if isfield(parnear,'ci')
                parnear.ci = repmat((1:kerd)',sum(far_idx_v)/kerd,1);
            end
            if isfield(parnear,'cj')
                parnear.cj = params.cj(I_box);
            end

            % Add near spectral contribution
            Y(I_nghv,:) = Y(I_nghv,:) + Ynear;

            % Add far-away interactions via direct Kernel_Eval
            if any(far_idx_v)
                ke_eval = Kernel_Eval(Xv(far_idx_v,:), Xv(I_box,:), parnear) * V(I_box,:);
                Y(far_idx_v,:) = Y(far_idx_v,:) + ke_eval;
            end
        else
            if nortrg
               parnear.nor = Nor(I_nghv,:);  
            else
               parnear.nor = Nor(I_box,:);
            end
            parnear.W2 = W2(I_box); 
            if isfield(parnear,'ci')
                parnear.ci = repmat((1:kerd)',sum(~far_idx_v)/3,1);
            end
            if isfield(parnear,'cj')
                parnear.cj = params.cj(I_box);
            end
            parnear.a = 0; 
            
            % Subtract neighbor Kernel_Eval; keep near spectral contribution
            Y(I_nghv,:) = (Y(I_nghv,:) + Ynear) - Kernel_Eval(Xv(I_nghv,:), Xv(I_box,:), parnear) * V(I_box,:);
        end
    end

    % If using FMM flag (dense==false), add far interactions via FMM
    if ~dense
        Nr = Nor(1:kerd:end,:);
        W = W2.';
        Y = Y + SSph_FMM_Eval(V, W, kerd, pot, X, X, Nr);
    end

    if out
        ct = 0.5;
    else
        ct = -0.5;
    end

    % Add jump-relation (i.e. params.a)
    if ~strcmp(pot(1:3),'SL_') && ~strcmp(pot(1:3),'dDL')
        if isfield(params,'a')
            if strcmp(pot(2:3),'SL') 
                if isfield(params,'eta') && strcmp(pot(1:3),'dSL') 
                    Y = params.eta*Y + (params.a+ct)*V; 
                else
                    Y = Y + (params.a+ct)*V; 
                end
            else
                Y = Y + (params.a-ct)*V; 
            end
        end
    end

    % Add nullspace correction term (block-diagonal)
    if ~isempty(L)
        Y = Y + L*V;
    end

elseif isempty(V)
    fprintf('\nGoing matrix free...\n')
    % Matrix-free
    Y = @(V) SSph_MatVec(V,L,params,modlap_opts);
else
    error('Invalid input passed into SSph_MatVec. Supposed to be a string "Mat" or empty or a numeric matrix.');
end

end %% END SSph_MatVec

function Ynear = LOCAL_spheroid_near_matrix( ...
    pot, u0, a, oblate, target_pts, target_normals, np, iprm, source_Gmatrix, lambda, modlap_opts, ...
    tsl_dealiasing_flag, tsl_dealiasing_pad)
    %{
    Builds the dense block for near-field interaction from one source spheroid by applying
    the near matvec to unit-basis densities.

    Inputs
    pot - (string) potential name
    u0 - (double)
    a - (double)
    oblate - (bool)
    target_pts - (double ntrg x 3) target points in the local frame
        if empty, implies self-evaluation on the source surface
    target_normals - (double ntrg x 3) target normals in the local frame
        required for traction and normal-derivative evaluations
    np - (int) number of discretization points
    iprm - (int 1 x 3*np) permutation from blocked to interleaved

    Output
    Ynear - (double 3*ntrg x 3*np) dense near-interaction block in interleaved format
    %}

    I = eye(np);
    Z = zeros(np);

    % Build parameters once and reuse across the three basis evaluations.
    params_i = SpheroidalParameters();
    params_i.sigma = I; % set p for geometry generation
    params_i.u0 = u0;
    params_i.a = a;
    params_i.oblate = oblate;
    params_i.centers = [0 0 0];
    params_i.thetas = 0;
    params_i.phis = 0;
    params_i.Rmat = eye(3);

    if isempty(target_pts)
        X_eval = [];
        Nu_eval = [];
    else
        X_eval = reshape(target_pts, size(target_pts,1), 3, 1);
        Nu_eval = reshape(target_normals, size(target_normals,1), 3, 1);
    end

    switch pot
        case {'SL_L_3D', 'dSL_L_3D', 'DL_L_3D', 'dDL_L_3D', ...
                'SL_LMOD_3D', 'dSL_LMOD_3D', 'DL_LMOD_3D', 'dDL_LMOD_3D'}
            Ynear = LOCAL_eval_scalar_spheroidal_potential(pot, params_i, X_eval, Nu_eval, I, source_Gmatrix, lambda, modlap_opts);
        case 'SL_Stk_3D'
            Yx = LOCAL_eval_l2stk_slp(params_i, X_eval, I, Z, Z, source_Gmatrix);
            Yy = LOCAL_eval_l2stk_slp(params_i, X_eval, Z, I, Z, source_Gmatrix);
            Yz = LOCAL_eval_l2stk_slp(params_i, X_eval, Z, Z, I, source_Gmatrix);
        case 'DL_Stk_3D'
            error('No reason to use this for the mobility solver.');
        case 'TSL_Stk_3D'
            Yx = LOCAL_eval_l2stk(params_i, X_eval, Nu_eval, I, Z, Z, source_Gmatrix, tsl_dealiasing_flag, tsl_dealiasing_pad);
            Yy = LOCAL_eval_l2stk(params_i, X_eval, Nu_eval, Z, I, Z, source_Gmatrix, tsl_dealiasing_flag, tsl_dealiasing_pad);
            Yz = LOCAL_eval_l2stk(params_i, X_eval, Nu_eval, Z, Z, I, source_Gmatrix, tsl_dealiasing_flag, tsl_dealiasing_pad);
        otherwise
            error('Incorrect potential passed in.');
    end

    if isempty(target_pts)
        ntrg = np;
    else
        ntrg = size(target_pts,1);
    end

    if LOCAL_is_scalar_laplace_potential(pot)
        Ynear = reshape(Ynear, ntrg, np);
        Ynear = Ynear(:, iprm);
        return;
    end

    Yx = reshape(Yx, 3*ntrg, np);
    Yy = reshape(Yy, 3*ntrg, np);
    Yz = reshape(Yz, 3*ntrg, np);
    Ynear = [Yx, Yy, Yz];

    % Go from block format to interleaved format.
    Ynear = Ynear(:, iprm);
end

function eval = LOCAL_eval_l2stk_slp(params_i, X_eval, sigma_x, sigma_y, sigma_z, source_Gmatrix)
    [vx, vy, vz] = L2StkSLPOptimized(X_eval, params_i, sigma_x, sigma_y, sigma_z, source_Gmatrix);

    % Interleave result
    eval = reshape([vx(:), vy(:), vz(:)].', [], 1);
end

function eval = LOCAL_eval_l2stk(params_i, X_eval, Nu_eval, sigma_x, sigma_y, sigma_z, source_Gmatrix, tsl_dealiasing_flag, tsl_dealiasing_pad)
    [vx, vy, vz] = L2StkTLPOptimized( ...
        X_eval, Nu_eval, params_i, sigma_x, sigma_y, sigma_z, source_Gmatrix, false, tsl_dealiasing_flag, tsl_dealiasing_pad);

    % Interleave result
    eval = reshape([vx(:), vy(:), vz(:)].', [], 1);
end

function eval = LOCAL_eval_scalar_spheroidal_potential(pot, params_i, X_eval, Nu_eval, sigma, source_Gmatrix, lambda, modlap_opts)
    switch pot
        case 'SL_L_3D'
            eval = LOCAL_eval_lslp(params_i, X_eval, sigma, source_Gmatrix);
        case 'dSL_L_3D'
            eval = LOCAL_eval_ldslp(params_i, X_eval, Nu_eval, sigma, source_Gmatrix);
        case 'SL_LMOD_3D'
            eval = LOCAL_eval_lmod_slp(params_i, X_eval, sigma, lambda, modlap_opts);
        case 'dSL_LMOD_3D'
            eval = LOCAL_eval_lmod_sp(params_i, X_eval, sigma, lambda, modlap_opts);
        case 'DL_LMOD_3D'
            eval = LOCAL_eval_lmod_dlp(params_i, X_eval, sigma, lambda, modlap_opts);
        case 'dDL_LMOD_3D'
            eval = LOCAL_eval_lmod_dp(params_i, X_eval, sigma, lambda, modlap_opts);
        otherwise
            error('Incorrect scalar Laplace potential passed in.');
    end
end

function eval = LOCAL_eval_lmod_slp(params_i, X_eval, sigma, lambda, modlap_opts)
    if isempty(lambda)
        error('Missing lambda for SL_LMOD_3D.');
    end

    params_eval = copy(params_i);
    params_eval.sigma = sigma;
    params_eval.get_shc();

    if isempty(X_eval)
        modSL = spheroidalModifiedSLP(params_eval, lambda, [], modlap_opts);
    else
        X_trg = reshape(X_eval, size(X_eval,1), 3, 1);
        modSL = spheroidalModifiedSLP(params_eval, lambda, X_trg, modlap_opts);
    end

    eval = reshape(modSL, size(modSL,1), size(modSL,2));
end

function eval = LOCAL_eval_lmod_dlp(params_i, X_eval, sigma, lambda, modlap_opts)
    if isempty(lambda)
        error('Missing lambda for DL_LMOD_3D.');
    end

    params_eval = copy(params_i);
    params_eval.sigma = sigma;
    params_eval.get_shc();

    if isempty(X_eval)
        modDL = spheroidalModifiedDLP(params_eval, lambda, [], modlap_opts);
    else
        X_trg = reshape(X_eval, size(X_eval,1), 3, 1);
        modDL = spheroidalModifiedDLP(params_eval, lambda, X_trg, modlap_opts);
    end

    eval = reshape(modDL, size(modDL,1), size(modDL,2));
end

function eval = LOCAL_eval_lmod_sp(params_i, X_eval, sigma, lambda, modlap_opts)
    if isempty(lambda)
        error('Missing lambda for dSL_LMOD_3D.');
    end

    params_eval = copy(params_i);
    params_eval.sigma = sigma;
    params_eval.get_shc();

    if isempty(X_eval)
        modSP = spheroidalModifiedSP(params_eval, lambda, [], modlap_opts);
    else
        X_trg = reshape(X_eval, size(X_eval,1), 3, 1);
        modSP = spheroidalModifiedSP(params_eval, lambda, X_trg, modlap_opts);
    end

    eval = reshape(modSP, size(modSP,1), size(modSP,2));
end

function eval = LOCAL_eval_lmod_dp(params_i, X_eval, sigma, lambda, modlap_opts)
    if isempty(lambda)
        error('Missing lambda for dDL_LMOD_3D.');
    end

    params_eval = copy(params_i);
    params_eval.sigma = sigma;
    params_eval.get_shc();

    if isempty(X_eval)
        modDP = spheroidalModifiedDP(params_eval, lambda, [], modlap_opts);
    else
        X_trg = reshape(X_eval, size(X_eval,1), 3, 1);
        modDP = spheroidalModifiedDP(params_eval, lambda, X_trg, modlap_opts);
    end

    eval = reshape(modDP, size(modDP,1), size(modDP,2));
end

function tf = LOCAL_is_scalar_laplace_potential(pot)
    tf = any(strcmp(pot, { ...
        'SL_L_3D', 'dSL_L_3D', 'DL_L_3D', 'dDL_L_3D', ...
        'SL_LMOD_3D', 'dSL_LMOD_3D', 'DL_LMOD_3D', 'dDL_LMOD_3D' ...
    }));
end

function eval = LOCAL_eval_lslp(params_i, X_eval, sigma, source_Gmatrix)
    shc = shAna(sigma);
    if isempty(source_Gmatrix)
        Gshc = shc;
    else
        Gshc = source_Gmatrix \ shc;
    end
    Gshc = reshape(Gshc, size(Gshc,1), size(Gshc,2), 1);

    if isempty(X_eval)
        X_trg = [];
    else
        X_trg = reshape(X_eval, size(X_eval,1), 3, 1);
    end

    SL = spheroidalSLOptimized(params_i.p, params_i.u0, params_i.a, params_i.oblate, Gshc, params_i.isReal, X_trg);
    eval = reshape(SL, size(SL,1), size(SL,2));
end

function eval = LOCAL_eval_ldslp(params_i, X_eval, Nu_eval, sigma, source_Gmatrix)
    shc = shAna(sigma);
    if isempty(source_Gmatrix)
        Gshc = shc;
    else
        Gshc = source_Gmatrix \ shc;
    end
    Gshc = reshape(Gshc, size(Gshc,1), size(Gshc,2), 1);

    if isempty(X_eval)
        X_trg = [];
        Nu = params_i.get_Norm();
    else
        X_trg = reshape(X_eval, size(X_eval,1), 3, 1);
        Nu = Nu_eval;
    end

    if isempty(Nu)
        error('Target normals required for dSL_L_3D evaluation.');
    end

    nu_x = reshape(Nu, size(Nu,1), 3, 1);
    nu_y = nu_x;
    nu_z = nu_x;

    [SPx, ~, ~] = spheroidalSPOptimized(params_i.p, params_i.u0, params_i.a, params_i.oblate, Gshc, params_i.isReal, nu_x, nu_y, nu_z, X_trg);
    eval = reshape(SPx, size(SPx,1), size(SPx,2));
end

function [u0, a, oblate] = LOCAL_calculate_u0_a(equ_radius, polar_radius, shape_type)
    %{
    TODO: Don't do this.
    %}
    [u0, a] = calculate_u0_and_a_from_radii(shape_type, equ_radius, polar_radius);
    oblate = strcmp(shape_type,'oblate');
end

function [tsl_dealiasing_flag, tsl_dealiasing_pad] = LOCAL_get_tsl_dealiasing_options(params)
    if isfield(params, 'tsl_dealiasing')
        tsl_dealiasing_flag = logical(params.tsl_dealiasing);
    else
        tsl_dealiasing_flag = false;
    end

    if isfield(params, 'tsl_dealiasing_pad')
        tsl_dealiasing_pad = params.tsl_dealiasing_pad;
    else
        tsl_dealiasing_pad = 4;
    end
end


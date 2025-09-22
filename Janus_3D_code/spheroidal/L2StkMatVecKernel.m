function M = L2StkMatVecKernel(pars, pot, p, nu_eval)
    %{
        Finds the matrix representation of the matvec for a given Stokes
        potential.
    %}

    if isempty(pars.u0)
        error("No surface parameter u_0 given")
    end
    if isempty(pars.a)
        error("No scale (a) given")
    end
    if isempty(pars.matvec_eta)
        error("No matvec_eta provided.");
    end

    if p < 2
        error('KernelEval (smooth quadrature) requires at least p=2.')
    end

    u0=pars.u0;
    a=pars.a;
    centers = pars.centers;
    thetas = pars.thetas;
    phis = pars.phis;
    if_oblate = pars.oblate;

    ns = size(centers,1);
    np = 2*p*(p+1);

    if length(u0)==1
        u0=u0*ones(1,ns); 
    end
    if length(a)==1
        a=a*ones(1,ns); 
    end

    % Do self evaluation of all particles (matrices on main block diagonal)
    % ------------------------------------------------------------
    IDparams = copy(pars);
    IDparams.sigma = repmat(eye(np),1,1,ns);

    if strcmp(pot, 'SLP') % Not yet implemented!
        M = SLPmatrix(IDparams, [], np, ns);
    elseif strcmp(pot, 'DLP')
        M = DLPmatrix_dev(IDparams, [], np, ns);
        M = blkdiag(M{:});
    elseif strcmp(pot, 'TLP')
        M = TSLmatrix(IDparams, p, ns);
    else
        error("Invalid potential given: should be 'SLP', 'DLP', or 'TLP'.");
    end

    separation = pars.separate_spheroids();

    if ns > 1
        fprintf("Generating off-diagonal blocks for multiple bodies (%d bodies)...\n", ns);
        for i=1:ns % Source particle i
            % This is necessary due to how the spheroidal LPs are setup.
            % Otherwise, we would have to do something similar to
            % spheroidalMatVecKernel's code.
            %
            % The idea is to make a SpheroidalParameters that only has
            % knowledge of the source particle. The other particles are
            % then just arbitrary information for the points and normals.
            params_i = IDparams.copy();
            params_i.u0 = IDparams.u0(i);
            params_i.a = IDparams.a(i);
            params_i.oblate = IDparams.oblate(i);
            params_i.sigma = eye(np);

            % Recall that we shift the entire coordinate system so that the
            % spheroid is upright and is centered at the origin. Thus, the
            % only complication are the target points.
            params_i.centers = [0 0 0];
            params_i.thetas = 0;
            params_i.phis = 0;
            params_i.Rmat = eye(3);

            % Get all coordinates of other spheroids with respect to
            % spheroid i
            [X_all_targets, X_self] = IDparams.get_X(i);
            [Nu_all_targets, Nu_self] = IDparams.get_Norm_rot(p,i);

            spheroid_index = 0;
            for j=1:ns % Target particle j
                if i==j % Self-interaction already handled above.
                    continue;
                end

                spheroid_index = spheroid_index + 1;
                row_indices = (spheroid_index - 1)*np+1 : spheroid_index*np;

                X_target = X_all_targets(row_indices,:);
                Nu_target = Nu_all_targets(row_indices,:);
                M_ji = zeros(3*np, 3*np); % Interaction of j on i

                is_smooth = separation(i, j);
                if ~is_smooth % Near
                    if strcmp(pot, 'SLP')
                        error("not implemented.");
                        M_cell = SLPmatrix(IDparams, X_target, np, 1);
                        M_ji = M_cell{1};
                    elseif strcmp(pot, 'DLP')
                        M_ji = DLP_matrix_at_target(params_i, X_target, np);
                    elseif strcmp(pot, 'TLP')
                        M_ji = TSLmatrix_at_target(params_i, X_target, Nu_target, p);
                    end
                else % Far
                    fprintf("Handling smooth points in L2StkMatVecKernel...\n");
                    %%% Setup Kernel_Eval parameters
                    Sns=SurfaceSph(X_self);
    
                    if strcmp(pot, 'SLP')
                        KEpot = 'SL_Stk_3D';
                    elseif strcmp(pot, 'DLP')
                        KEpot = 'DL_Stk_3D';
                    elseif strcmp(pot, 'TLP')
                        KEpot = 'TSL_Stk_3D';
                    else
                        error("Invalid potential given.");
                    end
        
                    KEparams = Kernel_Eval_parameters(KEpot,0,1,1,1,1e-8,2,400,1);
                    KEparams.dim = 3; KEparams.mu=1;
                    [~, gwt]=g_grid(p+1);
                    wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
                    wt = wt(:);
                    Wns = Sns.geoProp.W .* wt;
                    Wv = repmat(Wns,1,3)'; Wv=Wv(:);
                    Xv = reshape(repmat(X_self,1,3)',3,[])';
                    Xtrg_ii = reshape(repmat(X_target,1,3)',3,[])';
                    KEparams.X = Xv;
                    KEparams.W2 = Wv.';
                    KEparams.cj=repmat((1:3)',size(X_self,1),1);
                    KEparams.ci=repmat((1:3)',size(X_target,1),1);
    
                    if strcmp(pot, 'TLP')
                        KEparams.nor = reshape(repmat(Nu_target,1,3)',3,[])';
                    elseif strcmp(pot, 'DLP')
                        KEparams.nor = reshape(repmat(Nu_self,1,3)',3,[])';
                    end
        
                    M_ji = Kernel_Eval(Xtrg_ii,Xv,KEparams);

                    prm = zeros(1,3*np); 
                    prm(1:np) = 1:3:3*np; prm(np+1:2*np) = 2:3:3*np; prm(2*np+1:3*np)=3:3:3*np;
                    M_ji = M_ji(prm, prm);
                end

                % Rotation matrix: move from frame i to frame j
                % Note that it should not be kron(eye(np), R) as this
                % assumes a different format for the vectors (i.e. we want
                % the x, y, and z's separately blocked together).
                R = kron(IDparams.Rmat(:,:,j)' * IDparams.Rmat(:,:,i), eye(np));

                % Rows are targets, columns are source
                % Thus, column i represents the block from the perspective
                % of spheroid i. However, we don't want this: we want row i
                % to correspond to spheroid i. This is why we multiply by
                % R.
                target_block = (j-1)*(3*np) + (1:3*np);
                source_block = (i-1)*(3*np) + (1:3*np);
                M(target_block, source_block) = R*M_ji;
            end % End for target spheroid j
        end % End for source spheroid i
    end
end

%% DEVELOPMENT
function M = DLP_matrix_at_target(params_i, X_target, np)
    nu_x_all = repmat([1,0,0], np, 1);
    nu_y_all = repmat([0,1,0], np, 1);
    nu_z_all = repmat([0,0,1], np, 1);

    [SP_x, SP_y, SP_z] = spheroidalSP(params_i, X_target, nu_x_all, nu_y_all, nu_z_all);
    [DP_x, DP_y, DP_z] = spheroidalDP(params_i, X_target, nu_x_all, nu_y_all, nu_z_all);

    X_src = params_i.get_X - params_i.centers;
    N_src = params_i.get_Norm(params_i.p, 1);
    X_trg = X_target;

    % Extract DP sum term
    % [ x*DPx y*DPx z*DPx ] [ sigma_x ]
    % [ x*DPy y*DPy z*DPy ] [ sigma_y ]
    % [ x*DPz y*DPz z*DPz ] [ sigma_z ]
    M1_11 = diag(X_trg(:,1)) * DP_x; M1_12 = diag(X_trg(:,2)) * DP_x; M1_13 = diag(X_trg(:,3)) * DP_x;
    M1_21 = diag(X_trg(:,1)) * DP_y; M1_22 = diag(X_trg(:,2)) * DP_y; M1_23 = diag(X_trg(:,3)) * DP_y;
    M1_31 = diag(X_trg(:,1)) * DP_z; M1_32 = diag(X_trg(:,2)) * DP_z; M1_33 = diag(X_trg(:,3)) * DP_z;
    M1 = [
        M1_11, M1_12, M1_13;
        M1_21, M1_22, M1_23;
        M1_31, M1_32, M1_33
    ];

    % Extract individual DP term; i.e. DP[y \cdot \sigma]
    % [ DPx ]                        [ sigma_x ]
    % [ DPy ][ Xsrc_x Xsrc_y Xsrc_z ][ sigma_y ]
    % [ DPz ]                        [ sigma_z ]
    M2 = [DP_x; DP_y; DP_z] * [diag(X_src(:,1)), diag(X_src(:,2)), diag(X_src(:,3))];

    % Extract SP sum term
    % [ SPx ][ Nx Nx Nx ][ sigma_x ]
    % [ SPy ][ Ny Ny Ny ][ sigma_y ]
    % [ SPz ][ Nz Nz Nz ][ sigma_z ]
    M3_11 = SP_x * diag(N_src(:,1)); M3_12 = SP_y * diag(N_src(:,1)); M3_13 = SP_z * diag(N_src(:,1));
    M3_21 = SP_x * diag(N_src(:,2)); M3_22 = SP_y * diag(N_src(:,2)); M3_23 = SP_z * diag(N_src(:,2));
    M3_31 = SP_x * diag(N_src(:,3)); M3_32 = SP_y * diag(N_src(:,3)); M3_33 = SP_z * diag(N_src(:,3));
    M3 = [
        M3_11, M3_12, M3_13;
        M3_21, M3_22, M3_23;
        M3_31, M3_32, M3_33
    ];

    M = -M1 + M2 - M3;
end

function M_cells = DLPmatrix_dev(IDparams, X_spectral, np, ns)
    is_self_interaction = isempty(X_spectral);

    M_cells = cell(1, ns);
    
    nu_x_cells = cell(1, ns);
    nu_y_cells = cell(1, ns);
    nu_z_cells = cell(1, ns);

    for i = 1:ns
        if is_self_interaction
            num_trg = np;
        else
            num_trg = size(X_spectral{i}, 1);
        end
        
        if num_trg > 0
            nu_x_cells{i} = repmat([1,0,0], num_trg, 1);
            nu_y_cells{i} = repmat([0,1,0], num_trg, 1);
            nu_z_cells{i} = repmat([0,0,1], num_trg, 1);
        end
    end

    % Cell i corresponds to particle i acting on target points in X_spectral 
    [SP_x_cells, SP_y_cells, SP_z_cells] = spheroidalSP(IDparams, X_spectral, nu_x_cells, nu_y_cells, nu_z_cells);
    [DP_x_cells, DP_y_cells, DP_z_cells] = spheroidalDP(IDparams, X_spectral, nu_x_cells, nu_y_cells, nu_z_cells);

    for i=1:ns
        if ~is_self_interaction && isempty(X_spectral{i})
            continue;
        end

        [~, X_src_i] = IDparams.get_X(i);
        N_src_i = IDparams.get_Norm(IDparams.p, i); %maybe rotation is needed?
        
        % Get target points
        if is_self_interaction
            X_trg_i = X_src_i;
        else
            X_trg_i = X_spectral{i};
        end
        

        % Extract SP terms
        SP_x_i = SP_x_cells{i};
        SP_y_i = SP_y_cells{i};
        SP_z_i = SP_z_cells{i};

        % Extract DP terms
        DP_x_i = DP_x_cells{i};
        DP_y_i = DP_y_cells{i};
        DP_z_i = DP_z_cells{i};

        % Extract DP sum term
        % [ x*DPx y*DPx z*DPx ] [ sigma_x ]
        % [ x*DPy y*DPy z*DPy ] [ sigma_y ]
        % [ x*DPz y*DPz z*DPz ] [ sigma_z ]
        M1_11 = diag(X_trg_i(:,1)) * DP_x_i; M1_12 = diag(X_trg_i(:,2)) * DP_x_i; M1_13 = diag(X_trg_i(:,3)) * DP_x_i;
        M1_21 = diag(X_trg_i(:,1)) * DP_y_i; M1_22 = diag(X_trg_i(:,2)) * DP_y_i; M1_23 = diag(X_trg_i(:,3)) * DP_y_i;
        M1_31 = diag(X_trg_i(:,1)) * DP_z_i; M1_32 = diag(X_trg_i(:,2)) * DP_z_i; M1_33 = diag(X_trg_i(:,3)) * DP_z_i;
        M1 = [
            M1_11, M1_12, M1_13;
            M1_21, M1_22, M1_23;
            M1_31, M1_32, M1_33
        ];

        % Extract individual DP term; i.e. DP[y \cdot \sigma]
        % [ DPx ]                        [ sigma_x ]
        % [ DPy ][ Xsrc_x Xsrc_y Xsrc_z ][ sigma_y ]
        % [ DPz ]                        [ sigma_z ]
        M2 = [DP_x_i; DP_y_i; DP_z_i] * [diag(X_src_i(:,1)), diag(X_src_i(:,2)), diag(X_src_i(:,3))];

        % Extract SP sum term
        % [ SPx ][ Nx Nx Nx ][ sigma_x ]
        % [ SPy ][ Ny Ny Ny ][ sigma_y ]
        % [ SPz ][ Nz Nz Nz ][ sigma_z ]
        M3_11 = SP_x_i * diag(N_src_i(:,1)); M3_12 = SP_y_i * diag(N_src_i(:,1)); M3_13 = SP_z_i * diag(N_src_i(:,1));
        M3_21 = SP_x_i * diag(N_src_i(:,2)); M3_22 = SP_y_i * diag(N_src_i(:,2)); M3_23 = SP_z_i * diag(N_src_i(:,2));
        M3_31 = SP_x_i * diag(N_src_i(:,3)); M3_32 = SP_y_i * diag(N_src_i(:,3)); M3_33 = SP_z_i * diag(N_src_i(:,3));
        M3 = [
            M3_11, M3_12, M3_13;
            M3_21, M3_22, M3_23;
            M3_31, M3_32, M3_33
        ];

        M_cells{i} = -M1 + M2 - M3;
    end
end
function M_cells = SLPmatrix_dev(IDparams, X_spectral, np, ns)
    is_self_interaction = isempty(X_spectral);

    M_cells = cell(1, ns);
    
    nu_x_cells = cell(1, ns);
    nu_y_cells = cell(1, ns);
    nu_z_cells = cell(1, ns);

    for i = 1:ns
        if is_self_interaction
            num_trg = np;
        else
            num_trg = size(X_spectral{i}, 1);
        end
        
        if num_trg > 0
            nu_x_cells{i} = repmat([1,0,0], num_trg, 1);
            nu_y_cells{i} = repmat([0,1,0], num_trg, 1);
            nu_z_cells{i} = repmat([0,0,1], num_trg, 1);
        end
    end
    
    % Get all SP matrices. Each output is np x np x ns.
    SLM = spheroidalSL(IDparams);
    [SP_x_cells, SP_y_cells, SP_z_cells] = spheroidalSP(IDparams, X_spectral, nu_x_cells, nu_y_cells, nu_z_cells);

    for i=1:ns
        if ~is_self_interaction && isempty(X_spectral{i})
            continue;
        end

        if is_self_interaction
            [~, X_src_i] = IDparams.get_X(i);
            X_trg_i = X_src_i;
        else
            X_trg_i = X_spectral{i};
        end

        % Extract SL term
        % M1 for each spheroid: 3np x 3np
        % [ SL 0  0 ][ sigma_x ] = [ SL[sigma_x] ]
        % [ 0  SL 0 ][ sigma_y ]   [ SL[sigma_y] ]
        % [ 0  0  SL][ sigma_z ]   [ SL[sigma_z] ] 
        M1 = kron(eye(3), SLM(:,:,i));

        % Extract SP terms
        SP_x_i = SP_x_cells{i};
        SP_y_i = SP_y_cells{i};
        SP_z_i = SP_z_cells{i};

        % Surface SP but with arbitrary derivative vector, 
        % dSx: np x np; 
        % M2 for each spheroid: 3np x 3np
        % [ x*dSx y*dSx z*dSx ][ sigma_x ]
        % [ x*dSy y*dSy z*dSy ][ sigma_y ]
        % [ x*dSz y*dSz z*dSz ][ sigma_z ]
        Dx = diag(X_trg_i(:,1));
        Dy = diag(X_trg_i(:,2));
        Dz = diag(X_trg_i(:,3));
        M2 = [
            Dx*SP_x_i, Dy*SP_x_i, Dz*SP_x_i;
            Dx*SP_y_i, Dy*SP_y_i, Dz*SP_y_i;
            Dx*SP_z_i, Dy*SP_z_i, Dz*SP_z_i
        ];

        % SPM[Xsrc dot sigma]
        % dSx: np x np;
        % Xsrc_x: np x np;
        % [ dSx ]                        [ sigma_x ]
        % [ dSy ][ Xsrc_x Xsrc_y Xsrc_z ][ sigma_y ]
        % [ dSz ]                        [ sigma_z ]
        M3 = [SP_x_i; SP_y_i; SP_z_i] * [Dx, Dy, Dz];

        % Collect all of the results
        M_cells{i} = M1 - M2 + M3;
    end
end

function M = TSLmatrix_at_target(params_i, X_target, Nu_target, p)
    np = 2*p*(p+1);
    nu_x_all = repmat([1,0,0], np, 1);
    nu_y_all = repmat([0,1,0], np, 1);
    nu_z_all = repmat([0,0,1], np, 1);

    [SP_x, SP_y, SP_z] = spheroidalSP(params_i, X_target, nu_x_all, nu_y_all, nu_z_all);

    %% Build spheroidalgraddiv operator
    I = eye(np);
    Z = zeros(np);

    % Essentially, reconstruct the Hessian.
    [D1_U, D1_V, D1_PHI] = spheroidalgraddivSL(params_i, I, Z, Z, X_target);
    [D2_U, D2_V, D2_PHI] = spheroidalgraddivSL(params_i, Z, I, Z, X_target);
    [D3_U, D3_V, D3_PHI] = spheroidalgraddivSL(params_i, Z, Z, I, X_target);

    S_target = cart2spheroidal(X_target, params_i.a, params_i.oblate);
    u = S_target(:,1); v = S_target(:,2); phi = S_target(:,3);
    if params_i.oblate
        [DXX, DYX, DZX] = convert_from_oblate_basis(D1_U, D1_V, D1_PHI, u, v, phi);
        [DXY, DYY, DZY] = convert_from_oblate_basis(D2_U, D2_V, D2_PHI, u, v, phi);
        [DXZ, DYZ, DZZ] = convert_from_oblate_basis(D3_U, D3_V, D3_PHI, u, v, phi);
    else
        [DXX, DYX, DZX] = convert_from_prolate_basis(D1_U, D1_V, D1_PHI, u, v, phi);
        [DXY, DYY, DZY] = convert_from_prolate_basis(D2_U, D2_V, D2_PHI, u, v, phi);
        [DXZ, DYZ, DZZ] = convert_from_prolate_basis(D3_U, D3_V, D3_PHI, u, v, phi);
    end
    
    %{
        Note that the gradient of the divergence should be represented as
        follows. Let (F_x, F_y, F_z) denote the vector field of interest.
        [(d^2x F_x) + (dxdy F_y) + (dxdz F_z)]
        [(dydx F_x) + (d^2y F_y) + (dydz F_z)]
        [(dzdx F_x) + (dzdy F_y) + (d^2z F_z)]

        Alternatively,
        [d^2x dxdy dxdz][F_x]
        [dydx d^2y dydz][F_y]
        [dzdx dzdy d^2z][F_z]

        The matrix above is ddS_op.
    %}
    
    % TODO: for memory efficiency, one can use commutativity of second-order
    % derivatives (i.e. we should have DXY=DYX, DYZ=DZY, and DXZ = DZX).
    % However, not sure how this relationship holds in near-singular
    % evaluation.
    ddS_x = [DXX, DXY, DXZ];
    ddS_y = [DYX, DYY, DYZ];
    ddS_z = [DZX, DZY, DZZ];
    ddS_op = [ddS_x; ddS_y; ddS_z];

    %% Build TLP matvec kernel
    % Grab source points
    X_src = params_i.get_X - params_i.centers;
    Xtarget_dot_N = repmat(dot(X_target, Nu_target, 2), 3, 1); % 3np x 1 vector.

    %%%
    %%% Build matrix
    %%%
    % Extract SP sum term
    M1_diag_term = diag(Nu_target(:,1))*SP_x + diag(Nu_target(:,2))*SP_y + diag(Nu_target(:,3))*SP_z;
    M1 = kron(eye(3), M1_diag_term);

    % Extract ddS sum term
    % The idea is that we multiply on the left by the normal vector term, on the right
    % by the y_j term, and the ddS operator is represented by [ddS_x; ddS_y; ddS_z].
    M2_1 = kron(eye(3), diag(Nu_target(:,1))) * ddS_op * kron(eye(3), diag(X_src(:,1)));
    M2_2 = kron(eye(3), diag(Nu_target(:,2))) * ddS_op * kron(eye(3), diag(X_src(:,2)));
    M2_3 = kron(eye(3), diag(Nu_target(:,3))) * ddS_op * kron(eye(3), diag(X_src(:,3)));
    M2 = M2_1 + M2_2 + M2_3;

    % Extract ddS term
    M3 = diag(Xtarget_dot_N) * ddS_op;

    M = M1 + M2 - M3;
end

%% SELF TO SELF
function M = SLPmatrix(IDparams, X_spectral, np, ns)
    SLM = spheroidalSL(IDparams);
    M_cells = cell(1, ns);
    nu_x_all = repmat([1,0,0], np, 1, ns);
    nu_y_all = repmat([0,1,0], np, 1, ns);
    nu_z_all = repmat([0,0,1], np, 1, ns);
    
    % Get all SP matrices. Each output is np x np x ns.
    [SP_x_cells, SP_y_cells, SP_z_cells] = spheroidalSP(IDparams, X_spectral, nu_x_all, nu_y_all, nu_z_all);

    % A crutch to handle the ns=1 case...
    if ns == 1
        SP_x_cells = {SP_x_cells}; SP_y_cells = {SP_y_cells}; SP_z_cells = {SP_z_cells};
    end

    for i=1:ns
        [~, X_src_i] = IDparams.get_X(i);

        % Extract SL term
        % M1 for each spheroid: 3np x 3np
        % [ SL 0  0 ][ sigma_x ] = [ SL[sigma_x] ]
        % [ 0  SL 0 ][ sigma_y ]   [ SL[sigma_y] ]
        % [ 0  0  SL][ sigma_z ]   [ SL[sigma_z] ] 
        M1 = kron(eye(3), SLM);

        % Extract SP terms
        SP_x_i = SP_x_cells{i};
        SP_y_i = SP_y_cells{i};
        SP_z_i = SP_z_cells{i};

        % Surface SP but with arbitrary derivative vector, 
        % dSx: np x np; 
        % M2 for each spheroid: 3np x 3np
        % [ x*dSx y*dSx z*dSx ][ sigma_x ]
        % [ x*dSy y*dSy z*dSy ][ sigma_y ]
        % [ x*dSz y*dSz z*dSz ][ sigma_z ]
        Dx = diag(X_src_i(:,1));
        Dy = diag(X_src_i(:,2));
        Dz = diag(X_src_i(:,3));
        M2 = [
            Dx*SP_x_i, Dy*SP_x_i, Dz*SP_x_i;
            Dx*SP_y_i, Dy*SP_y_i, Dz*SP_y_i;
            Dx*SP_z_i, Dy*SP_z_i, Dz*SP_z_i
        ];

        % SPM[Xsrc dot sigma]
        % dSx: np x np;
        % Xsrc_x: np x np;
        % [ dSx ]                        [ sigma_x ]
        % [ dSy ][ Xsrc_x Xsrc_y Xsrc_z ][ sigma_y ]
        % [ dSz ]                        [ sigma_z ]
        M3 = [SP_x_i; SP_y_i; SP_z_i] * [Dx, Dy, Dz];

        % Collect all of the results
        M_cells{i} = M1 - M2 + M3;
    end

    M = 0.5 * blkdiag(M_cells{:});
end

function M = TSLmatrix(IDparams, p, ns)
    np = 2*p*(p + 1);
    M_cells = cell(1, ns);

    for i = 1:ns
        nu_x_cells{i} = repmat([1,0,0], np, 1);
        nu_y_cells{i} = repmat([0,1,0], np, 1);
        nu_z_cells{i} = repmat([0,0,1], np, 1);
    end

    [SP_x_cells, SP_y_cells, SP_z_cells] = spheroidalSP(IDparams, [], nu_x_cells, nu_y_cells, nu_z_cells);

    % Build spheroidalgraddiv operator
    I = repmat(eye(np),1,1,ns);
    Z = repmat(zeros(np),1,1,ns);
    ddS_x_cells = cell(1, ns); ddS_y_cells = ddS_x_cells; ddS_z_cells = ddS_x_cells;
    [D1_U, D1_V, D1_PHI] = spheroidalgraddivSL(IDparams, I, Z, Z, []);
    [D2_U, D2_V, D2_PHI] = spheroidalgraddivSL(IDparams, Z, I, Z, []);
    [D3_U, D3_V, D3_PHI] = spheroidalgraddivSL(IDparams, Z, Z, I, []);
    for i=1:ns
        if IDparams.oblate(i)
            S_self = oblate_spheroid_shape(p, IDparams.u0(i), IDparams.a(i), 'spheroidal');
            u = S_self(:,1); v = S_self(:,2); phi = S_self(:,3);
            [DXX, DYX, DZX] = convert_from_oblate_basis(D1_U(:,:,i), D1_V(:,:,i), D1_PHI(:,:,i), u, v, phi);
            [DXY, DYY, DZY] = convert_from_oblate_basis(D2_U(:,:,i), D2_V(:,:,i), D2_PHI(:,:,i), u, v, phi);
            [DXZ, DYZ, DZZ] = convert_from_oblate_basis(D3_U(:,:,i), D3_V(:,:,i), D3_PHI(:,:,i), u, v, phi);
        else
            S_self = prolate_spheroid_shape(p, IDparams.u0(i), IDparams.a(i), 'spheroidal');
            u = S_self(:,1); v = S_self(:,2); phi = S_self(:,3);
            [DXX, DYX, DZX] = convert_from_prolate_basis(D1_U(:,:,i), D1_V(:,:,i), D1_PHI(:,:,i), u, v, phi);
            [DXY, DYY, DZY] = convert_from_prolate_basis(D2_U(:,:,i), D2_V(:,:,i), D2_PHI(:,:,i), u, v, phi);
            [DXZ, DYZ, DZZ] = convert_from_prolate_basis(D3_U(:,:,i), D3_V(:,:,i), D3_PHI(:,:,i), u, v, phi);
        end

        ddS_x_cells{i} = [DXX, DXY, DXZ];
        ddS_y_cells{i} = [DYX, DYY, DYZ];
        ddS_z_cells{i} = [DZX, DZY, DZZ];
    end

    % Build TLP matvec kernel
    for i=1:ns
        %%%
        %%% Setup
        %%%
        % Grab source points
        [~, X_src_i] = IDparams.get_X(i);
        N_src_i = IDparams.get_Norm(IDparams.p, i);
        Xi_dot_Ni_vec = repmat(dot(X_src_i, N_src_i, 2), 3, 1); % 3np x 1 vector.

        % Extract SP terms
        SP_x_i = SP_x_cells{i};
        SP_y_i = SP_y_cells{i};
        SP_z_i = SP_z_cells{i};

        % Extract ddS terms
        ddS_x_i = ddS_x_cells{i};
        ddS_y_i = ddS_y_cells{i};
        ddS_z_i = ddS_z_cells{i};
        ddS_op = [ddS_x_i; ddS_y_i; ddS_z_i];

        %%%
        %%% Build matrix
        %%%
        % Extract SP sum term
        M1_diag_term = diag(N_src_i(:,1))*SP_x_i + diag(N_src_i(:,2))*SP_y_i + diag(N_src_i(:,3))*SP_z_i;
        M1 = kron(eye(3), M1_diag_term);

        % Extract ddS sum term
        % The idea is that we multiply on the left by the normal vector term, on the right
        % by the y_j term, and the ddS operator is represented by [ddS_x; ddS_y; ddS_z].
        M2_1 = kron(eye(3), diag(N_src_i(:,1))) * ddS_op * kron(eye(3), diag(X_src_i(:,1)));
        M2_2 = kron(eye(3), diag(N_src_i(:,2))) * ddS_op * kron(eye(3), diag(X_src_i(:,2)));
        M2_3 = kron(eye(3), diag(N_src_i(:,3))) * ddS_op * kron(eye(3), diag(X_src_i(:,3)));
        M2 = M2_1 + M2_2 + M2_3;

        % Extract ddS term
        M3 = diag(Xi_dot_Ni_vec) * ddS_op;

        M_cells{i} = M1 + M2 - M3;
    end

    M = blkdiag(M_cells{:});
end

%% Helper code
function kerneldSmtx = generate_kerneldS_mtx(X_self, np, permute_flag)
    %{
        Verify that the diagonal blocks are indeed correct.
    %}
    kerneldSmtx = kerneldS([], SurfaceSph(X_self));
    % Permute
    if permute_flag
        prm = zeros(1,3*np); 
        prm(1:np) = 1:3:3*np; prm(np+1:2*np) = 2:3:3*np; prm(2*np+1:3*np)=3:3:3*np;
        kerneldSmtx = kerneldSmtx(prm, prm);
    end
end
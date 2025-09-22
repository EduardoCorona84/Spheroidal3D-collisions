function [soln, fluxsoln, truesoln, truefluxsoln, sigma_vec, condK] = stokes_bie_problem(p, eta, ns, u0, target_distances, plt, neumann, interior)
    %{
        Inputs:
            - p -> order
            - eta -> matvec eta (used to distinguish points that are close/far)
            - ns -> number of spheroids
            - u0 -> eccentricities of spheroids
            - target_distances -> how far away from the surface to evaluate
            target points
            - plt (true/false) -> whether to plot the points; a holdover 
            from old code.
            - neumann (true/false) -> whether to do a neumann problem or 
            not; note that we get an integral equation of the first kind, 
            which is known to have bad conditioning
            - interior (true/false) -> whether to do an exterior problem or 
            an interior problem

        Outputs:
            - soln -> computed SL/DL result based on Dirichlet/Neumann
            problem selection
            - fluxsoln -> computed SP/DP result based on Dirichlet/Neumann
            problem selection
            - truesoln -> analytical solution (potential induced by a 
            number of point charges)
            - trueflux -> analytical flux
            - sigma_vec -> computed sigma vector
            - condK -> condition number of BIE matrix
    %}

    fluxsoln = []; truefluxsoln = [];

    %% Flags for debugging
    mix_obl = false; % Decide whether to mix oblates (or make it an oblate for the case of 1 spheroid)
    useS = true; % For exterior Dirichlet
    usekerneldS = false; % For Neumann problems
    usekernelD = false; % For interior Dirichlet problems
    
    % Flag for plotting error in the exterior Neumann case.
    PLOT_ERROR_FLAG = false;
   
    % Flag for rotating the entire coordinate system (used to debug the
    % failure of rotation for multiple bodies)
    GLOBAL_ROTATION_FLAG = false;

    %% Set up spheroid system
    np=2*p*(p+1);
    
    if length(u0)~=ns
        error("The input u0 length does not match number of spheroids indicated.");
    end
    
    pars=SpheroidalParameters;
    pars.matvec_eta = eta;
    pars.isReal = true;
    pars.u0=u0;

    if ns==3
        if mix_obl
            obl = [1,0,1];
        else
            obl = [0,0,0];
        end
        pars.oblate = obl;
        a = 1./u0;
        a(obl==1) = 1./sqrt(u0(obl==1).^2+1);
        pars.a = a;
        
        pars.centers = [0 0 0; 5 0 0; -5 0 3];
        thetas = [0 pi/10 5*pi/3];
        phis = [0 0 pi/5];
        % thetas = [0 0 0];
        % phis = [0 0 0];
    elseif ns == 2
        if mix_obl
            obl = [1,0];
        else
            obl = [0,0];
        end
        pars.oblate = obl;
        a = 1./u0;
        a(obl==1) = 1./sqrt(u0(obl==1).^2+1);
        pars.a = a;
        
        pars.centers = [0 0 0; 5 0 -5];
        % thetas = [0 0];
        % phis = [0 0];
        thetas = [0 pi/6];
        phis = [0 pi/10];
    else
        if ns~=1
            error("Number of spheroids given not implemented here.");
        end

        if mix_obl
            pars.a = 1./sqrt(u0^2+1);
            pars.oblate=1;
        else
            pars.a=1./u0;
            pars.oblate=0;
        end
        
        pars.centers=[0 0 0];
        thetas=0;
        phis=0;
    end

    Ri=zeros(3,3,ns);
    for i=1:ns
        thetai=thetas(i);
        phii=phis(i);
        Riy=[cos(thetai) 0 sin(thetai); 0 1 0; -sin(thetai) 0 cos(thetai)];
        Riz=[cos(phii) -sin(phii) 0; sin(phii) cos(phii) 0; 0 0 1];
        Ri(:,:,i)=Riz*Riy;
    end
    pars.Rmat=Ri;
    pars.thetas=thetas;
    pars.phis=phis;
    pars.sigma=zeros(np,1,ns);

    %% Global rotation (if enabled)
    % The entire idea is to rotate the entire global coordinate system.
    % Then, the result should stay the same--except that the resulting
    % answer will differ by the (global) rotation.
    if GLOBAL_ROTATION_FLAG
        fprintf("-- Global rotation flag enabled. --\n");
        theta_g = -pi/4;
        phi_g = 0;
        Riy = [cos(theta_g) 0 sin(theta_g); 0 1 0; -sin(theta_g) 0 cos(theta_g)];
        Riz = [cos(phi_g) -sin(phi_g) 0; sin(phi_g) cos(phi_g) 0; 0 0 1];
        R_global = Riz*Riy;

        fprintf("Rotating centers...\n");
        pars.centers = pars.centers * R_global';

        fprintf("Re-orienting each spheroid...\n");
        for i = 1:ns
            Ri(:,:,i) = R_global * Ri(:,:,i);
        end
        pars.Rmat = Ri;

        fprintf("-- Done re-orienting coordinate system. --\n");
    end

    %% Placement of point forces
    num_pf = 3; % Number of point forces per spheroid
    c = pars.centers;
    F_pos_vec = reshape(repmat(reshape(c',3,1,[]),1,num_pf),3,[],1)';

    switch (interior)
        case true % Interior problem (place points randomly outside)
            placement_scale = 2;
            if ns>1
                % Some scaling factor times the major radius of the prolate 
                % spheroids guarantees that we shall be outside the
                % spheroids.
                major_radii = repelem(sqrt(pars.u0.^2 + 1), num_pf)' .* repelem(pars.a, num_pf)';
                d = placement_scale * major_radii .* generate_random_unit_vec(ns*num_pf, 3);
            else
                d = [];
                for i=1:num_pf
                    d = [d ; placement_scale.*pars.a.*sqrt(pars.u0.^2 + 1).*generate_random_unit_vec(1, 3)];
                end
            end
        case false % Exterior problem (place points near center)
            if ns>1
                % A multiple of the minor radius of the prolate spheroids
                scale = repmat(0.3 .* pars.a .*sqrt(pars.u0.^2-1), num_pf, 1);
                scale = scale(:);
                d = repmat(scale,1,3) .* (rand(size(F_pos_vec))-.5);
            else
                d_z=0.2 .* (rand(num_pf,1)-.5);
                d_xy=0.05.*(rand(num_pf,2)-0.5);
                d = apply_random_rotation([d_xy d_z]);
            end
    end

    F_pos_vec = F_pos_vec + d;
    F_vec = rand(ns*num_pf, 3); % Point forces

    % Rotate everything accordingly since the above was constructed from
    % the perspective of an unrotated spheroid.
    F_pos_vec = rotate_point_forces(F_pos_vec, pars.centers, num_pf, ns, Ri);
    F_vec = rotate_point_forces(F_vec, pars.centers, num_pf, ns, Ri);

    if plt
        pars.plot(F_pos_vec);
        title("Spheroids and locations of point forces");
    end

    %% Setup BIE
    Y = pars.get_X;
    [NrY, ~] = pars.get_Norm_rot(p); % Necessary to account of the spheroid's rotation

    if ~neumann % Dirichlet
        % Find boundary condition (surface velocity induced by point forces)
        truesolnSurf = stokeslet_velocity(F_vec, F_pos_vec, Y);

        % Construct DLP on-surface matrices
        DM = L2StkMatVecKernel(pars, 'DLP', p);

        if usekernelD
            DMkernelD = kernelD([], SurfaceSph(Y));
        end

        if interior % Interior problem
            CM = nu_completion(pars, p, ns);
            K = -0.5*eye(3*ns*np) + DM + CM;
        else % Exterior problem
            % See Pozrikidis, Chapter 4.7.
            if useS
                CM = L2StkMatVecKernel(pars, 'SLP', p);
            else
                CM = RBM_completion(pars, Y, ns);
            end
            K = 0.5*eye(3*ns*np) + DM + CM;
        end

        if usekernelD
            prm = zeros(1,3*np); 
            prm(1:3:3*np) = 1:np; prm(2:3:3*np) = np+1:2*np; prm(3:3:3*np)=2*np+1:3*np;
            if interior
                interior_factor = -1;
            else
                interior_factor = 1;
            end
            KkernelDM = interior_factor*0.5*eye(3*ns*np) + DMkernelD;
        end
    else % Neumann problem
        truesolnSurf = stresslet_traction(F_vec, F_pos_vec, Y, NrY);
        if usekerneldS
            TMkerneldS = kerneldS([], SurfaceSph(Y));
        end
        if interior
            % See Pozrikidis, Chapter 4.2. The completion term should be
            % the projection onto the space representing the 6 rigid body
            % motion operators.
            TM = L2StkMatVecKernel(pars, 'TLP', p, NrY);
            CM = RBM_completion(pars, Y, ns);
            K = 0.5*eye(3*ns*np) + TM + CM;
        else
            % See Hsiao and Wedland, Chapter 2.3. The completion term
            % should be the projection onto the space of surface normals
            % associated with every body.
            TM = L2StkMatVecKernel(pars, 'TLP', p, NrY);
            CM = nu_completion(pars, p, ns);
            K = -0.5*eye(3*ns*np) + TM + CM;
        end

        if usekerneldS
            prm = zeros(1,3*np); 
            prm(1:3:3*np) = 1:np; prm(2:3:3*np) = np+1:2*np; prm(3:3:3*np)=2*np+1:3*np;
            CMkdS = CM(prm, prm);
            if interior
                interior_factor = 1;
            else
                interior_factor = -1;
            end
            KkerneldS = interior_factor*0.5*eye(3*ns*np) + TMkerneldS + CMkdS;
        end
    end

    %% Processing to setup linear equation
    condK = cond(K);

    if (neumann && usekerneldS) || (~neumann && usekernelD)
        kernel_truesolnSurf = reshape([truesolnSurf(:,1), truesolnSurf(:,2), truesolnSurf(:,3)].', [], 1);
    end

    % Move data (i.e. surface velocity/traction) to local frame
    truesolnSurf = rotate_data(truesolnSurf, np, ns, Ri);

    tsSurf = [];
    for i=1:ns
        truesolnSurf_particle_i = truesolnSurf((i-1)*np+1 : i*np,:);
        tsSurf = [
            tsSurf;
            truesolnSurf_particle_i(:,1);
            truesolnSurf_particle_i(:,2);
            truesolnSurf_particle_i(:,3)
        ];
    end

    %% Solve for density
    sigma_vec = gmres(K, tsSurf, 1000, 1e-14);
    % sigma_vec = K \ tsSurf;

    if neumann && usekerneldS
        kerneldS_sigma_vec = KkerneldS \ kernel_truesolnSurf;
    elseif ~neumann && usekernelD
        kernelD_sigma_vec = KkernelDM \ kernel_truesolnSurf;
    end

    % Unpack density
    sigma_x = zeros(np, 1, ns); sigma_y = sigma_x; sigma_z = sigma_x;
    for i=1:ns
        sigma_particle_i = sigma_vec(3*(i-1)*np+1 : 3*i*np);
        sigma_x(:,:,i) = sigma_particle_i(1:np);
        sigma_y(:,:,i) = sigma_particle_i(np+1:2*np);
        sigma_z(:,:,i) = sigma_particle_i(2*np+1:3*np);
    end

    N = ns*np;
    if neumann && usekerneldS && ns==1
        idx = 1:3:3*N; idy = 2:3:3*N; idz = 3:3:3*N;
        kdS_sigma_x = kerneldS_sigma_vec(idx); kdS_sigma_y = kerneldS_sigma_vec(idy); kdS_sigma_z = kerneldS_sigma_vec(idz);
    elseif ~neumann && usekernelD && ns==1
        idx = 1:3:3*N; idy = 2:3:3*N; idz = 3:3:3*N;
        kD_sigma_x = kernelD_sigma_vec(idx); kD_sigma_y = kernelD_sigma_vec(idy); kD_sigma_z = kernelD_sigma_vec(idz);
    elseif ns==1
        sigma_x = sigma_vec(1:N); sigma_y = sigma_vec(N+1:2*N); sigma_z = sigma_vec(2*N+1:end);
    end

    %%% Target points to test at
    Xcell=cell(1,ns);

    if interior
        interior_factor = -1;
    else
        interior_factor = 1;
    end

    for i=1:ns
        if ~pars.oblate(i)
            Yi = prolate_spheroid_shape(p, pars.u0(i), pars.a(i), 'cart');
        else
            Yi = oblate_spheroid_shape(p, pars.u0(i), pars.a(i), 'cart');
        end
        Xcell{i} = Yi + interior_factor * target_distances(i).*get_norm_vecs(p, pars.u0(i), pars.oblate(i));
    end
    Xeval = pars.set_X_targets(Xcell); % Re-orient target points

    if plt
        pars.plot(Xeval);
        title('spheroids and target points')
    end

    %% Now, compute the solution.
    if ~neumann % Dirichlet problem
        if interior % Interior problem
            [Dterm_x, Dterm_y, Dterm_z] = L2StkMatVec(pars, 'DLP', sigma_x, sigma_y, sigma_z, Xeval);
            soln = [Dterm_x, Dterm_y, Dterm_z];

            if usekernelD
                [kD_Sterm_x, kD_Sterm_y, kD_Sterm_z] = L2StkMatVec(pars, 'DLP', kD_sigma_x, kD_sigma_y, kD_sigma_z, Xeval);
                solnkD = [kD_Sterm_x, kD_Sterm_y, kD_Sterm_z];
                truesoln = stokeslet_velocity(F_vec, F_pos_vec, Xeval);
                fprintf('p=%d: kernelD eval comparison = %.6e\n',p, norm(truesoln - solnkD) ./ norm(truesoln));

                sigma_rel_err = norm([sigma_x, sigma_y, sigma_z] - [kD_sigma_x, kD_sigma_y, kD_sigma_z])/norm([sigma_x, sigma_y, sigma_z]);
                fprintf('p=%d: difference between densities: %.6e\n', p, sigma_rel_err);
            end
        else % Exterior problem
            [Dterm_x, Dterm_y, Dterm_z] = L2StkMatVec(pars, 'DLP', sigma_x, sigma_y, sigma_z, Xeval);
            if useS
                [Sterm_x, Sterm_y, Sterm_z] = L2StkMatVec(pars, 'SLP', sigma_x, sigma_y, sigma_z, Xeval);
                soln = [Dterm_x + Sterm_x, Dterm_y + Sterm_y, Dterm_z + Sterm_z];
            else
                soln = [Dterm_x, Dterm_y, Dterm_z];
            end
        end
    else % Neumann problem
        [Sterm_x, Sterm_y, Sterm_z] = L2StkMatVec(pars, 'SLP', sigma_x, sigma_y, sigma_z, Xeval);
        soln = [Sterm_x, Sterm_y, Sterm_z];

        [TSLterm_x, TSLterm_y, TSLterm_z] = L2StkMatVec(pars, 'TLP', sigma_x, sigma_y, sigma_z, Xeval, NrY);
        fluxsoln = [TSLterm_x, TSLterm_y, TSLterm_z];

        truefluxsoln = stresslet_traction(F_vec, F_pos_vec, Xeval, NrY);

        fprintf('p=%d: flux comparison = %.6e\n',p, norm(truefluxsoln - fluxsoln) ./ norm(truefluxsoln)); 

        if usekerneldS
            [kdS_Sterm_x, kdS_Sterm_y, kdS_Sterm_z] = L2StkMatVec(pars, 'SLP', kdS_sigma_x, kdS_sigma_y, kdS_sigma_z, Xeval);
            solnkdS = [kdS_Sterm_x, kdS_Sterm_y, kdS_Sterm_z];
    
            [kdS_TSLterm_x, kdS_TSLterm_y, kdS_TSLterm_z] = L2StkMatVec(pars, 'TLP', kdS_sigma_x, kdS_sigma_y, kdS_sigma_z, Xeval, NrY);
            fluxsolnkdS = [kdS_TSLterm_x, kdS_TSLterm_y, kdS_TSLterm_z];

            fprintf('p=%d: kerneldS flux comparison = %.6e\n',p, norm(truefluxsoln - fluxsolnkdS) ./ norm(truefluxsoln));
            truesoln = stokeslet_velocity(F_vec, F_pos_vec, Xeval);
            fprintf('p=%d: kerneldS velocity comparison = %.6e\n',p, norm(truesoln - solnkdS) ./ norm(truesoln));
        end
    end

    %% Finally, we compare our computed solutions to the true solutions.
    truesoln = stokeslet_velocity(F_vec, F_pos_vec, Xeval);

    % Rotate solution back to correct coordinate frame
    soln = rotate_soln(soln, np, ns, Ri);

    fprintf('p=%d: velocity comparison = %.6e\n',p, norm(truesoln - soln) ./ norm(truesoln));

    if ns > 1
        return;
    end

    %% Spurious RBM -- fails for rotated bodies
   fprintf('Examining whether the computed velocity differs by an RBM...\n');
   
    v_RBM = project_onto_RBM(pars, soln - truesoln, Xeval);
    if ~neumann && usekernelD
        % to fix for rotation, need to pass Xeval rotated back in canonical
        % frame?
        v_RBM_kernel = project_onto_RBM(pars, solnkD - truesoln, Xeval);
        rel_vel_error = norm(solnkD - truesoln - v_RBM_kernel);
        fprintf('p=%d: eval comparison for kernelD after RBM adjustment = %.6e\n',p, rel_vel_error);
    elseif neumann && usekerneldS
        v_RBM_kernel = project_onto_RBM(pars, solnkdS - truesoln, Xeval);
        rel_vel_error = norm(solnkdS - truesoln - v_RBM_kernel);
        fprintf('p=%d: eval comparison for kerneldS after RBM adjustment = %.6e\n',p, rel_vel_error);
    end

    error_aligned = norm(soln - v_RBM - truesoln) / norm(truesoln);
    
    fprintf('p=%d: velocity comparison after RBM adjustment = %.6e\n',p, error_aligned);

    %% Plot debug code
    % Plot heat map of traction error
    if neumann
        scatter3(Xeval(:,1), Xeval(:,2), Xeval(:,3), 40, log10(abs(vecnorm(truefluxsoln - fluxsoln, 2, 2))), 'filled');
        colorbar; colormap(jet);
    end

    % Plot error
    if PLOT_ERROR_FLAG && ns==1 && (~interior && neumann) % Only do this for exterior Neumann.
        distances = [10.^(0:-1:-5)];
        errs_vel = cell(1, numel(distances));
        errs_traction = cell(1, numel(distances));
        aspect_ratio = abs(pars.u0/sqrt(pars.u0.^2 - 1));

        %% Computation
        for i=1:numel(distances)
            d = distances(i);
            if pars.oblate
                X_trg = oblate_spheroid_shape(p, pars.u0, pars.a);
            else
                X_trg = prolate_spheroid_shape(p, pars.u0, pars.a);
            end

            nu_trg = get_norm_vecs(p, pars.u0, pars.oblate);
            X_trg = X_trg + d*nu_trg;

            [vx, vy, vz] = L2StkMatVec(pars, 'SLP', sigma_x, sigma_y, sigma_z, X_trg);
            true_velocity = stokeslet_velocity(F_vec, F_pos_vec, X_trg);

            [tx, ty, tz] = L2StkMatVec(pars, 'TLP', sigma_x, sigma_y, sigma_z, X_trg, NrY);
            true_traction = stresslet_traction(F_vec, F_pos_vec, X_trg, NrY);

            errs_vel{i} = log10(abs([vx, vy, vz] - true_velocity) ./ abs(true_velocity));
            errs_traction{i} = log10(abs([tx, ty, tz] - true_traction) ./ abs(true_traction));
        end

        %% Plot errors
        % Find correct padding
        min_y = inf;
        max_y = -inf;
        for j = 1:numel(distances)
            all_data_at_dist = [errs_vel{j}(:); errs_traction{j}(:)];
            if ~isempty(all_data_at_dist)
                min_y = min(min_y, min(all_data_at_dist));
                max_y = max(max_y, max(all_data_at_dist));
            end
        end
        padding = (max_y - min_y) * 0.05;
        y_limits = [min_y - padding, max_y + padding];
    
        % Create a figure for the current p-value
        figure('Name', sprintf('p = %d', p), 'Position', [100, 100, 1200, 500]);
        sgtitle(sprintf('evaluation for exterior Neumann problem with p=%d, aspect ratio='+string(aspect_ratio), p));
        
        num_plots = numel(distances);
        
        % Loop through distances to create a subplot for each
        for j = 1:numel(distances)
            d = distances(j);
            subplot(1, num_plots, j);
            
            data1 = errs_vel{j}(:);
            data2 = errs_traction{j}(:);
            combined_data = [data1; data2];
            grouping_variable = [
                repmat({'SLP'}, numel(data1), 1) ;
                repmat({'TSL'}, numel(data2), 1)
            ];
            
            % Create the box chart
            boxchart(ones(size(combined_data)), combined_data, 'GroupByColor', categorical(grouping_variable));
            
            % --- Formatting ---
            xticklabels('')
            ylim(y_limits);
            title(sprintf('d=%.0d', d));
            grid on;
    
            if j == 1
                ylabel('log10(relative error)');
            end
        end
    
        legend
    end
end

function u = stokeslet_velocity(F_vec, F_pos_vec, x)
    %{
        Calculates the resulting velocity due to a number of point forces.
        Inputs
            -   F_vec : strength of point forces
            -   F_pos_vec : 3D positions of point forces
            -   x : evaluation point
        Outputs
            -   u : resulting velocity vector (np x 3)
    %}
    num_pf = size(F_vec, 1); % number of point forces
    np = size(x, 1);
    u = zeros(np, 3);

    % Note that (r \oplus r) F = (r \cdot F) * r
    for i=1:num_pf
        F = F_vec(i,:);
        y = F_pos_vec(i,:);
        
        r = x - y;
        normR = sqrt(sum(r.^2, 2));

        u = u + (1/(8*pi)) * ( F./normR + (sum(r .* F, 2).*r)./(normR.^3) );
    end
end

function t = stresslet_traction(F_vec, F_pos_vec, x, n_x)
    %{
        Sum up Stresslets at evaluation points due to given point forces.

        Inputs
            - F_vec : strength of point forces
            - F_pos_vec : 3D positions of point forces
            - x : evaluation point
            - n_x : surface normals at each eval point
        Outputs  
            - flux : resulting traction on surface (np x 3)
    %}
    num_pf = size(F_vec, 1); % number of point forces
    np = size(x, 1);
    t = zeros(np, 3);

    for i=1:num_pf
        F = F_vec(i,:);
        y = F_pos_vec(i,:);
        
        r = x - y;
        normR = sqrt(sum(r.^2, 2));
        t = t + -(3/(4*pi)) * (sum(r .* F, 2) .* sum(r .* n_x, 2)) ./ (normR.^5) .* r;
    end
end

function unit_vectors = generate_random_unit_vec(num_vectors, vector_length)
    unit_vectors = zeros(num_vectors, vector_length);
    for i=1:num_vectors
        random_vector = randn(1, vector_length);
        unit_vectors(i,:) = random_vector / norm(random_vector);
    end
end

function res_vec = apply_random_rotation(pos_vecs)
    %{
        Apply a random rotation to each row of pos_vecs (nr x 3)
    %}
    [Q, ~]  = qr(randn(3));
    res_vec = pos_vecs * Q;
end

function v_RBM = project_onto_RBM(pars, v, Xeval)
    %{
        Project onto the 6-dimensional subspace formed by the rigid body
        motion vectors.
    %}
    N = size(Xeval, 1);
    
    % Center of the spheroid (for rotation)
    center = pars.centers(1,:);
    R = Xeval - center;
    
    % The columns represent the basis for the 6 RBMs.
    RBM_mtx = zeros(3*N, 6);
    
    % Translation
    RBM_mtx(1:N, 1) = 1;
    RBM_mtx(N+1:2*N, 2) = 1;
    RBM_mtx(2*N+1:end, 3) = 1;
    
    % Rotation
    % [0, -R_z, R_y]
    RBM_mtx(N+1:2*N, 4) = -R(:,3);
    RBM_mtx(2*N+1:end, 4) = R(:,2);
    % [R_z, 0, -R_x]
    RBM_mtx(1:N, 5) = R(:,3);
    RBM_mtx(2*N+1:end, 5) = -R(:,1);
    % [-R_y, R_x, 0]
    RBM_mtx(1:N, 6) = -R(:,2);
    RBM_mtx(N+1:2*N, 6) = R(:,1);
    
    b = [v(:,1); v(:,2); v(:,3)];
    rbm_coeffs = RBM_mtx \ b;
    trans_coeffs = rbm_coeffs(1:3)';
    rot_coeffs = rbm_coeffs(4:6)';
    
    v_RBM = RBM_mtx * rbm_coeffs;
    v_RBM = reshape(v_RBM, N, 3);
end

function [total_force, total_torque] = calculate_force_and_torque(pars, traction, Y, p, ns)
    % Note that the input "traction" here is assumed to be np x 3 (for 1
    % spheroid only).

    np = 2*p*(p+1);
    [~, gwt]=g_grid(p+1);
    wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
    wt = wt(:);

    total_torque = zeros(1,3);
    total_force = zeros(1,3);
    for i=1:ns
        % Grab information for current spheroid
        source_indices = (i-1)*np+1 : i*np;
        Yi = Y(source_indices, :);
        traction_i = traction(source_indices, :);
        
        % Surface area element
        Sns = SurfaceSph(Yi(:));
        W = Sns.geoProp.W;
        
        % Scaled quadrature weights
        Wns = W.*wt;
        
        % Calculate net torque for this spheroid (relative to its center)
        ci = pars.centers(i,:);
        R = Yi - ci;
        net_torque_i = sum(cross(R, traction_i) .* Wns, 1);
        total_torque = total_torque + net_torque_i;
        
        % Calculate net force for this spheroid
        total_force = total_force + Wns.' * traction_i;
    end

    fprintf("Norm of force on surface: %e\n", norm(total_force));
    fprintf('Norm of net torque on surface: %e\n', norm(total_torque));
end

function plot_singular_values(mtx)
    [U, S, V] = svd(mtx);
    s_vals = diag(S);

    rank_no_tol = rank(mtx);
    rank_tol = rank(mtx, 1e-6);
    fprintf("Rank (no given tolerance): %d -- nullspace dimension is %d\n", rank(mtx), size(mtx,1) - rank_no_tol);
    fprintf("Rank (tolerance of 1e-6): %d -- nullspace dimension is %d\n", rank_tol, size(mtx,1) - rank_tol);
    
    figure;
    semilogy(s_vals, 'o-');
    title('Singular Values of Matrix');
    xlabel('Index');
    ylabel('Singular Value');
    grid on;
end

function plot_eigenspace(mtx)
end

function res = rotate_point_forces(vec, centers, num_pf, ns, R)
    %{
        Helper function to apply rotation to a list of vectors (this is
        mostly to handle the case of multiple spheroids.)

        The point forces are respectively rotated around each spheroid's
        center.
    %}
    res = zeros(size(vec));
    for i=1:ns
        pf_indices = (i-1)*num_pf+1:i*num_pf;
        res(pf_indices,:) = (vec(pf_indices,:) - centers(i))*R(:,:,i)' + centers(i);
    end
end

function res = rotate_data(vec, np, ns, R)
    %{
        In the case of rotated bodies, one must also rotate the data (i.e.
        the velocity and traction) to the local frame of the associated
        spheroid (i.e. upright).

        This is why we multiply on the right by R, and not the transpose
        R^T.
    %}
    res = zeros(size(vec));
    for i=1:ns
        spheroid_block = (i-1)*np+1:i*np;
        res(spheroid_block,:) = vec(spheroid_block,:)*R(:,:,i);
    end
end

function res = rotate_soln(vec, np, ns, R)
    %{
        The solution that is calculated is in the local frame, so we need
        to bring it back to the correct frame by rotating the solutions
        back.
    %}
    res = zeros(size(vec));
    for i=1:ns
        spheroid_block = (i-1)*np+1:i*np;
        res(spheroid_block,:) = vec(spheroid_block,:)*R(:,:,i)';
    end
end

function CM = nu_completion(pars, p, ns)
    %{
        Returns the completion term that corresponds to the space spanned
        by the normals on every spheroid.
    %}
    np = 2*p*(p+1);
    CM_cells = cell(1, ns);
    for i=1:ns
        NrY_i = get_norm_vecs(p, pars.u0(i), pars.oblate(i));
        NrY_i_stacked = [NrY_i(:,1) ; NrY_i(:,2) ; NrY_i(:,3)];
        CM_cells{i} = 1/(norm(NrY_i_stacked)^2) * (NrY_i_stacked*NrY_i_stacked.');
    end
    CM = blkdiag(CM_cells{:});
end
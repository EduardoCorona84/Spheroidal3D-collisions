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

    %%% Backwards compatibility
    Xeval_p = p;
    mix_obl = true;
    useS = false;

    %%% Set up spheroid system
    np=2*p*(p+1);
    
    if length(u0)~=ns
        fprintf("\n input u0 length does not match number of spheroids indicated.");
        ns=length(u0);
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
        
        pars.centers = [0 0 0; 5 0 0; 3.2 3.2 3.2];
        thetas = [0 pi/10 5*pi/3];
        phis = [0 0 pi/5];
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
        
    %%% Point charges
    % 2 point charges per spheroid
    num_pf = 2;
    c = pars.centers;
    F_pos_vec = reshape(repmat(reshape(c',3,1,[]),1,num_pf),3,[],1)';

    %%% Placement of point forces
    switch (interior)
        case true % Interior problem (place points randomly outside)
            if ns==3
                % Some scaling factor times the major radius of the prolate 
                % spheroids guarantees that we shall be outside the
                % spheroids.
                placement_scale = 2;
                major_radii = repelem(sqrt(pars.u0.^2 + 1), num_pf)' .* repelem(pars.a, num_pf)';
                d = placement_scale * major_radii .* generate_random_unit_vec(6, 3);
            else
                placement_scale = 2;
                v1 = placement_scale.*pars.a.*sqrt(pars.u0.^2 + 1).*generate_random_unit_vec(1, 3);
                v2 = placement_scale.*pars.a.*sqrt(pars.u0.^2 + 1).*generate_random_unit_vec(1, 3);
                d = [v1 ; v2];
            end
        case false % Exterior problem (place points near center)
            if ns==3
                % 1/2 of the minor radius of the prolate spheroids
                scale = repmat(0.3.*pars.a .*sqrt(pars.u0.^2-1),num_pf,1);
                scale = scale(:);
                d = repmat(scale,1,3) .* (rand(size(F_pos_vec))-.5);
            else
                d_z=0.5.*(rand(num_pf,1)-.5);
                d_xy=0.01.*(rand(num_pf,2)-0.5);
                d = [d_xy d_z];
            end
    end

    F_pos_vec = F_pos_vec + d;
    F_vec = rand(ns*num_pf, 3); % Point forces

    if plt
        pars.plot(F_pos_vec);
        title("Spheroids and locations of point forces");
    end

    Y = pars.get_X;
    [NrY, ~] = pars.get_Norm_rot(p); % Necessary to account of the spheroid's rotation

    if ~neumann
        % Find boundary condition (surface velocity induced by point forces)
        truesolnSurf = stokeslet_velocity(F_vec, F_pos_vec, Y);

        % Construct DL on-surface matrices
        DM = L2StkMatVecKernel(pars, 'DLP', p);

        if interior % Interior problem -> need to account for compatibility condition?
            N_stacked = [NrY(:,1) ; NrY(:,2) ; NrY(:,3)];
            CM = 1/(norm(N_stacked)^2) * (N_stacked*N_stacked.'); % In this case, the norm squared should be 3*np.
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
    else % Neumann problem
        truesolnSurf = stokeslet_traction(F_vec, F_pos_vec, Y, NrY);
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
            % assocaited with every body.
            TM = L2StkMatVecKernel(pars, 'TLP', p, NrY);
            N_stacked = [NrY(:,1) ; NrY(:,2) ; NrY(:,3)];
            CM = 1/(norm(N_stacked)^2) * (N_stacked*N_stacked.'); % In this case, the norm squared should be 3*np.
            K = -0.5*eye(3*ns*np) + TM + CM;
        end
    end

    condK = cond(K);
    
    truesolnSurf = [truesolnSurf(:,1) ; truesolnSurf(:,2) ; truesolnSurf(:,3)];
    sigma_vec = gmres(K, truesolnSurf, 1000, 1e-12);
    N = ns*np;
    sigma_x = sigma_vec(1:N); sigma_y = sigma_vec(N+1:2*N); sigma_z = sigma_vec(2*N+1:end);
    
    %%% Target points to test at
    Xcell=cell(1,ns);

    if interior
        interior_factor = -1;
    else
        interior_factor = 1;
    end

    for i=1:ns
        if ~pars.oblate(i)
            Yi = prolate_spheroid_shape(Xeval_p, pars.u0(i), pars.a(i), 'cart');
        else
            Yi = oblate_spheroid_shape(Xeval_p, pars.u0(i), pars.a(i), 'cart');
        end
        Xcell{i} = Yi + interior_factor * target_distances(i).*get_norm_vecs(Xeval_p, pars.u0(i), pars.oblate(i));
    end
    Xeval = pars.set_X_targets(Xcell);

    if plt
        pars.plot(Xeval);
        title('spheroids and target points')
    end
    
    %%% Now, compute the solution.
    if ~neumann % Dirichlet problem
        if interior % Interior problem
            [Dterm_x, Dterm_y, Dterm_z] = L2StkMatVec(pars, 'DLP', sigma_x, sigma_y, sigma_z, Xeval);
            soln = [Dterm_x, Dterm_y, Dterm_z];
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
        if interior
            [Sterm_x, Sterm_y, Sterm_z] = L2StkMatVec(pars, 'SLP', sigma_x, sigma_y, sigma_z, Xeval);
            soln = [Sterm_x, Sterm_y, Sterm_z];

            [TSLterm_x, TSLterm_y, TSLterm_z] = L2StkMatVec(pars, 'TLP', sigma_x, sigma_y, sigma_z, Xeval, NrY);
            fluxsoln = [TSLterm_x, TSLterm_y, TSLterm_z];

            truefluxsoln = stokeslet_traction(F_vec, F_pos_vec, Xeval, NrY);

            fprintf('p=%d: flux comparison = %.6e\n',p, norm(truefluxsoln - fluxsoln) ./ norm(truefluxsoln)); 
        else
            [Sterm_x, Sterm_y, Sterm_z] = L2StkMatVec(pars, 'SLP', sigma_x, sigma_y, sigma_z, Xeval);
            soln = [Sterm_x, Sterm_y, Sterm_z];

            [TSLterm_x, TSLterm_y, TSLterm_z] = L2StkMatVec(pars, 'TLP', sigma_x, sigma_y, sigma_z, Xeval, NrY);
            fluxsoln = [TSLterm_x, TSLterm_y, TSLterm_z];

            truefluxsoln = stokeslet_traction(F_vec, F_pos_vec, Xeval, NrY);

            fprintf('p=%d: flux comparison = %.6e\n',p, norm(truefluxsoln - fluxsoln) ./ norm(truefluxsoln)); 
        end
    end

    %%% Finally, we compare our computed solutions to the true solutions.
    truesoln = stokeslet_velocity(F_vec, F_pos_vec, Xeval);
    fprintf('p=%d: eval comparison = %.6e\n',p, norm(truesoln - soln) ./ norm(truesoln));
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

function t = stokeslet_traction(F_vec, F_pos_vec, x, n_x)
    %{
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

function CM = RBM_completion(pars, Y, ns)
    %{
        Calculates the completion term that corresponds to rigid body
        motion.
    %}
    p = pars.p; np = 2*p*(p+1);

    % Unit quadrature weights
    [~, gwt]=g_grid(p+1);
    wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
    wt = wt(:);

    CM_cells = cell(1, ns);
    for i=1:ns
        if ~pars.oblate(i)
            Yi = prolate_spheroid_shape(p,pars.u0(i),pars.a(i),'cart');
        else
            Yi = oblate_spheroid_shape(p,pars.u0(i),pars.a(i),'cart');
        end
        ci = pars.centers(i,:);

        % Surface area weights
        Sns = SurfaceSph(Yi(:));
        W = Sns.geoProp.W;

        % Scaled quadrature weights
        Wns = W.*wt;

        % Calculate moment of inertia tensor \tau
        R = Yi - ci;
        M = R' * (Wns .* R);
        tau = trace(M) * eye(3) - M;

        % Rotational portion
        rotational_integral_term = [
            zeros(1, np),       -1*(Wns.*R(:,3)).', (Wns.*R(:,2)).';
            (Wns .* R(:,3)).',    zeros(1, np),     -1*(Wns.*R(:,1)).';
            -1*(Wns.*R(:,2)).',   (Wns.*R(:,1)).',    zeros(1, np)
        ]; % 3 x 3np (input is density)

        % This represents the cross-product with the integral (accounting 
        % for the anti-commutativity),
        cross_prod_term = -1*[
            zeros(np, 1), -1*R(:,3),    R(:,2);
            R(:,3),       zeros(np, 1), -1*R(:,1);
            -1*R(:,2),    R(:,1),       zeros(np, 1)
        ]; % 3np x 3

        CM_rot = cross_prod_term * (tau \ rotational_integral_term); % 3np x 3np

        % Calculate translational portion
        surf_area = sum(Wns);
        T = (1/surf_area) * ones(np, 1) * Wns';
        CM_trans = kron(eye(3), T); % 3np x 3np

        CM_cells{i} = CM_trans + CM_rot;
    end
    CM = blkdiag(CM_cells{:});
end
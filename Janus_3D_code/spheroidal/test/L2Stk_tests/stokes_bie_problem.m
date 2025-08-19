function [soln, truesoln, sigma_vec, condK] = stokes_bie_problem(p, eta, ns, u0, target_distances, plt, neumann, interior)
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

    %%% Backwards compatibility
    Xeval_p = p;
    mix_obl = true;

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
        % pars.centers = [0 0 0; 2.5 0 0; 1.6 1.6 1.6];
        thetas = [0 pi/10 5*pi/3];
        phis = [0 0 pi/5];
        % thetas = [0 0 0];
        % phis = [0 0 0];
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
    nc = 2;
    c = pars.centers;
    F_pos_vec = reshape(repmat(reshape(c',3,1,[]),1,nc),3,[],1)';

    %%% Placement of point forces
    switch (interior)
        case true % Interior problem (place points randomly outside)
            if ns==3
                % Some scaling factor times the major radius of the prolate 
                % spheroids guarantees that we shall be outside the
                % spheroids.
                placement_scale = 2;
                major_radii = repelem(sqrt(pars.u0.^2 + 1), nc)' .* repelem(pars.a, nc)';
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
                scale = repmat(0.3.*pars.a .*sqrt(pars.u0.^2-1),nc,1);
                scale = scale(:);
                d = repmat(scale,1,3) .* (rand(size(F_pos_vec))-.5);
            else
                d_z=0.5.*(rand(nc,1)-.5);
                d_xy=0.01.*(rand(nc,2)-0.5);
                d = [d_xy d_z];
            end
        F_pos_vec = F_pos_vec + d;
    end

    F_vec = rand(ns*nc, 3); % Point forces

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

        if interior % Interior problem
            % -1/2 I + D has a 1-dimensional kernel
            K = -0.5*eye(3*ns*np) + DM;
        else % Exterior problem
            S = L2StkMatVecKernel(pars, 'SLP', p);
            K = 0.5*eye(3*ns*np) + DM + S;
        end
    else % Neumann problem
        error("not implemented.");
    end

    condK = cond(K);
    
    % truesolnSurf = reshape(truesolnSurf.', [], 1);
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
    soln = zeros(np, 3);
    if ~neumann % Dirichlet problem
        if interior % Interior problem
            [Dterm_x, Dterm_y, Dterm_z] = L2StkMatVec(pars, 'DLP', sigma_x, sigma_y, sigma_z, Xeval);
            soln = [Dterm_x, Dterm_y, Dterm_z];
        else % Exterior problem
            [Dterm_x, Dterm_y, Dterm_z] = L2StkMatVec(pars, 'DLP', sigma_x, sigma_y, sigma_z, Xeval);
            [Sterm_x, Sterm_y, Sterm_z] = L2StkMatVec(pars, 'SLP', sigma_x, sigma_y, sigma_z, Xeval);
            soln = [Dterm_x + Sterm_x, Dterm_y + Sterm_y, Dterm_z + Sterm_z];
        end
    else % Neumann problem
        error("not implemented.");
    end

    %%% Finally, we compare our computed solutions to the true solutions.
    truesoln = stokeslet_velocity(F_vec, F_pos_vec, Xeval);
    fprintf('p=%d: surface comparison = %.6e\n',p, norm(truesoln - soln) ./ norm(soln)); 
end

function u = stokeslet_velocity(F_vec, F_pos_vec, x)
    %{
        Calculates the resulting velocity due to a number of point forces.
        Inputs
            -   F_vec : strength of point forces
            -   F_pos_vec : 3D positions of point forces
            -   x : evaluation point
        Outputs
            -   u : resulting velocity vector (1 x 3)
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

function flux = stokeslet_traction(F_vec, F_pos_vec, x, n_x)
    %{
        Inputs
            - F_vec : strength of point forces
            - F_pos_vec : 3D positions of point forces
            - x : evaluation point
            - n_x :
        Outputs  
    %}
    num_pf = size(F_vec, 1); % number of point forces
    np = size(x, 1);
    R = zeros(np, num_pf);
    % \|x - y\|
    for i=1:num_pf
        R(:,i) = sqrt(sum((F_pos_vec(i,:) - x).^2, 2));
    end
    flux=EdotN./(4*pi*R.^2)*ptch;
end

function unit_vectors = generate_random_unit_vec(num_vectors, vector_length)
    unit_vectors = zeros(num_vectors, vector_length);
    for i=1:num_vectors
        random_vector = randn(1, vector_length);
        unit_vectors(i,:) = random_vector / norm(random_vector);
    end
end

function CM = calculate_completion_term(pars, Y, ns)
    %{
        Calculates the completion term that is necessary for the exterior
        Dirichlet problems (the completion term here should be the matrix
        that represents rigid body motion).
        Inputs
        Outputs
    %}
end
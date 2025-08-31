%{
    Helper code to test spheroidalDP.
    This code tests the code on a simple charge problem.
%}
function [soln, fluxsoln, truesoln, trueflux, sigma_vec,condK] = spheroidalDP_charge_problem(p, eta, ns, u0, target_distances, plt, neumann, interior)
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
        
        The idea is:
        (1) set up system of spheroids
        (2) put a few point charges in appropriate locations
              - Calculate the potential explicitly.
              - Determine potential on the boundary of each spheroid for BCs
        (3) Use MatVec self-evaluation (no X) to construct matrix operator D.
        (4) Solve (1/2*I + D + Completion)*sigma = f
              * Completion term will need some care, since it's multiple particles.
              For each particle it will be something like sum_i[ 1/|| x- c_i||
              integral(sigma_i) ]
        (5) Compare sigma with true point charge potential at targets that
            are 'target_distances' away from surface of each spheroid.
    %}

    %%% Backwards compatibility
    S_scale = 1;
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
    Xptch = reshape(repmat(reshape(c',3,1,[]),1,nc),3,[],1)';

    %%% Placement of point charges
    switch (interior)
        case true % Interior problem
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

            Xptch = Xptch + d; % Point charge locations, near centers of spheroids
        case false % Exterior problem
            if ns==3
                % 1/2 of the minor radius of the prolate spheroids
                scale = repmat(0.3.*pars.a .*sqrt(pars.u0.^2-1),nc,1);
                scale = scale(:);
                d = repmat(scale,1,3) .* (rand(size(Xptch))-.5);
            else
                d_z=0.5.*(rand(nc,1)-.5);
                d_xy=0.01.*(rand(nc,2)-0.5);
                d = [d_xy d_z];
            end

            Xptch = Xptch + d; % Point charge locations, near centers of spheroids
    end

    ptch = (2.*rand(ns*nc,1)-1); % Charge value, between -1 and 1

    if neumann || (~neumann & ~interior)
        % For the compatibility condition associated for Neumann problems
        ptch = ptch - sum(ptch)/(ns*nc);
    end
    
    if plt
        pars.plot(Xptch);
        title("spheroids and point charges")
    end

    Y = pars.get_X;
    [NrY,~]=pars.get_Norm_rot(p); % Necessary to account of the spheroid's rotation
    
    if ~neumann
        % Find boundary condition (potential induced by point charges)
        truesolnSurf=PtChargePotential(ptch,Xptch,Y);

        % Construct DL on-surface matrices
        DM = spheroidalMatVecKernel(pars,'DL',p);

        if interior % Interior problem
            % -1/2 I + D (should be well-conditioned!)
            K = -0.5*eye(ns*np) + DM;
        else % Exterior problem
            % 1/2 I + D + C
            CM = calculate_completion_term(pars, Y, ns, useS);
            K = .5*eye(ns*np) + DM + CM;
        end
    else % Neumann problem
        error("not implemented.");
        if FIRST_KIND_FLAG
            truesolnSurf=PtChargeFlux(ptch, Xptch, Y, NrY);
            K = spheroidalMatVecKernel(pars,'DP',p);
        else
            truesolnSurf=PtChargeFlux(ptch, Xptch,Y,NrY);

            % Construct SP on-surface matrices
            SPM= spheroidalMatVecKernel(pars,'SP',p);
            
            % Final operator: -1/2 I + SPM
            K = -.5*eye(ns*np) + SPM;
        end
    end

    condK = cond(K);
    
    sigma_vec=K\truesolnSurf;
    sigma = reshape(sigma_vec,np,[],ns);
    pars.sigma = sigma;
    
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
            soln = spheroidalMatVec(pars,'DL',Xeval);
        else % Exterior problem
            if useS
                C = S_scale.*spheroidalMatVec(pars,'SL',Xeval);
            else
                C=zeros(size(Xeval,1),1);
                for i=1:ns
                
                    if ~pars.oblate(i)
                        Yi = prolate_spheroid_shape(p,pars.u0(i),pars.a(i),'cart');
                    else
                        Yi = oblate_spheroid_shape(p,pars.u0(i),pars.a(i),'cart');
                    end
                    ci=c(i,:);
                    
                    % Completion matrix operator: sum of 1/||x-c_i|| * integral of each sigma over
                    % respective spheroid surface
                    Sns = SurfaceSph(Yi(:));
                    sigmaSurfInt = integrateOverS(Sns,sigma(:,:,i));
                
                    RXeval = vecnorm(Xeval-repmat(ci,size(Xeval,1),1),2,2);
                    
                    C = C + sigmaSurfInt./RXeval;
                end
            end
            soln = spheroidalMatVec(pars,'DL',Xeval) + C;
        end
    else
        if FIRST_KIND_FLAG
            soln = spheroidalMatVec(pars,'DL',Xeval);
        else
            soln = spheroidalMatVec(pars,'SL',Xeval);
        end
    end

    %%% Finally, we compare our computed solutions to the true solutions.
    truesoln = PtChargePotential(ptch,Xptch,Xeval);
    fprintf('p=%d: DL comparison = %.6e\n',p, norm(truesoln - soln) ./ norm(soln)); 

    %%% Compare flux solution evaluated away from the surface
    trueflux=PtChargeFlux(ptch, Xptch, Xeval, NrY);
    fluxsoln = spheroidalMatVec(pars, 'DP', Xeval, NrY);

    % For exterior Dirichlet problems, we must account for the completion 
    % term.
    if ~interior && ~neumann % Exterior Dirichlet
        completion_flux = zeros(size(Xeval,1), 1);
        c = pars.centers; % Get spheroid centers
    
        for i=1:ns
            if ~pars.oblate(i)
                Yi = prolate_spheroid_shape(p, pars.u0(i), pars.a(i), 'cart');
            else
                Yi = oblate_spheroid_shape(p, pars.u0(i), pars.a(i), 'cart');
            end
            Sns = SurfaceSph(Yi(:));
    
            % Calculate the total charge Qi on the i-th spheroid
            Qi = integrateOverS(Sns, sigma(:,:,i));
            
            % Calculate the flux contribution from this charge at all target points
            ci = c(i,:);
            r_vec = Xeval - ci; % Vectors from center to each target point
            r_norm = vecnorm(r_vec, 2, 2); % Distances ||x - c_i||
            
            % ((x-ci)/||x-ci||^3) . n_x
            E_dot_n = dot(r_vec, NrY, 2) ./ (r_norm.^3);
            
            % Add the contribution from this spheroid's completion charge
            completion_flux = completion_flux + Qi * E_dot_n;
        end
        
        % Add the completion flux to the total computed flux (accounting 
        % for the sign change due to the gradient).
        fluxsoln = fluxsoln - completion_flux;
    end

    fprintf('p=%d: DP comparison -- relative error = %.6e\n',p, norm(trueflux-fluxsoln) / norm(trueflux));
end

function pcp = PtChargePotential(ptch,Xptch,Y)
    %{
        Calculates the electric potential due to a set of point charges.
    %}
    M=length(ptch);
    np=size(Y,1);
    Rptch=zeros(np,M);

    % \|x - y\|
    for i=1:M
        Rptch(:,i)=sqrt(sum((Xptch(i,:)-Y).^2,2));
    end

    % \sum_j q_j/\|x - y\|
    pcp=1./(4*pi*Rptch)*ptch; 
end

function flux=PtChargeFlux(ptch,Xptch,Y,NrY)
    M=length(ptch);
    np=length(Y);
    Rptch=zeros(np,M);
    EdotN=zeros(np,M);
    for i=1:M
        % \|x - y\|
        r=Xptch(i,:)-Y;
        Rptch(:,i)=sqrt(sum(r.^2,2));
        
        EdotN(:,i)=dot(NrY,r./Rptch(:,i),2);
    end
    flux=EdotN./(4*pi*Rptch.^2)*ptch;
end

function unit_vectors = generate_random_unit_vec(num_vectors, vector_length)
    unit_vectors = zeros(num_vectors, vector_length);
    for i=1:num_vectors
        random_vector = randn(1, vector_length);
        unit_vectors(i,:) = random_vector / norm(random_vector);
    end
end

function CM = calculate_completion_term(pars, Y, ns, useS)
    p = pars.p; np = 2*p*(p + 1);
    if useS
        CM = S_scale.*spheroidalMatVecKernel(pars,'SL',p);
    else
        CM = [];
        [~, gwt]=g_grid(p+1);
        wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
        wt = wt(:)';
        for i=1:ns
            if ~pars.oblate(i)
                Yi = prolate_spheroid_shape(p,pars.u0(i),pars.a(i),'cart');
            else
                Yi = oblate_spheroid_shape(p,pars.u0(i),pars.a(i),'cart');
            end
            ci = pars.centers(i,:);
            
            % Completion matrix operator: sum of 1/||x-c_i|| * integral of each sigma over
            % respective spheroid surface
            Sns = SurfaceSph(Yi(:));
            W = Sns.geoProp.W;
        
            RY = vecnorm(Y-repmat(ci,np*ns,1),2,2);
            CM = [CM (1./RY)*(W'.*wt)];
        end
    end
end
classdef TEST_stokes_bie_problem < matlab.unittest.TestCase
    properties
    end

    methods (TestClassSetup)
        function setup(testCase)
            clc();

            addpath(genpath('../../../spheroidal'));
            addpath(genpath('../../../support'))
        end
    end

    methods (Test)
        function testInteriorDirichletProblemOneSpheroid(testCase)
            rng(42);
            p = 16;
            eta = 10;
            ns = 1;
            u0 = 2/sqrt(3);
            target_distances = 0.2; % Distance from the surface to evaluate the potential
            plt = false;
            neumann = false;
            interior = true;
            
            [soln, ~, truesoln, ~, sigma_vec, condK] =  stokes_bie_problem(p, eta, ns, u0, target_distances, plt, neumann, interior);
        end

        function testInteriorDirichletProblemTwoSpheroids(testCase)
            rng(2024);
            p = 16;
            eta = 10;
            ns = 2;
            u0 = [5 2/sqrt(3)];
            target_distances = 0.2*ones(1, ns);
            plt = false;
            neumann = false;
            interior = true;
            
            [soln, ~, truesoln, ~, sigma_vec, condK] =  stokes_bie_problem(p, eta, ns, u0, target_distances, plt, neumann, interior);
        end

        function testInteriorDirichletProblemThreeSpheroids(testCase)
            p = 16;
            eta = 30;
            ns = 3;
            u0 = [1.5 2 1.7];
            target_distances = 0.2*ones(1, ns);
            plt = true;
            neumann = false;
            interior = true;
            
            [soln, ~, truesoln, ~, sigma_vec, condK] =  stokes_bie_problem(p, eta, ns, u0, target_distances, plt, neumann, interior);
        end

        function testExteriorDirichletProblemOneSpheroid(testCase)
            rng(42);
            p = 16;
            eta = 10;
            ns = 1;
            u0 = 2/sqrt(3);
            target_distances = 1; % Distance from the surface to evaluate the potential
            plt = true;
            neumann = false;
            interior = false;
            
            [soln, ~, truesoln, ~, sigma_vec, condK] =  stokes_bie_problem(p, eta, ns, u0, target_distances, plt, neumann, interior);
        end

        function testExteriorDirichletProblemTwoSpheroids(testCase)
            rng(42);
            p = 16;
            eta = 10;
            ns = 2;
            u0 = [21/sqrt(41) 1.8];
            target_distances = 1*ones(1, ns);
            plt = false;
            neumann = false;
            interior = false;
            
            [soln, ~, truesoln, ~, sigma_vec, condK] =  stokes_bie_problem(p, eta, ns, u0, target_distances, plt, neumann, interior);
        end

        function testInteriorNeumannProblemOneSpheroid(testCase)
            rng(42);
            p = 16;
            eta = 10;
            ns = 1;
            u0 = 3;
            target_distances = 1e-1; % Distance from the surface to evaluate the potential
            plt = false;
            neumann = true;
            interior = true;
            
            [soln, fluxsoln, truesoln, truefluxsoln, sigma_vec, condK] =  stokes_bie_problem(p, eta, ns, u0, target_distances, plt, neumann, interior);
        end

        function testInteriorNeumannProblemTwoSpheroids(testCase)
            rng(42);
            p = 16;
            eta = 10;
            ns = 2;
            u0 = [500 500];
            target_distances = 1e-2*ones(ns, 1); % Distance from the surface to evaluate the potential
            plt = false;
            neumann = true;
            interior = true;
            
            [soln, fluxsoln, truesoln, truefluxsoln, sigma_vec, condK] =  stokes_bie_problem(p, eta, ns, u0, target_distances, plt, neumann, interior);
        end

        function testExteriorNeumannProblemOneSpheroid(testCase)
            rng(42);
            p = 16;
            eta = 10;
            ns = 1;
            % u0 = 11;
            % u0 = 11/sqrt(21); % AR = 1.1
            % u0 = 2/sqrt(3); % AR = 2
            % u0 = 4/sqrt(15); % AR = 4
            u0 = 500;
            target_distances = 1e-5; % Distance from the surface to evaluate the potential
            plt = false;
            neumann = true;
            interior = false;
            
            [soln, fluxsoln, truesoln, truefluxsoln, sigma_vec, condK] =  stokes_bie_problem(p, eta, ns, u0, target_distances, plt, neumann, interior);
        end

        function testExteriorNeumannProblemTwoSpheroids(testCase)
            rng(42);
            p = 16;
            eta = 5000;
            ns = 2;
            u0 = [5 5];
            target_distances = 0.5*ones(1, ns);
            plt = false;
            neumann = true;
            interior = false;
            
            [soln, fluxsoln, truesoln, truefluxsoln, sigma_vec, condK] =  stokes_bie_problem(p, eta, ns, u0, target_distances, plt, neumann, interior);
        end

        function testExteriorNeumannProblemThreeSpheroids(testCase)
            rng(42);
            p = 16;
            eta = 5000;
            ns = 3;
            u0 = [2/sqrt(3) 5 3];
            target_distances = 3*ones(1, ns);
            plt = false;
            neumann = true;
            interior = false;
            
            [soln, fluxsoln, truesoln, truefluxsoln, sigma_vec, condK] =  stokes_bie_problem(p, eta, ns, u0, target_distances, plt, neumann, interior);
        end

        %% Test helper functions
        function testStokesletVelocityRotationalCovariance(testCase)
            %{
                In an unrotated setting, one can calculate the velocity due
                to a number of point forces. Then, if one moves the entire 
                system (including the point forces), the velocity should 
                only differ by a rotation.
            %}
            rng(42);
            p = 16;
            np = 2*p*(p+1);
            params = SpheroidalParameters();
            params.matvec_eta = 10;
            params.u0 = 2/sqrt(3);
            params.a = 1/params.u0;
            params.oblate = false; 
            params.centers = [0 0 0];
            params.isReal = true;

            params.sigma = ones(np, 1); % A hack to update p...

            %% Unrotated system
            params.thetas = 0;
            params.phis = 0;
            [X_src_unrotated, ~] = params.get_X;

            % Place the point forces along the poles of the spheroid in the
            % exterior
            F_pos_vec = [0 0 2; 0 0 -2];
            F_vec = randn(2, 3); % Random point forces 

            surf_velocity_unrotated = stokeslet_velocity(F_vec, F_pos_vec, X_src_unrotated);

            %% Rotated system
            params.thetas = pi/6;
            params.phis = 0;
            [X_src_rotated, ~] = params.get_X;

            thetai = params.thetas;
            phii = params.phis;
            Riy=[cos(thetai) 0 sin(thetai); 0 1 0; -sin(thetai) 0 cos(thetai)];
            Riz=[cos(phii) -sin(phii) 0; sin(phii) cos(phii) 0; 0 0 1];
            Ri=Riz*Riy;

            surf_velocity_rotated = stokeslet_velocity(F_vec*Ri', F_pos_vec*Ri', X_src_rotated);

            %% Comparison
            abs_err = norm(surf_velocity_rotated - surf_velocity_unrotated*Ri');
            testCase.verifyLessThan(abs_err, 1e-15, ...
                'Velocity generated by point forces does not satisfy expected rotational relationship.');
        end

        function testStokesletVelocityRotationalCovarianceMultipleSpheroids(testCase)
            %{
                Like the other tests, rotating the entire global coordinate
                system should only affect the true velocity by the global
                rotation.
            %}
            rng(42);
            p = 16;
            np = 2*p*(p+1);
            ns = 2;
            params = SpheroidalParameters();
            params.matvec_eta = 10;
            params.u0 = [2/sqrt(3) 2/sqrt(3)];
            params.a = 1./params.u0;
            params.oblate = [0 0]; 
            params.centers = [0 0 0; 5 0 0];
            params.isReal = true;
            params.sigma = ones(np, 1);

            %% Unrotated system
            params.thetas = [0 0];
            params.phis = [0 0];

            Ri=zeros(3,3,ns);
            for i=1:ns
                thetai=params.thetas(i);
                phii=params.phis(i);
                Riy=[cos(thetai) 0 sin(thetai); 0 1 0; -sin(thetai) 0 cos(thetai)];
                Riz=[cos(phii) -sin(phii) 0; sin(phii) cos(phii) 0; 0 0 1];
                Ri(:,:,i)=Riz*Riy;
            end
            params.Rmat=Ri;

            [X_src_unrotated, ~] = params.get_X;

            % Place the point forces along the poles of the spheroid in the
            % exterior
            F_pos_vec = [
                0 0 2; 0 0 -2;
                5 0 2; 5 0 -2;
                ];
            F_vec = randn(2*ns, 3); % Random point forces 

            surf_velocity_unrotated = stokeslet_velocity(F_vec, F_pos_vec, X_src_unrotated);

            %% Global rotation
            fprintf("Global rotation flag enabled...\n");
            theta_g = -pi/4;
            phi_g = 0;
            Riy = [cos(theta_g) 0 sin(theta_g); 0 1 0; -sin(theta_g) 0 cos(theta_g)];
            Riz = [cos(phi_g) -sin(phi_g) 0; sin(phi_g) cos(phi_g) 0; 0 0 1];
            R_global = Riz*Riy;
    
            fprintf("Rotating centers...\n");
            params.centers = params.centers * R_global';
    
            fprintf("Re-orienting each spheroid...\n");
            for i = 1:ns
                Ri(:,:,i) = R_global * Ri(:,:,i);
            end
            params.Rmat = Ri;
    
            fprintf("Done re-orienting coordinate system.\n");
          
            %% Rotated system
            [X_src_rotated, ~] = params.get_X;

            surf_velocity_rotated = stokeslet_velocity(F_vec*R_global', F_pos_vec*R_global', X_src_rotated);

            %% Comparison
            abs_err = norm(surf_velocity_rotated - surf_velocity_unrotated*R_global');
            testCase.verifyLessThan(abs_err, 1e-15, ...
                'Velocity generated by point forces (for multiple spheroids) does not satisfy expected rotational relationship.');
        end

        function testStressletTractionRotationalCovariance(testCase)
            rng(42);
            p = 16;
            np = 2*p*(p+1);
            params = SpheroidalParameters();
            params.matvec_eta = 10;
            params.u0 = 2/sqrt(3);
            params.a = 1/params.u0;
            params.oblate = false; 
            params.centers = [0 0 0];
            params.isReal = true;

            params.sigma = ones(np, 1); % A hack to update p...

            %% Unrotated system
            params.thetas = 0;
            params.phis = 0;
            [X_src_unrotated, ~] = params.get_X;
            N_src_unrotated = params.get_Norm_rot(p);

            % Place the point forces along the poles of the spheroid in the
            % exterior
            F_pos_vec = [0 0 2; 0 0 -2];
            F_vec = randn(2, 3); % Random point forces 

            surf_traction_unrotated = stresslet_traction(F_vec, F_pos_vec, X_src_unrotated, N_src_unrotated);
            
            %% Rotated system
            params.thetas = pi/6;
            params.phis = 0;
            [X_src_rotated, ~] = params.get_X;
            N_src_rotated = params.get_Norm_rot(p);

            thetai = params.thetas;
            phii = params.phis;
            Riy=[cos(thetai) 0 sin(thetai); 0 1 0; -sin(thetai) 0 cos(thetai)];
            Riz=[cos(phii) -sin(phii) 0; sin(phii) cos(phii) 0; 0 0 1];
            Ri=Riz*Riy;

            surf_traction_rotated = stresslet_traction(F_vec*Ri', F_pos_vec*Ri', X_src_rotated, N_src_rotated);

            %% Comparison
            abs_err = norm(surf_traction_rotated - surf_traction_unrotated*Ri');
            testCase.verifyLessThan(abs_err, 1e-15, ...
                'Traction generated by point forces does not satisfy expected rotational relationship.');
        end
        
        %% Test completion term
        function testRBMCompletionTerm(testCase)
            p = 16;
            np = 2*p*(p+1);
            ns = 1;
            params = SpheroidalParameters();
            params.matvec_eta = 10;
            params.u0 = 2/sqrt(3);
            params.a = 1/params.u0;
            params.oblate = false; 
            params.centers = [0 0 0];
            params.isReal = true;

            params.sigma = ones(np, 1); % A hack to update p...
            params.thetas = 0;
            params.phis = 0;
            [X_src, ~] = params.get_X;

            %% Quadrature weights
            Sns = SurfaceSph(X_src);
            [~, gwt_gl] = g_grid(p + 1);
            wt_gl = pi/p * repmat(gwt_gl', 2*p, 1) ./ sin(gl_grid(p));
            wt_gl = wt_gl(:);
            W_src = Sns.geoProp.W .* wt_gl;
            Wv = repmat(W_src,1,3)'; Wv=Wv(:);

            rbm_completion_term = RBM_completion(params, X_src, ns);
            [C, B, D, L] = Build_SpheroidalAuxMats(W_src, X_src, params.centers, np, ns);
            L = full(L);

            %% Permute L since it is in interleaved format
            prm = zeros(1,3*np); 
            prm(1:np) = 1:3:3*np; prm(np+1:2*np) = 2:3:3*np; prm(2*np+1:3*np)=3:3:3*np;
            L = L(prm, prm);
            
            %% Compare results
            abs_err = norm(rbm_completion_term - L);
        end
    end
end

%% Helper functions
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
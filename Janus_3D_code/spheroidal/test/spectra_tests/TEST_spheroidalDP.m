%{
    Test code for spheroidalDP.m.
%}


classdef TEST_spheroidalDP < matlab.unittest.TestCase
    properties
        u0_prolate = 2/sqrt(3);
        u0_oblate = 2/sqrt(3);
        a_prolate;
        a_oblate;
        p_max = 32; % Maximum order for convergence tests
        conv_tol = 1e-5;
        tol = 1e-6;

        % Consistency tolerance
        consistency_tol = 1e-14;

        % Tolerance for gradient checks
        gradient_check_tol = 1e-6;
        fd_eps = 1e-5;

        % Non-trivial density function: chosen so that it is smooth
        % and does not allow the convergence tests to hit machine precision 
        % (i.e. 1e-16) too quickly for low orders.
        density_func = @(u,v) sin(u).^4 .* cos(v);
    end

    methods (TestClassSetup)
        function setup(testCase)
            clc();
            testCase.a_prolate = 1/testCase.u0_prolate;
            testCase.a_oblate = 1/sqrt(1 + testCase.u0_oblate^2);
        end
    end

    methods (Test)
        %%% Convergence tests
        function testProlateConvergenceOrderOffSurface(testCase)
            p_orders = [8, 16, 24];

            % Setup reference solution with provided maximum order.
            params_ref = SpheroidalParameters;
            params_ref.u0 = testCase.u0_prolate;
            params_ref.a = testCase.a_prolate;
            params_ref.oblate = false;
            [u_ref, v_ref] = gl_grid(testCase.p_max);
            params_ref.sigma = testCase.density_func(u_ref, v_ref);
            params_ref.get_shc();

            % Off-surface target points
            nu_trg = get_norm_vecs(testCase.p_max, params_ref.u0, params_ref.oblate);
            X_trg = prolate_spheroid_shape(testCase.p_max, params_ref.u0, params_ref.a) + nu_trg;

            DP_ref = spheroidalDP(params_ref, X_trg, nu_trg);

            % Loop over lower orders
            errors = zeros(size(p_orders));
            for i = 1:length(p_orders)
                p = p_orders(i);
                params_p = SpheroidalParameters;
                params_p.u0 = testCase.u0_prolate;
                params_p.a = testCase.a_prolate;
                params_p.oblate = false;
                [u_p, v_p] = gl_grid(p);
                params_p.sigma = testCase.density_func(u_p, v_p);
                params_p.get_shc();

                DP_p = spheroidalDP(params_p, X_trg, nu_trg);

                errors(i) = norm(DP_p - DP_ref, inf) / norm(DP_ref, inf);
            end

            % Ensures that the errors are at least decreasing by a factor
            % of 10 everytime we increase the order
            testCase.verifyTrue(all(errors(2:end) < errors(1:end-1) / 10), ...
                'Spectral convergence should be observed as the order p is increased.');
            testCase.verifyLessThan(errors(end), testCase.conv_tol, ...
                'Error for prolate off-surface should be below tolerance.');
        end

        function testOblateConvergenceOrderOffSurface(testCase)
            p_orders = [8, 16, 24];

            % Setup reference solution with provided maximum order.
            params_ref = SpheroidalParameters;
            params_ref.u0 = testCase.u0_oblate;
            params_ref.a = testCase.a_oblate;
            params_ref.oblate = true;
            [u_ref, v_ref] = gl_grid(testCase.p_max);
            params_ref.sigma = testCase.density_func(u_ref, v_ref);
            params_ref.get_shc();

            % Off-surface target points
            nu_trg = get_norm_vecs(testCase.p_max, params_ref.u0, params_ref.oblate);
            X_trg = prolate_spheroid_shape(testCase.p_max, params_ref.u0, params_ref.a) + nu_trg;

            DP_ref = spheroidalDP(params_ref, X_trg, nu_trg);

            % Loop over lower orders
            errors = zeros(size(p_orders));
            for i = 1:length(p_orders)
                p = p_orders(i);
                params_p = SpheroidalParameters;
                params_p.u0 = testCase.u0_oblate;
                params_p.a = testCase.a_oblate;
                params_p.oblate = true;
                [u_p, v_p] = gl_grid(p);
                params_p.sigma = testCase.density_func(u_p, v_p);
                params_p.get_shc();

                DP_p = spheroidalDP(params_p, X_trg, nu_trg);

                errors(i) = norm(DP_p - DP_ref, inf) / norm(DP_ref, inf);
            end

            testCase.verifyTrue(all(errors(2:end) < errors(1:end-1) / 10), ...
                'Spectral convergence should be observed as the order p is increased.');
            testCase.verifyLessThan(errors(end), testCase.conv_tol, ...
                'Error for oblate off-surface should be below tolerance.');
        end

        function testProlateConvergenceOnSurface(testCase)
            p_orders = [8, 12, 16, 20, 24, 28, 32]; % Finer grid needed for on-surface evaluations

            params_ref = SpheroidalParameters;
            params_ref.u0 = testCase.u0_prolate;
            params_ref.a = testCase.a_prolate;
            params_ref.oblate = false;
            [u_ref, v_ref] = gl_grid(testCase.p_max);
            params_ref.sigma = testCase.density_func(u_ref, v_ref);
            params_ref.get_shc();

            % By default, we have outward normal vectors
            DP_ref = spheroidalDP(params_ref);

            errors = zeros(size(p_orders));
            for i = 1:length(p_orders)
                p = p_orders(i);
                params_p = SpheroidalParameters;
                params_p.u0 = testCase.u0_prolate;
                params_p.a = testCase.a_prolate;
                params_p.oblate = false;
                [u_p, v_p] = gl_grid(p);
                params_p.sigma = testCase.density_func(u_p, v_p);
                params_p.get_shc();

                DP_p_full = spheroidalDP(params_p);

                % We use spherical transforms for spheroidal transforms, so this should be OK.
                DP_p_interp = interpsh(DP_p_full, testCase.p_max);
                errors(i) = norm(DP_p_interp - DP_ref, inf) / norm(DP_ref, inf);
            end

            testCase.verifyLessThan(errors(end), testCase.conv_tol, ...
                'Error for prolate on-surface should be below tolerance.');
        end

        function testOblateConvergenceOnSurface(testCase)
            p_orders = [8, 12, 16, 20, 24, 28, 32]; % Finer grid needed for on-surface evaluations

            params_ref = SpheroidalParameters;
            params_ref.u0 = testCase.u0_oblate;
            params_ref.a = testCase.a_oblate;
            params_ref.oblate = true;
            [u_ref, v_ref] = gl_grid(testCase.p_max);
            params_ref.sigma = testCase.density_func(u_ref, v_ref);
            params_ref.get_shc();

            % By default, we have outward normal vectors
            DP_ref = spheroidalDP(params_ref);

            errors = zeros(size(p_orders));
            for i = 1:length(p_orders)
                p = p_orders(i);
                params_p = SpheroidalParameters;
                params_p.u0 = testCase.u0_oblate;
                params_p.a = testCase.a_oblate;
                params_p.oblate = true;
                [u_p, v_p] = gl_grid(p);
                params_p.sigma = testCase.density_func(u_p, v_p);
                params_p.get_shc();

                DP_p_full = spheroidalDP(params_p);

                DP_p_interp = interpsh(DP_p_full, testCase.p_max);
                errors(i) = norm(DP_p_interp - DP_ref, inf) / norm(DP_ref, inf);
            end

            testCase.verifyLessThan(errors(end), testCase.conv_tol, ...
                'Error for oblate on-surface should be below tolerance.');
        end

        %%% Consistency checks
        function testCompareProlateOnSurfaceImplementations(testCase)
            % Compares all method to calculate on the surface: they should all
            % match each other.
            p = 16;

            params = SpheroidalParameters;
            params.isReal = true;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = false;
            [u_p, v_p] = gl_grid(p);
            params.sigma = testCase.density_func(u_p, v_p);
            params.get_shc();

            % These normal vectors are in Cartesian coordinates, but in
            % spheroidal coordinates, it is just e_u.
            nu_src = get_norm_vecs(p, params.u0, params.oblate);
            [X_src, ~] = params.get_X(); 

            % Method 1: use off-surface code
            DP_method1 = spheroidalDP(params, X_src, nu_src);

            % Method 2: handle []
            DP_method2 = spheroidalDP(params, [], nu_src);

            % Method 3: no inputs
            DP_method3 = spheroidalDP(params);

            testCase.verifyEqual(DP_method1, DP_method2, 'AbsTol', testCase.consistency_tol, ...
                'Off-surface code and [] handling does not match each other.');

            testCase.verifyEqual(DP_method2, DP_method3, 'AbsTol', testCase.consistency_tol, ...
                'On-surface code and [] handling does not match each other.');
        end

        function testCompareOblateOnSurfaceImplementations(testCase)
            % Compares all method to calculate on the surface: they should all
            % match each other.
            p = 16;

            params = SpheroidalParameters;
            params.u0 = testCase.u0_oblate;
            params.a = testCase.a_oblate;
            params.oblate = true;
            [u_p, v_p] = gl_grid(p);
            params.sigma = testCase.density_func(u_p, v_p);
            params.get_shc();

            % These normal vectors are in Cartesian coordinates, but in
            % spheroidal coordinates, it is just e_u.
            nu_src = get_norm_vecs(p, params.u0, params.oblate);
            [X_src, ~] = params.get_X(); 

            % Method 1: use off-surface code
            DP_method1 = spheroidalDP(params, X_src, nu_src);

            % Method 2: handle []
            DP_method2 = spheroidalDP(params, [], nu_src);

            % Method 3: no inputs
            DP_method3 = spheroidalDP(params);

            testCase.verifyEqual(DP_method1, DP_method2, 'AbsTol', testCase.consistency_tol, ...
                'Off-surface code and [] handling does not match each other.');

            testCase.verifyEqual(DP_method2, DP_method3, 'AbsTol', testCase.consistency_tol, ...
                'On-surface code and [] handling does not match each other.');
        end

        %%% Gradient checks
        function testGradientCheckProlateOffSurface(testCase)
            p = 16;
            eps = testCase.fd_eps;

            params = SpheroidalParameters;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = false;
            [u_p, v_p] = gl_grid(p);
            params.sigma = testCase.density_func(u_p, v_p);
            params.get_shc();

            nu_trg = get_norm_vecs(p, params.u0, params.oblate);
            X_trg = oblate_spheroid_shape(p, params.u0, params.a) + nu_trg;

            DP_spectral = spheroidalDP(params, X_trg, nu_trg);

            X_plus = X_trg + eps * nu_trg;
            X_minus = X_trg - eps * nu_trg;

            DL_plus = spheroidalDL(params, X_plus);
            DL_minus = spheroidalDL(params, X_minus);
            DP_fd = (DL_plus - DL_minus) / (2 * eps);

            rel_err = norm(DP_spectral - DP_fd, inf) / norm(DP_spectral, inf);
            testCase.verifyLessThan(rel_err, testCase.gradient_check_tol, ...
                'Gradient check failed for prolate case: should be below tolerance.');
        end

        function testGradientCheckOblateOffSurface(testCase)
            p = 16;
            eps = testCase.fd_eps;

            params = SpheroidalParameters;
            params.u0 = testCase.u0_oblate;
            params.a = testCase.a_oblate;
            params.oblate = true;
            [u_p, v_p] = gl_grid(p);
            params.sigma = testCase.density_func(u_p, v_p);
            params.get_shc();

            nu_trg = get_norm_vecs(p, params.u0, params.oblate);
            X_trg = oblate_spheroid_shape(p, params.u0, params.a) + nu_trg;

            DP_spectral = spheroidalDP(params, X_trg, nu_trg);

            DL_plus = spheroidalDL(params, X_trg + eps * nu_trg);
            DL_minus = spheroidalDL(params, X_trg - eps * nu_trg);
            DP_fd = (DL_plus - DL_minus) / (2 * eps);

            rel_err = norm(DP_spectral - DP_fd, inf) / norm(DP_spectral, inf);
            testCase.verifyLessThan(rel_err, testCase.gradient_check_tol, ...
                'Gradient check failed for oblate case: should be below tolerance.');
        end

        %%% Spectral coefficient tests
        function testProlateSpectralCoefficients(testCase)
            %{
                Use that the spheroidal harmonics Y_n^m form an orthogonal
                family to test the spectral coefficients.
            %}

        end

        %%% Kernel_Eval check
        function testProlateWithKernelEvalOffSurface(testCase)
            %{
                Verify that off-surface evaluations coincide with the
                implementation in Kernel_Eval for the prolate case.
            %}
            p = 16;
            
            params = SpheroidalParameters;
            params.isReal = true;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = false;
            [u_p, v_p] = gl_grid(p);
            params.sigma = testCase.density_func(u_p, v_p);
            params.get_shc();

            nu_trg = get_norm_vecs(p, params.u0, params.oblate);
            X_trg = prolate_spheroid_shape(p, params.u0, params.a) + nu_trg;

            DP_spectral = spheroidalDP(params, X_trg, nu_trg);

            % Get source geometry and weights
            [X_src_orig, ~] = params.get_X(); % Cartesian coordinates of source points
            N_src_orig = params.get_Norm(p, 1);
            
            % Quadrature weights for Kernel_Eval
            Sns = SurfaceSph(X_src_orig);
            [~, gwt_gl] = g_grid(p + 1);
            wt_gl = pi/p * repmat(gwt_gl', 2*p, 1) ./ sin(gl_grid(p));
            wt_gl = wt_gl(:);
            W_src_orig = Sns.geoProp.W .* wt_gl;

            pot = 'dDL_L_3D';
            KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-12,2,400,1);
            KEparams.dim = 3;
            KEparams.X = X_src_orig;
            KEparams.W2 = W_src_orig.';
            KEparams.nor = N_src_orig;
            KEparams.targnor = nu_trg;

            dDL_mat = Kernel_Eval(X_trg, X_src_orig, KEparams);
            DP_kernel_eval = dDL_mat * params.sigma;

            rel_err = norm(DP_spectral - DP_kernel_eval, inf) / norm(DP_spectral, inf);
            testCase.verifyLessThan(rel_err, testCase.gradient_check_tol, ...
                'spheroidalDP vs Kernel_Eval for prolate off-surface failed.');
        end

        function testOblateWithKernelEvalOffSurface(testCase)
            %{
                Verify that off-surface evaluations coincide with the
                implementation in Kernel_Eval for the oblate case.

                Note that the oblate case requires higher order (i.e. 'p') 
                for the desired convergence.
            %}
            p = 16;
            
            params = SpheroidalParameters;
            params.isReal = true;
            params.u0 = testCase.u0_oblate;
            params.a = testCase.a_oblate;
            params.oblate = true;
            [u_p, v_p] = gl_grid(p);
            params.sigma = testCase.density_func(u_p, v_p);
            params.get_shc();

            nu_trg = get_norm_vecs(p, params.u0, params.oblate);
            X_trg = oblate_spheroid_shape(p, params.u0, params.a) + 2*nu_trg;

            DP_spectral = spheroidalDP(params, X_trg, nu_trg);

            % Get source geometry and weights
            [X_src_orig, ~] = params.get_X(); % Cartesian coordinates of source points
            N_src_orig = params.get_Norm(p, 1);
            
            % Quadrature weights for Kernel_Eval
            Sns = SurfaceSph(X_src_orig);
            [~, gwt_gl] = g_grid(p + 1);
            wt_gl = pi/p * repmat(gwt_gl', 2*p, 1) ./ sin(gl_grid(p));
            wt_gl = wt_gl(:);
            W_src_orig = Sns.geoProp.W .* wt_gl;

            pot = 'dDL_L_3D';
            KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-12,2,400,1);
            KEparams.dim = 3;
            KEparams.X = X_src_orig;
            KEparams.W2 = W_src_orig.';
            KEparams.nor = N_src_orig;
            KEparams.targnor = nu_trg;

            dDL_mat = Kernel_Eval(X_trg, X_src_orig, KEparams);
            DP_kernel_eval = dDL_mat * params.sigma;

            rel_err = norm(DP_spectral - DP_kernel_eval, inf) / norm(DP_spectral, inf);
            testCase.verifyLessThan(rel_err, testCase.gradient_check_tol, ...
                'spheroidalDP vs Kernel_Eval for prolate off-surface failed.');
        end

        function testRandomNormalVectors(testCase)
            %{
                Verifies that spheroidalDP still functions with random
                normal vectors (i.e. the choice of normal vectors should
                not matter). Compares with Kernel_Eval and handles both the
                prolate and oblate cases.
            %}
            p = 16;
            
            %%% Prolate
            params = SpheroidalParameters;
            params.isReal = true;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = false;
            [u_p, v_p] = gl_grid(p);
            params.sigma = testCase.density_func(u_p, v_p);
            params.get_shc();

            nu_trg = get_norm_vecs(p, params.u0, params.oblate);
            X_trg = oblate_spheroid_shape(p, params.u0, params.a) + 2*nu_trg;
            nu_trg = nu_trg + 2*rand(size(nu_trg)); % Random perturbation

            DP_spectral = spheroidalDP(params, X_trg, nu_trg);

            % Get source geometry and weights
            [X_src_orig, ~] = params.get_X(); % Cartesian coordinates of source points
            N_src_orig = params.get_Norm(p, 1);
            
            % Quadrature weights for Kernel_Eval
            Sns = SurfaceSph(X_src_orig);
            [~, gwt_gl] = g_grid(p + 1);
            wt_gl = pi/p * repmat(gwt_gl', 2*p, 1) ./ sin(gl_grid(p));
            wt_gl = wt_gl(:);
            W_src_orig = Sns.geoProp.W .* wt_gl;

            pot = 'dDL_L_3D';
            KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-12,2,400,1);
            KEparams.dim = 3;
            KEparams.X = X_src_orig;
            KEparams.W2 = W_src_orig.';
            KEparams.nor = N_src_orig;
            KEparams.targnor = nu_trg;

            dDL_mat = Kernel_Eval(X_trg, X_src_orig, KEparams);
            DP_kernel_eval = dDL_mat * params.sigma;

            rel_err = norm(DP_spectral - DP_kernel_eval, inf) / norm(DP_spectral, inf);
            testCase.verifyLessThan(rel_err, testCase.gradient_check_tol, ...
                'Randomized spheroidalDP vs Kernel_Eval for prolate off-surface failed.');

            %%% Oblate
            params = SpheroidalParameters;
            params.isReal = true;
            params.u0 = testCase.u0_oblate;
            params.a = testCase.a_oblate;
            params.oblate = true;
            [u_p, v_p] = gl_grid(p);
            params.sigma = testCase.density_func(u_p, v_p);
            params.get_shc();

            nu_trg = get_norm_vecs(p, params.u0, params.oblate);
            X_trg = oblate_spheroid_shape(p, params.u0, params.a) + 2*nu_trg;
            nu_trg = nu_trg + 2*rand(size(nu_trg)); % Random perturbation

            DP_spectral = spheroidalDP(params, X_trg, nu_trg);

            % Get source geometry and weights
            [X_src_orig, ~] = params.get_X(); % Cartesian coordinates of source points
            N_src_orig = params.get_Norm(p, 1);
            
            % Quadrature weights for Kernel_Eval
            Sns = SurfaceSph(X_src_orig);
            [~, gwt_gl] = g_grid(p + 1);
            wt_gl = pi/p * repmat(gwt_gl', 2*p, 1) ./ sin(gl_grid(p));
            wt_gl = wt_gl(:);
            W_src_orig = Sns.geoProp.W .* wt_gl;

            pot = 'dDL_L_3D';
            KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-12,2,400,1);
            KEparams.dim = 3;
            KEparams.X = X_src_orig;
            KEparams.W2 = W_src_orig.';
            KEparams.nor = N_src_orig;
            KEparams.targnor = nu_trg;

            dDL_mat = Kernel_Eval(X_trg, X_src_orig, KEparams);
            DP_kernel_eval = dDL_mat * params.sigma;

            rel_err = norm(DP_spectral - DP_kernel_eval, inf) / norm(DP_spectral, inf);
            testCase.verifyLessThan(rel_err, testCase.gradient_check_tol, ...
                'Randomized spheroidalDP vs Kernel_Eval for oblate off-surface failed.');
        end

        %%% On-surface checks
        function testNullspaceIsConstant(testCase)
            %{
                Verifies that the nullspace of the Laplace double-layer
                potential is indeed the space generated by the constant
                function.
            %}
            p = 16;
            np = 2*p*(p+1);
            
            %%% Prolate
            params = SpheroidalParameters;
            params.isReal = true;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = false;
            params.sigma = ones(np, 1, 1);
            params.get_shc();

            % On-surface
            nu_arbitrary = repmat(rand(1,3), np, 1, 1);
            DP_onsurface = spheroidalDP(params, [], nu_arbitrary);
            testCase.verifyLessThan(norm(DP_onsurface, inf), 1e-10, ...
                'Prolate: On-surface DP for constant density with arbitrary normal should be near zero.');

            % Off-surface
            nu_trg = get_norm_vecs(p, params.u0, params.oblate);
            X_trg = prolate_spheroid_shape(p, params.u0, params.a) + 2*nu_trg;
            DP_offsurface = spheroidalDP(params, X_trg, nu_trg);
            testCase.verifyLessThan(norm(DP_offsurface, inf), testCase.consistency_tol, ...
                'Prolate: Off-surface DP for constant density should be near zero.');

            %%% Oblate
            params = SpheroidalParameters;
            params.isReal = true;
            params.u0 = testCase.u0_oblate;
            params.a = testCase.a_oblate;
            params.oblate = true;
            params.sigma = ones(np, 1, 1);
            params.get_shc();

            % On-surface
            nu_arbitrary = repmat(rand(1,3), np, 1, 1);
            DP_onsurface = spheroidalDP(params, [], nu_arbitrary);
            testCase.verifyLessThan(norm(DP_onsurface, inf), 1e-10, ...
                'Oblate: On-surface DP for constant density with arbitrary normal should be near zero.');

            % Off-surface
            nu_trg = get_norm_vecs(p, params.u0, params.oblate);
            X_trg = oblate_spheroid_shape(p, params.u0, params.a) + 2*nu_trg;
            DP_offsurface = spheroidalDP(params, X_trg, nu_trg);
            testCase.verifyLessThan(norm(DP_offsurface, inf), testCase.consistency_tol, ...
                'Oblate: Off-surface DP for constant density should be near zero.');
        end
        
        %%% Near-surface tests
        function testJumpRelationProlateConvergenceTest(testCase)
            %{
                As we approach the surface of the spheroid, we should
                satisfy the jump relation. However, since the normal
                derivative of the Laplace DLP is continuous across the
                surface, we should just converge to the on-surface
                evaluation.

                For reference for this jump relation fact, see Hsiao and
                Wedland.
            %}
            p = 16;
            
            %%% Prolate
            params = SpheroidalParameters;
            params.isReal = true;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = false;
            [u_p, v_p] = gl_grid(p);
            params.sigma = testCase.density_func(u_p, v_p);
            params.get_shc();

            X_self = params.get_X();
            nu_self = get_norm_vecs(p, params.u0, params.oblate);
            DP_onsurface = spheroidalDP(params, [], nu_self);

            distances_from_surface = 10.^(-4:-1:-10);
            errors = zeros(size(distances_from_surface));
            for i = 1:length(distances_from_surface)
                DP_offsurface = spheroidalDP(params, X_self + nu_self .* distances_from_surface(i), nu_self);
                errors(i) = norm(DP_onsurface - DP_offsurface) / norm(DP_offsurface);
            end

            testCase.verifyLessThan(errors(2:end), errors(1:end-1) / 10, ...
                'Spectral convergence should be observed for the prolate case as the order p is increased.');
        end

        function testJumpRelationOblateConvergenceTest(testCase)
            p = 16;

            params = SpheroidalParameters;
            params.isReal = true;
            params.u0 = testCase.u0_oblate;
            params.a = testCase.a_oblate;
            params.oblate = true;
            [u_p, v_p] = gl_grid(p);
            params.sigma = testCase.density_func(u_p, v_p);
            params.get_shc();

            X_self = params.get_X();
            nu_self = get_norm_vecs(p, params.u0, params.oblate);
            DP_onsurface = spheroidalDP(params, [], nu_self);

            distances_from_surface = 10.^(-4:-1:-10);
            errors = zeros(size(distances_from_surface));
            for i = 1:length(distances_from_surface)
                DP_offsurface = spheroidalDP(params, X_self + nu_self .* distances_from_surface(i), nu_self);
                errors(i) = norm(DP_onsurface - DP_offsurface) / norm(DP_offsurface);
            end

            testCase.verifyLessThan(errors(2:end), errors(1:end-1) / 10, ...
                'Spectral convergence should be observed for the oblate case as the order p is increased.');
        end

        %%% Dirichlet problems
        function testExteriorDirichletProblemOneSpheroid(testCase)
            p = 16;
            eta = 10;
            ns = 1; % Number of spheroids
            u0 = 1.8; % 1/eccentricity of a prolate spheroidal surface
            target_distances = 1e-6; % Distance from the surface to evaluate the potential
            plt = false;
            neumann = false;
            interior = false;
            
            [soln, fluxsoln, truesoln, truefluxSurf, ~, ~] = spheroidalDP_charge_problem(p, eta, ns, u0, target_distances, plt, neumann, interior);
            testCase.verifyLessThan(norm(fluxsoln - truefluxSurf)/norm(truefluxSurf), 1e-7, ...
                'spheroidalDP should match true flux.');
        end

        function testExteriorDirichletProblemMultipleSpheroids(testCase)
            p = 16;
            eta = 10;
            ns = 3; % Number of spheroids
            u0 = [1.8 1.8 1.8]; % 1/eccentricity of a prolate spheroidal surface
            target_distances = 1e-4*ones(1, ns); % Distance from the surface to evaluate the potential
            plt = false;
            neumann = false;
            interior = false;
            
            [soln, fluxsoln, truesoln, truefluxSurf, ~, ~] = spheroidalDP_charge_problem(p, eta, ns, u0, target_distances, plt, neumann, interior);
            testCase.verifyLessThan(norm(fluxsoln - truefluxSurf)/norm(truefluxSurf), 1e-7, ...
                'spheroidalDP should match true flux.');
        end

        function testInteriorDirichletProblem(testCase)
            %%% One spheroid
            p = 16;
            eta = 10;
            ns = 1; % Number of spheroids
            u0 = 1.8; % 1/eccentricity of a prolate spheroidal surface
            target_distances = 1e-6; % Distance from the surface to evaluate the potential
            plt = false;
            neumann = false;
            interior = true;
            
            [soln, truesoln, truefluxSurf, sigma_vec, condK] =  spheroidalDP_charge_problem(p, eta, ns, u0, target_distances, plt, neumann, interior);

            %%% Multiple spheroids
        end
    end
end

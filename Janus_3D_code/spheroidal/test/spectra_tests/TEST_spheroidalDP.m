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
            target_u0 = params_ref.u0 * 1.2;
            X_trg = prolate_spheroid_shape(testCase.p_max, target_u0, params_ref.a);
            % Radial normal vectors out of spheroid
            nu_trg = get_norm_vecs(testCase.p_max, target_u0, params_ref.oblate);

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
            target_u0 = params_ref.u0 * 1.2;
            X_trg = oblate_spheroid_shape(testCase.p_max, target_u0, params_ref.a);
            % Radial normal vectors out of spheroid
            nu_trg = get_norm_vecs(testCase.p_max, target_u0, params_ref.oblate);

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
            p = 1;

            params = SpheroidalParameters;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
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

            target_u0 = params.u0 * 3;
            X_trg = prolate_spheroid_shape(p, target_u0, params.a);
            % Radial normal vectors out of spheroid
            nu_trg = get_norm_vecs(p, target_u0, params.oblate);

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
            p = 1;
            eps = testCase.fd_eps;

            params = SpheroidalParameters;
            params.u0 = testCase.u0_oblate;
            params.a = testCase.a_oblate;
            params.oblate = true;
            [u_p, v_p] = gl_grid(p);
            params.sigma = testCase.density_func(u_p, v_p);
            params.get_shc();

            target_u0 = params.u0 * 3;
            X_trg = oblate_spheroid_shape(p, target_u0, params.a);
            % Radial normal vectors out of spheroid
            nu_trg = get_norm_vecs(p, target_u0, params.oblate);

            DP_spectral = spheroidalDP(params, X_trg, nu_trg);

            DL_plus = spheroidalDL(params, X_trg + eps * nu_trg);
            DL_minus = spheroidalDL(params, X_trg - eps * nu_trg);
            DP_fd = (DL_plus - DL_minus) / (2 * eps);

            rel_err = norm(DP_spectral - DP_fd, inf) / norm(DP_spectral, inf);
            testCase.verifyLessThan(rel_err, testCase.gradient_check_tol, ...
                'Gradient check failed for oblate case: should be below tolerance.');
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

            target_u0 = params.u0 * 3;
            X_trg = prolate_spheroid_shape(p, target_u0, params.a);
            nu_trg = get_norm_vecs(p, target_u0, params.oblate);

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

                Note that this requires higher order (i.e. 'p') for the 
                desired convergence.
            %}
            p = 16;
            
            params = SpheroidalParameters;
            params.isReal = true;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = true;
            [u_p, v_p] = gl_grid(p);
            params.sigma = testCase.density_func(u_p, v_p);
            params.get_shc();

            target_u0 = params.u0 * 3;
            X_trg = prolate_spheroid_shape(p, target_u0, params.a);
            nu_trg = get_norm_vecs(p, target_u0, params.oblate);

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

        %%% On-surface checks
        function testProlateOnSurface(testCase)
        end
    end
end

%{
    Test code for spheroidalSP.m.

    Tests spectral convergence on-surface and off-surface, and has
    a gradient check using finite differences.

    The convergence tests are quite slow on this one, so some tinkering
    with the orders is probably needed.
%}


classdef TEST_spheroidalSP < matlab.unittest.TestCase
    properties
        u0_prolate = 2/sqrt(3);
        u0_oblate = 2/sqrt(3);
        a_prolate;
        a_oblate;
        p_max = 32; % Maximum order for convergence tests
        conv_tol = 1e-5;
        tol = 1e-6;

        % Tolerance for consistency checks
        consistency_tol = 1e-12;

        % Tolerance for gradient checks
        gradient_check_tol = 1e-6;
        fd_eps = 1e-8;

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

            SP_ref = spheroidalSP(params_ref, X_trg, nu_trg);

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

                SP_p = spheroidalSP(params_p, X_trg, nu_trg);

                errors(i) = norm(SP_p - SP_ref, inf) / norm(SP_ref, inf);
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

            SP_ref = spheroidalSP(params_ref, X_trg, nu_trg);

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

                SP_p = spheroidalSP(params_p, X_trg, nu_trg);

                errors(i) = norm(SP_p - SP_ref, inf) / norm(SP_ref, inf);
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
            SP_ref = spheroidalSP(params_ref);

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

                SP_p_full = spheroidalSP(params_p);

                % We use spherical transforms for spheroidal transforms, so this should be OK.
                SP_p_interp = interpsh(SP_p_full, testCase.p_max);
                errors(i) = norm(SP_p_interp - SP_ref, inf) / norm(SP_ref, inf);
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
            SP_ref = spheroidalSP(params_ref);

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

                SP_p_full = spheroidalSP(params_p);

                SP_p_interp = interpsh(SP_p_full, testCase.p_max);
                errors(i) = norm(SP_p_interp - SP_ref, inf) / norm(SP_ref, inf);
            end

            testCase.verifyLessThan(errors(end), testCase.conv_tol, ...
                'Error for oblate on-surface should be below tolerance.');
        end

        %%% Consistency checks
        function testCompareProlateOnSurfaceImplementations(testCase)
            % Compares the on-surface implementation with the general 
            % implementation evaluated on the surface (X=[], nu=nu_outward).
            p = 8;

            params = SpheroidalParameters;
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
            SP_method1 = spheroidalSP(params, X_src, nu_src);

            % Method 2: handle []
            SP_method2 = spheroidalSP(params, [], nu_src);

            % Method 3: no inputs
            SP_method3 = spheroidalSP(params);

            testCase.verifyEqual(SP_method1, SP_method2, 'AbsTol', testCase.consistency_tol, ...
                'Off-surface code and [] handling does not match each other.');

            testCase.verifyEqual(SP_method2, SP_method3, 'AbsTol', testCase.consistency_tol, ...
                'On-surface code and [] handling does not match each other.');
        end

        function testCompareOblateOnSurfaceImplementations(testCase)
            % Compares the on-surface implementation with the general 
            % implementation evaluated on the surface (X=[], nu=nu_outward).
            p = 8;

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
            SP_method1 = spheroidalSP(params, X_src, nu_src);

            % Method 2: handle []
            SP_method2 = spheroidalSP(params, [], nu_src);

            % Method 3: no inputs
            SP_method3 = spheroidalSP(params);

            testCase.verifyEqual(SP_method1, SP_method2, 'AbsTol', testCase.consistency_tol, ...
                'Off-surface code and [] handling does not match each other.');

            testCase.verifyEqual(SP_method2, SP_method3, 'AbsTol', testCase.consistency_tol, ...
                'On-surface code and [] handling does not match each other.');
        end

        %%% Gradient checks
        function testGradientCheckProlateOffSurface(testCase)
            p = 4;
            eps = testCase.fd_eps;

            params = SpheroidalParameters;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = false;
            [u_p, v_p] = gl_grid(p);
            params.sigma = testCase.density_func(u_p, v_p);
            params.get_shc();

            target_u0 = params.u0 * 1.5;
            X_trg = prolate_spheroid_shape(testCase.p_max, target_u0, params.a);
            % Radial normal vectors out of spheroid
            nu_trg = get_norm_vecs(testCase.p_max, target_u0, params.oblate);

            SP_spectral = spheroidalSP(params, X_trg, nu_trg);

            X_plus = X_trg + eps * nu_trg;
            X_minus = X_trg - eps * nu_trg;

            SL_plus = spheroidalSL(params, X_plus);
            SL_minus = spheroidalSL(params, X_minus);
            SP_fd = (SL_plus - SL_minus) / (2 * eps);

            rel_err = norm(SP_spectral - SP_fd, inf) / norm(SP_spectral, inf);
            testCase.verifyLessThan(rel_err, testCase.gradient_check_tol, ...
                'Gradient check failed for prolate case.');
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

            target_u0 = params.u0 * 1.5;
            X_trg = prolate_spheroid_shape(testCase.p_max, target_u0, params.a);
            % Radial normal vectors out of spheroid
            nu_trg = get_norm_vecs(testCase.p_max, target_u0, params.oblate);

            SP_spectral = spheroidalSP(params, X_trg, nu_trg);

            SL_plus = spheroidalSL(params, X_trg + eps * nu_trg);
            SL_minus = spheroidalSL(params, X_trg - eps * nu_trg);
            SP_fd = (SL_plus - SL_minus) / (2 * eps);

            rel_err = norm(SP_spectral - SP_fd, inf) / norm(SP_spectral, inf);
            testCase.verifyLessThan(rel_err, testCase.gradient_check_tol, ...
                'Gradient check failed for oblate case.');
        end

        function testGradientCheckProlateOnSurface(testCase)
            p = 16;
            eps = testCase.fd_eps;

            params = SpheroidalParameters;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = false;
            [u_p, v_p] = gl_grid(p);
            params.sigma = testCase.density_func(u_p, v_p);
            params.get_shc();

            SP_spectral = spheroidalSP(params);

            X_src = prolate_spheroid_shape(p, params.u0, params.a);
            nu_src = get_norm_vecs(p, params.u0, params.oblate);

            X_plus = X_src + eps * nu_src;
            X_minus = X_src - eps * nu_src;

            SL_plus = spheroidalSL(params, X_plus);
            SL_minus = spheroidalSL(params, X_minus);
            SP_fd = (SL_plus - SL_minus) / (2 * eps);

            rel_err = norm(SP_spectral - SP_fd, inf) / norm(SP_spectral, inf);
            testCase.verifyLessThan(rel_err, testCase.gradient_check_tol, ...
                'Prolate on-surface gradient check failed to meet tolerance.');
        end

        function testGradientCheckOblateOnSurface(testCase)
            p = 16;
            eps = testCase.fd_eps;

            params = SpheroidalParameters;
            params.u0 = testCase.u0_oblate;
            params.a = testCase.a_oblate;
            params.oblate = true;
            [u_p, v_p] = gl_grid(p);
            params.sigma = testCase.density_func(u_p, v_p);
            params.get_shc();

            SP_spectral = spheroidalSP(params);

            X_src = oblate_spheroid_shape(p, params.u0, params.a);
            nu_src = get_norm_vecs(p, params.u0, params.oblate);

            X_plus = X_src + eps * nu_src;
            X_minus = X_src - eps * nu_src;

            SL_plus = spheroidalSL(params, X_plus);
            SL_minus = spheroidalSL(params, X_minus);
            SP_fd = (SL_plus - SL_minus) / (2 * eps);

            rel_err = norm(SP_spectral - SP_fd, inf) / norm(SP_spectral, inf);
            testCase.verifyLessThan(rel_err, testCase.gradient_check_tol, ...
                'Prolate on-surface gradient check failed to meet tolerance.');
        end
    end
end

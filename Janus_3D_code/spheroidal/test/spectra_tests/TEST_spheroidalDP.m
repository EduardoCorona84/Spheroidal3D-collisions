%{
    Test code for spheroidalDP.m.

    Tests spectral convergence on-surface and off-surface, and has
    a gradient check using finite differences.
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

        % Tolerance for gradient checks
        gradient_check_tol = 1e-6;
        fd_eps = 1e-4;

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
            X_trg = oblate_spheroid_shape(testCase.p_max, target_u0, params.a);
            % Radial normal vectors out of spheroid
            nu_trg = get_norm_vecs(testCase.p_max, target_u0, params.oblate);

            DP_spectral = spheroidalDP(params, X_trg, nu_trg);

            DL_plus = spheroidalDL(params, X_trg + eps * nu_trg);
            DL_minus = spheroidalDL(params, X_trg - eps * nu_trg);
            DP_fd = (DL_plus - DL_minus) / (2 * eps);

            rel_err = norm(DP_spectral - DP_fd, inf) / norm(DP_spectral, inf);
            testCase.verifyLessThan(rel_err, testCase.gradient_check_tol, ...
                'Gradient check failed for oblate case: should be below tolerance.');
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

            DP_spectral = spheroidalDP(params);

            X_src = prolate_spheroid_shape(p, params.u0, params.a);
            nu_src = get_norm_vecs(p, params.u0, params.oblate);

            X_plus = X_src + eps * nu_src;
            X_minus = X_src - eps * nu_src;

            DL_plus = spheroidalDL(params, X_plus);
            DL_minus = spheroidalDL(params, X_minus);
            DP_fd = (DL_plus - DL_minus) / (2 * eps);

            rel_err = norm(DP_spectral - DP_fd, inf) / norm(DP_spectral, inf);
            testCase.verifyLessThan(rel_err, testCase.gradient_check_tol, ...
                'Prolate on-surface gradient check failed to meet tolerance.');
        end

        function testGradientCheckOblateOnSurface(testCase)
            p = 16;
            eps = testCase.fd_eps;

            params = SpheroidalParameters;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = true;
            [u_p, v_p] = gl_grid(p);
            params.sigma = testCase.density_func(u_p, v_p);
            params.get_shc();

            DP_spectral = spheroidalDP(params);

            X_src = oblate_spheroid_shape(p, params.u0, params.a);
            nu_src = get_norm_vecs(p, params.u0, params.oblate);

            X_plus = X_src + eps * nu_src;
            X_minus = X_src - eps * nu_src;

            DL_plus = spheroidalDL(params, X_plus);
            DL_minus = spheroidalDL(params, X_minus);
            DP_fd = (DL_plus - DL_minus) / (2 * eps);

            rel_err = norm(DP_spectral - DP_fd, inf) / norm(DP_spectral, inf);
            testCase.verifyLessThan(rel_err, testCase.gradient_check_tol, ...
                'Prolate on-surface gradient check failed to meet tolerance.');
        end
    end
end

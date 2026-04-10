classdef TEST_spheroidalModifiedDP < matlab.unittest.TestCase
    methods (TestClassSetup)
        function setupPaths(testCase)
            test_dir = fileparts(mfilename('fullpath'));
            repo_root = fullfile(test_dir, '..', '..', '..');

            addpath(genpath(fullfile(repo_root, 'spheroidal')));
            addpath(genpath(fullfile(repo_root, 'support')));
            addpath(fullfile(repo_root, 'support'));
        end
    end

    methods (Test)
        function testExteriorAnalyticalEigenvalue(testCase)
            lambda = 1;
            p = 12;
            u0 = 2 / sqrt(3);
            a = 1 / u0;
            n = 4;
            m = 2;

            gamma = 1j * lambda * a;
            [theta, phi] = gl_grid(p);
            v = cos(theta);
            params = LOCAL_build_params(p, u0, a, false, ASWFnm(n, m, v, phi, gamma, p, 0));

            X_src = prolate_spheroid_shape(p, u0, a);
            nu_src = get_norm_vecs(p, u0, false);
            X_trg = X_src + 1 * nu_src;
            modDP = spheroidalModifiedDP(params, lambda, X_trg);

            S = cart2spheroidal(X_trg, a, false);
            u_trg = S(:, 1);
            v_trg = S(:, 2);
            phi_trg = S(:, 3);
            Snm_trg = ASWFnm(n, m, v_trg, phi_trg, gamma, p, 0);

            [~, dR1_u0, ~, ~] = ...
                modified_laplace_test_radial_mode(p, n, m, u0, gamma, false);
            [~, ~, ~, dR3_u] = ...
                modified_laplace_test_radial_mode(p, n, m, u_trg, gamma, false);
            anm = (1i * gamma * (u0^2 - 1) / a) .* sqrt((u_trg.^2 - 1) ./ (u_trg.^2 - v_trg.^2));
            expected = anm .* dR1_u0 .* dR3_u .* Snm_trg;

            rel_err = norm(modDP - expected) / norm(expected);
            testCase.verifyLessThan(rel_err, 1e-11, ...
                'Exterior analytical eigenvalue check failed for spheroidalModifiedDP.');
        end

        function testInteriorAnalyticalEigenvalue(testCase)
            lambda = 1;
            p = 16;
            u0 = 2 / sqrt(3);
            a = 1 / u0;
            n = 6;
            m = 2;

            gamma = 1j * lambda * a;
            [theta, phi] = gl_grid(p);
            v = cos(theta);
            params = LOCAL_build_params(p, u0, a, false, ASWFnm(n, m, v, phi, gamma, p, 0));

            X_src = prolate_spheroid_shape(p, u0, a);
            nu_src = get_norm_vecs(p, u0, false);
            X_trg = X_src - 0.1 * nu_src;
            modDP = spheroidalModifiedDP(params, lambda, X_trg);
            modDP = reshape(modDP, [], 1);

            S = cart2spheroidal(X_trg, a, false);
            u_trg = S(:, 1);
            v_trg = S(:, 2);
            phi_trg = S(:, 3);
            Snm_trg = ASWFnm(n, m, v_trg, phi_trg, gamma, p, 0);

            [~, ~, ~, dR3_u0] = ...
                modified_laplace_test_radial_mode(p, n, m, u0, gamma, false);
            [~, dR1_u, ~, ~] = ...
                modified_laplace_test_radial_mode(p, n, m, u_trg, gamma, false);
            anm = (1i * gamma * (u0^2 - 1) / a) .* sqrt((u_trg.^2 - 1) ./ (u_trg.^2 - v_trg.^2));
            expected = anm .* dR3_u0 .* dR1_u .* Snm_trg;

            rel_err = norm(modDP - expected) / norm(expected);
            testCase.verifyLessThan(rel_err, 1e-12, ...
                'Interior analytical eigenvalue check failed for spheroidalModifiedDP.');
        end

        function testExteriorEvaluationMatchesKernelEval(testCase)
            %{
            dDP is a hypersingular operator, so it's more sensitive to high
            frequencies (hence the 1.5 scale factor instead of 1) 
            %}
            lambda = 1;
            p = 16;
            u0 = 2 / sqrt(3);
            a = 1 / u0;
            gamma = 1j * lambda * a;

            [theta, phi] = gl_grid(p);
            v = cos(theta);
            sigma = ASWFnm(5, 2, v, phi, gamma, p, 0);

            params = LOCAL_build_params(p, u0, a, false, sigma);
            X_src = prolate_spheroid_shape(p, u0, a);
            nu_src = get_norm_vecs(p, u0, false);
            X_trg = X_src + 1.5 * nu_src;
            S_trg = cart2spheroidal(X_trg, a, false);
            nu_trg = spheroidalNu2cart(repmat([1 0 0], size(X_trg, 1), 1), S_trg, a, false);

            modDP = spheroidalModifiedDP(params, lambda, X_trg);

            [~, gwt_gl] = g_grid(p + 1);
            wt_gl = pi / p * repmat(gwt_gl', 2 * p, 1) ./ sin(gl_grid(p));
            wt_gl = wt_gl(:);
            W_src = SurfaceSph(X_src).geoProp.W .* wt_gl;

            KEparams = Kernel_Eval_parameters('dDL_LMOD_3D', 0, 1, 1, 1, 1e-12, 2, 400, 1);
            KEparams.dim = 3;
            KEparams.X = X_src;
            KEparams.W2 = W_src.';
            KEparams.nor = get_norm_vecs(p, u0, false);
            KEparams.targnor = nu_trg;
            KEparams.lambda = lambda;

            modDP_KE = Kernel_Eval(X_trg, X_src, KEparams) * sigma;

            rel_err = norm(modDP - modDP_KE) / norm(modDP_KE);
            testCase.verifyLessThan(rel_err, 1e-10, ...
                'Exterior spheroidalModifiedDP does not match Kernel_Eval.');
        end
    end
end

function params = LOCAL_build_params(p, u0, a, oblate, sigma)
    params = SpheroidalParameters;
    params.p = p;
    params.isReal = false;
    params.u0 = u0;
    params.a = a;
    params.oblate = oblate;
    params.sigma = sigma;
    params.get_shc();
end

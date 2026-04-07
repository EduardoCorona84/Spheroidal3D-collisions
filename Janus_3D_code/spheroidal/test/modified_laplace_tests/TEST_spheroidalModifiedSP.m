classdef TEST_spheroidalModifiedSP < matlab.unittest.TestCase
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
            p = 16;
            u0 = 2 / sqrt(3);
            a = 1 / u0;
            n = 4;
            m = 2;

            c = 1j * lambda * a;
            [theta, phi] = gl_grid(p);
            v = cos(theta);
            sigma = LOCAL_weighted_Snm(p, u0, c, n, m, v, phi, false);
            params = LOCAL_build_params(p, u0, a, false, sigma);

            X_src = prolate_spheroid_shape(p, u0, a);
            nu_src = get_norm_vecs(p, u0, false);
            X_trg = X_src + 1 * nu_src;
            modSP = spheroidalModifiedSP(params, lambda, X_trg);
            modSP = reshape(modSP, [], 1);

            S = cart2spheroidal(X_trg, a, false);
            u_trg = S(:, 1);
            v_trg = real(S(:, 2));
            phi_trg = S(:, 3);
            Snm_trg = ASWFnm(n, m, v_trg, phi_trg, c, p, 0);

            [R1_u0, ~, ~, ~] = ...
                modified_laplace_test_radial_mode(p, n, m, u0, c, false);
            [~, ~, ~, dR3_u] = ...
                modified_laplace_test_radial_mode(p, n, m, u_trg, c, false);
            anm = (1i * c * sqrt(u0^2 - 1)) .* sqrt((u_trg.^2 - 1) ./ (u_trg.^2 - v_trg.^2));
            expected = anm .* R1_u0 .* dR3_u .* Snm_trg;

            rel_err = norm(modSP - expected) / norm(expected);
            testCase.verifyLessThan(rel_err, 1e-11, ...
                'Exterior analytical eigenvalue check failed for spheroidalModifiedSP.');
        end

        function testInteriorAnalyticalEigenvalue(testCase)
            lambda = 1;
            p = 16;
            u0 = 2 / sqrt(3);
            a = 1 / u0;
            n = 4;
            m = 2;

            c = 1j * lambda * a;
            [theta, phi] = gl_grid(p);
            v = cos(theta);
            sigma = LOCAL_weighted_Snm(p, u0, c, n, m, v, phi, false);
            params = LOCAL_build_params(p, u0, a, false, sigma);

            X_src = prolate_spheroid_shape(p, u0, a);
            nu_src = get_norm_vecs(p, u0, false);
            X_trg = X_src - 0.1 * nu_src;
            modSP = spheroidalModifiedSP(params, lambda, X_trg);
            modSP = reshape(modSP, [], 1);

            S = cart2spheroidal(X_trg, a, false);
            u_trg = S(:, 1);
            v_trg = real(S(:, 2));
            phi_trg = S(:, 3);
            Snm_trg = ASWFnm(n, m, v_trg, phi_trg, c, p, 0);

            [~, ~, R3_u0, ~] = ...
                modified_laplace_test_radial_mode(p, n, m, u0, c, false);
            [~, dR1_u, ~, ~] = ...
                modified_laplace_test_radial_mode(p, n, m, u_trg, c, false);
            anm = (1i * c * sqrt(u0^2 - 1)) .* sqrt((u_trg.^2 - 1) ./ (u_trg.^2 - v_trg.^2));
            expected = anm .* R3_u0 .* dR1_u .* Snm_trg;

            rel_err = norm(modSP - expected) / norm(expected);
            testCase.verifyLessThan(rel_err, 1e-11, ...
                'Interior analytical eigenvalue check failed for spheroidalModifiedSP.');
        end

        function testExteriorEvaluationMatchesKernelEval(testCase)
            lambda = 1;
            p = 16;
            u0 = 2 / sqrt(3);
            a = 1 / u0;
            c = 1j * lambda * a;

            [theta, phi] = gl_grid(p);
            v = cos(theta);
            sigma = ASWFnm(0, 0, v, phi, c, p, 0);
            sigma = sigma + 0.4 * ASWFnm(2, 1, v, phi, c, p, 0);
            sigma = sigma - 0.2 * ASWFnm(3, 2, v, phi, c, p, 0);
            sigma = sigma(:) / max(1, norm(sigma));

            params = LOCAL_build_params(p, u0, a, false, sigma);
            X_src = prolate_spheroid_shape(p, u0, a);
            nu_src = get_norm_vecs(p, u0, false);
            X_trg = X_src + 1.5 * nu_src;
            S_trg = cart2spheroidal(X_trg, a, false);
            nu_trg = spheroidalNu2cart(repmat([1 0 0], size(X_trg, 1), 1), S_trg, a, false);

            modSP = spheroidalModifiedSP(params, lambda, X_trg);
            modSP = reshape(modSP, [], 1);

            S_src = SurfaceSph(X_src);
            [~, gwt_gl] = g_grid(p + 1);
            wt_gl = pi / p * repmat(gwt_gl', 2 * p, 1) ./ sin(gl_grid(p));
            wt_gl = wt_gl(:);
            W_src = S_src.geoProp.W .* wt_gl;

            KEparams = Kernel_Eval_parameters('dSL_LMOD_3D', 0, 1, 1, 1, 1e-12, 2, 400, 1);
            KEparams.dim = 3;
            KEparams.X = X_src;
            KEparams.W2 = W_src.';
            KEparams.nor = nu_trg;
            KEparams.lambda = lambda;

            modSP_KE = Kernel_Eval(X_trg, X_src, KEparams) * sigma;

            rel_err = norm(modSP - modSP_KE) / norm(modSP_KE);
            testCase.verifyLessThan(rel_err, 1e-6, ...
                'Exterior spheroidalModifiedSP does not match Kernel_Eval.');
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

function sigma = LOCAL_weighted_Snm(p, u0, c, n, m, v, phi, oblate)
    sp = (p + 1)^2;
    ASWFmtx = LOCAL_build_swf_matrix(p, v, phi, c);
    idx = LOCAL_geti(n, m);
    
    coeffs = zeros(sp, 1);
    coeffs(idx) = 1;
    
    G = ASWF_Gmatrix(p, u0, c, 0, oblate, 50, 0);
    Gswfc = G * coeffs;
    sigma = ASWFmtx * Gswfc;
end

function S = LOCAL_build_swf_matrix(p, v, phi, c)
    sp = (p + 1)^2;
    v = v(:);
    phi = phi(:);
    S = zeros(numel(v), sp);
    for n = 0:p
        Sn = ASWFnm(n, [], v, phi, c, p, 0);
        S(:, n^2 + 1:(n + 1)^2) = Sn;
    end
end

function idx = LOCAL_geti(n, m)
    idx = m + n^2 + n + 1;
end

classdef TEST_spheroidalModifiedSLP < matlab.unittest.TestCase
%{
Tests for spheroidalModifiedSLP with weighted ASWF_Gmatrix basis transform.
%}
    methods (Test)
        function testOnSurfaceAnalyticalEigenvalue(testCase)
            lambda = 1;
            p = 8;
            n = 4;
            m = 2;
            u0 = 2 / sqrt(3);
            a = 1 / u0;
            gamma = 1j * lambda * a;

            [theta, phi] = gl_grid(p);
            v = cos(theta);
            Ssurf = LOCAL_build_swf_matrix(p, v, phi, gamma);

            sp = (p + 1)^2;
            idx = LOCAL_geti(n, m);
            weighted_coeff = zeros(sp, 1);
            weighted_coeff(idx) = 1;

            G = ASWF_Gmatrix(p, u0, gamma, 0, false, 50, 0);
            swfc = G * weighted_coeff;
            sigma = Ssurf * swfc;

            params = SpheroidalParameters;
            params.p = p;
            params.isReal = false;
            params.u0 = u0;
            params.a = a;
            params.oblate = false;
            params.sigma = sigma;
            params.get_shc();

            modSL = spheroidalModifiedSLP(params, lambda, []);

            [~, lambda_surf, ~] = LOCAL_modSLPspectrum(p, u0, a, gamma);
            expected = lambda_surf(idx) * Ssurf(:, idx);

            rel_err = norm(modSL - expected) / norm(expected);
            testCase.verifyLessThan(rel_err, 1e-12, ...
                'On-surface analytical eigenvalue check failed for modified SLP.');
        end

        function testFarEvaluationMatchesKernelEval(testCase)
            lambda = 2;
            p = 16;

            params = SpheroidalParameters;
            params.p = p;
            params.isReal = false;
            params.u0 = 2/sqrt(3);
            params.a = 1 / params.u0;
            params.oblate = false;

            gamma = 1j * lambda * params.a;
            [theta, phi] = gl_grid(p);
            v = cos(theta);
            % Note that 2p is needed here to get desired order.
            sigma = ASWFnm(5, 2, v, phi, gamma, 2*p, 0);;
            params.sigma = sigma;
            params.get_shc();

            X_src = prolate_spheroid_shape(p, params.u0, params.a);
            nu_src = params.get_Norm(p, 1);
            X_trg = X_src + 1*nu_src;

            modSL_spectral = spheroidalModifiedSLP(params, lambda, X_trg);

            Sns = SurfaceSph(X_src);
            [~, gwt_gl] = g_grid(p + 1);
            wt_gl = pi / p * repmat(gwt_gl', 2 * p, 1) ./ sin(gl_grid(p));
            wt_gl = wt_gl(:);
            W_src_orig = Sns.geoProp.W .* wt_gl;
            
            KEparams = Kernel_Eval_parameters('SL_LMOD_3D', 0, 1, 1, 1, 1e-12, 2, 400, 1);
            KEparams.dim = 3;
            KEparams.X = X_src;
            KEparams.W2 = W_src_orig.';
            KEparams.lambda = lambda;
            modSL_KE = Kernel_Eval(X_trg, X_src, KEparams) * sigma;

            rel_err = norm(modSL_spectral - modSL_KE) / norm(modSL_KE);
            testCase.verifyLessThan(rel_err, 1e-8, ...
                'Mismatch between spheroidalModifiedSLP and Kernel_Eval.');
        end

        function testOnSurfaceEvaluationMatchesRBS(testCase)
            lambda = 1;
            p = 16;

            params = SpheroidalParameters;
            params.p = p;
            params.isReal = false;
            params.u0 = 2;
            params.a = 1 / params.u0;
            params.oblate = false;

            gamma = 1j * lambda * params.a;
            [theta, phi] = gl_grid(p);
            v = cos(theta);
            sigma = ASWFnm(5, 3, v, phi, gamma, 2*p, 0);
            params.sigma = sigma;
            params.get_shc();

            modSL_spectral = spheroidalModifiedSLP(params, lambda, []);

            X_src = prolate_spheroid_shape(p, params.u0, params.a);
            [modSL_mat, ~, ~] = kernelModifiedLap(SurfaceSph(X_src), 'SMat', lambda);
            modSL_direct = modSL_mat * sigma;

            rel_err = norm(modSL_spectral - modSL_direct) / norm(modSL_direct);
            testCase.verifyLessThan(rel_err, 1e-4, ...
                'Mismatch between spheroidalModifiedSLP and kernelModifiedLap.');
        end
    end
end

function idx = LOCAL_geti(n, m)
    idx = m + n.^2 + n + 1;
end

function S = LOCAL_build_swf_matrix(p, v, phi, gamma)
    sp = (p + 1)^2;
    v = v(:);
    phi = phi(:);
    S = zeros(numel(v), sp);
    for n = 0:p
        Sn = ASWFnm(n, [], v, phi, gamma, p, 0);
        S(:, n^2 + 1:(n + 1)^2) = Sn;
end
end

function [lambda_int, lambda_surf, lambda_ext] = LOCAL_modSLPspectrum(p, u0, a, c)
    sp = (p + 1)^2;
    cnm = 1j * a * c * sqrt(u0^2 - 1);

    Rnm1_vec = zeros(1, sp);
    Rnm3_vec = zeros(1, sp);

    for j = 0:p
        idx = j^2 + 1:(j + 1)^2;
        R1_j = Rnm1(j, [], u0, c);
        R3_j = Rnm3(j, [], u0, c);
        Rnm1_vec(idx) = R1_j;
        Rnm3_vec(idx) = R3_j;
    end

    lambda_int = cnm .* Rnm3_vec;
    lambda_surf = cnm .* (Rnm1_vec .* Rnm3_vec);
    lambda_ext = cnm .* Rnm1_vec;
end

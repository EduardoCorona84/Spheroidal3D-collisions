classdef TEST_ASWFnm < matlab.unittest.TestCase
%{
Tests for angular spheroidal wave function ASWFnm.
%}
    methods (Test)
        function testMatchesManualLegendreSum(testCase)
            n = 5;
            m = 2;
            c = 1.2;
            p = 24;
            v = [-0.8; -0.3; 0.1; 0.6];
            phi = 0.4;

            S = ASWFnm(n, m, v, phi, c, p, 1);
            S_manual = manual_legendre_sum(n, m, v, phi, c, p);

            rel_err = norm(S - S_manual) / norm(S_manual);
            testCase.verifyLessThan(rel_err, 1e-10, ...
                'ASWFnm does not match Legendre reconstruction.');
        end

        function testSmallCMatchesLegendre(testCase)
            n = 4;
            m = 2;
            c = 1e-8;
            p = 16;
            v = linspace(-0.7, 0.7, 20).';
            phi = 0.9;

            S = ASWFnm(n, m, v, phi, c, p, 1);
            u = real(acos(v));
            expected = sqrt(2*pi) * Ynm(n, m, u, phi);
            rel_err = norm(S - expected) / norm(expected);
            testCase.verifyLessThan(rel_err, 1e-8, ...
                'Small-c limit does not match Ynm up to scale.');
        end

        function testEmptyMMatchesPerM(testCase)
            n = 3;
            c = 0.7;
            p = 16;
            v = linspace(-0.6, 0.6, 12).';
            phi = 0.2;

            S_all = ASWFnm(n, [], v, phi, c, p, 0);
            for m = -n:n
                S_m = ASWFnm(n, m, v, phi, c, p, 0);
                idx = m + n + 1;
                rel_err = norm(S_all(:, idx) - S_m) / max(1, norm(S_m));
                testCase.verifyLessThan(rel_err, 1e-12, ...
                    sprintf('m=[] column mismatch for m=%d.', m));
            end
        end

        function testNormalizationConstant(testCase)
            n = 4;
            m = 2;
            gamma = -1j*0.3;
            p = 24;

            [eta, w] = g_grid(p);
            eta = eta(:);
            w = w(:);

            S_norm0 = ASWFnm(n, m, eta, 0, gamma, p, 0);
            S_norm1 = ASWFnm(n, m, eta, 0, gamma, p, 1);

            S_norm0_L2 = sum(w .* (abs(S_norm0) .^ 2));
            S_norm1_L2 = sum(w .* (abs(S_norm1) .^ 2));

            expected_normalization = 2 / (2 * n + 1) * factorial(n + abs(m)) / factorial(n - abs(m));

            testCase.verifyLessThan(abs(S_norm0_L2 - expected_normalization), 1e-12, ...
                'Normalization for iopnorm = 0 is incorrect.');
            testCase.verifyLessThan(abs(S_norm1_L2 - 1), 1e-12, ...
                'Normalization for iopnorm = 1 is incorrect.');
        end

        function testNegativeImagGammaBranch(testCase)
            n = 7;
            m = 3;
            gamma = -0.6i;
            p = 48;
            v = linspace(-0.8, 0.8, 31).';
            phi = linspace(0, 2*pi, numel(v)).';

            S = ASWFnm(n, m, v, phi, gamma, p, 1);
            S_manual = manual_legendre_sum(n, m, v, phi, gamma, p);

            rel_err = norm(S - S_manual) / max(1, norm(S_manual));
            testCase.verifyLessThan(rel_err, 1e-10, ...
                'ASWFnm imag(gamma)<0 branch convention is inconsistent.');
        end

        function testODEResidualRealGamma(testCase)
            n = 4;
            m = 2;
            gamma = 1.1;
            p = 32;
            N = 100;
            phi = 0;

            rel_err = LOCAL_eval_aswf_ode(n, m, gamma, p, N, phi);
            testCase.verifyLessThan(rel_err, 1e-10, ...
                'ASWFnm does not satisfy the angular ODE for real gamma.');
        end

        function testODEResidualComplexGamma(testCase)
            n = 6;
            m = 2;
            gamma = 3-0.6i;
            p = 32;
            N = 100;
            phi = 0;

            rel_err = LOCAL_eval_aswf_ode(n, m, gamma, p, N, phi);
            testCase.verifyLessThan(rel_err, 1e-10, ...
                'ASWFnm does not satisfy the angular ODE for complex gamma.');
        end
    end
end

function S = manual_legendre_sum(n, m, v, phi, c, p)
    mm = abs(m);
    parity = mod(n - mm, 2);
    nu = mm + parity;
    j = (n - nu) / 2 + 1;

    [D, ~] = leg_to_pswf_mtx(mm, p, c, parity);
    coeffs = D(:, j);
    degs = (nu:2:p).';

    S = zeros(numel(v), 1);
    for k = 1:numel(degs)
        Pk = legendre2(degs(k), v);
        S = S + coeffs(k) * Pk(mm + 1, :).';
    end

    S = S .* exp(1i * m * phi);
end

function rel_err = LOCAL_eval_aswf_ode(n, m, gamma, p, N, phi)
    [eta, D] = LOCAL_cheb(N);
    D2 = D * D;

    S = ASWFnm(n, m, eta, phi, gamma, p, 1);

    mm = abs(m);
    parity = mod(n - mm, 2);
    nu = mm + parity;
    j = (n - nu) / 2 + 1;
    [~, lambda] = leg_to_pswf_mtx(mm, p, gamma, parity);
    lambda_nm = lambda(j);

    idx = 2:N;
    eta_i = eta(idx);
    res = (1 - eta_i.^2) .* (D2(idx, :) * S) - 2 * eta_i .* (D(idx, :) * S) + ...
          (lambda_nm + gamma^2 .* (1 - eta_i.^2) - (m^2 ./ (1 - eta_i.^2))) .* S(idx);

    rel_err = norm(res) / max(1, norm(S(idx)));
end

function [x, D] = LOCAL_cheb(N)
% Chebyshev derivative matrix from Trefthen
    if N == 0
        x = 1;
        D = 0;
        return;
    end
    k = (0:N).';
    x = cos(pi * k / N);
    c = [2; ones(N-1,1); 2] .* (-1).^k;
    X = repmat(x, 1, N+1);
    dX = X - X.';
    D = (c * (1./c).') ./ (dX + eye(N+1));
    D = D - diag(sum(D, 2));
end

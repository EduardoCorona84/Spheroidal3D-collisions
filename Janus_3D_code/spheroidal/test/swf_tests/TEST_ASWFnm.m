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

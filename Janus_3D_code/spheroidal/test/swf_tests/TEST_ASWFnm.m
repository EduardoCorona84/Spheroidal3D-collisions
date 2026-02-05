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

            S = ASWFnm(n, m, v, phi, c, p);
            S_manual = manual_legendre_sum(n, m, v, phi, c, p);

            relErr = norm(S - S_manual) / norm(S_manual);
            testCase.verifyLessThan(relErr, 1e-10, ...
                'ASWFnm does not match Legendre reconstruction.');
        end

        function testSmallCMatchesLegendre(testCase)
            n = 4;
            m = 2;
            c = 1e-8;
            p = 16;
            v = linspace(-0.7, 0.7, 20).';
            phi = 0.9;

            S = ASWFnm(n, m, v, phi, c, p);
            u = real(acos(v));
            expected = Ynm(n, m, u, phi);

            mm = abs(m);
            normalization = 2 / (2 * n + 1) * factorial(n + mm + 1) / factorial(n - mm + 1);
            expected = expected / sqrt(normalization);
            relErr = norm(S - expected) / norm(expected);
            testCase.verifyLessThan(relErr, 1e-8, ...
                'Small-c limit does not match Ynm up to scale.');
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

    normalization = @(l) sqrt((2*l + 1) / (4*pi) * factorial(l-mm) / factorial(l+mm));

    S = zeros(numel(v), 1);
    for k = 1:numel(degs)
        Pk = legendre2(degs(k), v);
        S = S + coeffs(k) * (normalization(degs(k)) * Pk(mm + 1, :).');
    end

    S = S .* exp(1i * m * phi);
end

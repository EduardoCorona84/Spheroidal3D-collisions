%{
    Code to test the calculation of associated Legendre functions--namely 
    second derivatives. For old code that tests the zeroth and first
    orders, check legendre_otc_test.m.
%}

classdef TEST_legendre_otc < matlab.unittest.TestCase
    properties
        p = 16;
        %x_vals = linspace(1.1,3,100);
        x_vals = [1j*0.5]
        tol = 1e-6;
        L = {};
    end

    methods (TestClassSetup)
        function setup(testCase)
            clc();
            testCase.L = legendre_otc(testCase.p, testCase.x_vals, 1, 2, 2);
        end
    end

    methods (Test)
        %%% Gradient check
        function testSecondDerivativesByGradientCheck(testCase)
            geti = @(n,m) m+n^2+n+1;
            x_val = 1.4;
            n_eps = 6; % Number of unique epsilons to check (should be even)

            ddL = legendre_otc(testCase.p, x_val, 1, 2, 2);
            dPP = ddL{5}; dQQ = ddL{6};

            P_errors = zeros(1, n_eps); Q_errors = P_errors;
            for i = 1:n_eps
                epsilon = 10^(-i);
                L = legendre_otc(testCase.p, [x_val - epsilon x_val + epsilon], 1, 1, 1);
                dP = L{3}; dQ = L{4};

                P_errors(i) = norm((dP(:,2) - dP(:,1))/(2*epsilon) - dPP) ./ norm(dPP);
                Q_errors(i) = norm((dQ(:,2) - dQ(:,1))/(2*epsilon) - dQQ) ./ norm(dQQ);
            end

            testCase.verifyLessThan(P_errors(n_eps), 1e-10, "Calculation of P'' does not agree with central finite difference.");
            testCase.verifyLessThan(Q_errors(n_eps), 1e-10, "Calculation of P'' does not agree with central finite difference.");
        end

        %%% Formula check
        function testPSecondDerivatveFormula(testCase)
            %{
                Verifies that the second derivative of the Associated
                Legendre functions satisfy the following recursion
                relation:
                    (1 - x^2)P_n^{m \prime \prime}(x) - 2xP_n^{m \prime}(x) 
                        = (m-n-1) P_{n+1}^{m \prime}(x)+(n+1) P_n^m(x) 
                                + (n+1) P_{n}^{m \prime}(x)

                Reference:
                    https://dlmf.nist.gov/14.10
            %}
            debug_flag = false;
            x_vals = randn(1,1,like=1i) + 1;
            L = legendre_otc(testCase.p, x_vals, 1, 2, 2);
            P = L{1}; Q = L{2};
            dP = L{3}; dQ = L{4};
            ddP = L{5}; ddQ = L{6};

            geti = @(n,m) m+n^2+n+1;

            errors=zeros(size(P));
            for k = 1:(testCase.p)^2-1
                n = floor(sqrt(k-1)); m=k-n.^2-n-1;
                LHS_term = ddP(geti(n,m),:);
                RHS_term = 2 .* x_vals .* dP(geti(n,m),:) + (m - n - 1) .* dP(geti(n+1,m),:) + (n+1).*P(geti(n,m),:) + (n+1).*x_vals.*dP(geti(n,m),:);
                expression = (1 - x_vals.^2) .* LHS_term - RHS_term;
                errors(k) = min(norm(expression) / norm(RHS_term), norm(expression)); % Minimum of relative and absolute errors
                if (debug_flag && log10(expression) > -5)
                    disp(sprintf("potential issue on (n, m): (%d, %d); absolute: %f; however, relative error is %f", n, m, expression, errors(k)));
                end
            end

            testCase.verifyLessThan(errors(2:end,:), 1e-12, "P'' does not satisfy the known recurrence relation.");
        end

        function testQSecondDerivatveFormula(testCase)
            debug_flag = false;
            x_vals = randn(1,1,like=1i) + 1;
            L = legendre_otc(testCase.p, x_vals, 1, 2, 2);
            P = L{1}; Q = L{2};
            dP = L{3}; dQ = L{4};
            ddP = L{5}; ddQ = L{6};

            geti = @(n,m) m+n^2+n+1;

            errors=zeros(size(P));
            for k = 1:(testCase.p)^2-1
                n = floor(sqrt(k-1)); m=k-n.^2-n-1;
                LHS_term = ddQ(geti(n,m),:);
                RHS_term = 2 .* x_vals .* dQ(geti(n,m),:) + (m - n - 1) .* dQ(geti(n+1,m),:) + (n+1).*Q(geti(n,m),:) + (n+1).*x_vals.*dQ(geti(n,m),:);
                expression = (1 - x_vals.^2) .* LHS_term - RHS_term;
                errors(k) = min(norm(expression) / norm(RHS_term), norm(expression)); % Minimum of relative and absolute errors
                if (debug_flag && log10(expression) > -5)
                    disp(sprintf("potential issue on (n, m): (%d, %d); absolute: %f; however, relative error is %f", n, m, expression, errors(k)));
                end
            end

            testCase.verifyLessThan(errors(2:end,:), 1e-12, "Q'' does not satisfy the known recurrence relation.");
        end
    end
end

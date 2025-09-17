%{
    Scratchwork to develop tests for S''.
%}

classdef TEST_spheroidalMatVecKernel < matlab.unittest.TestCase
    properties
        u0_prolate = 2/sqrt(3);
        u0_oblate = 8/sqrt(3);
        a_prolate;
        a_oblate;
        tol = 1e-6;
    end

    methods (TestClassSetup)
        function setup(testCase)
            clc();
            testCase.a_prolate = 1/testCase.u0_prolate;
            testCase.a_oblate = 1/sqrt(1 + testCase.u0_oblate^2);
        end
    end

    methods (Test)
        function testEffectOfRotation(testCase)
            params = SpheroidalParameters;
            params.isReal = true;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = false;
            
            p = 16;
            np = 2*p*(p+1);

            params.sigma = eye(np);
            LPunrotated = spheroidalDL(params);

            params.thetas = pi/6;
            params.phis = pi/10;
            LProtated = spheroidalDL(params);

            testCase.verifyLessThan(norm(LPunrotated - LProtated), 1e-12, ...
                "In the case of one spheroid, rotation should have no effect.");
        end
    end
end

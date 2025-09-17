%{
    Scratchwork to develop tests for S''.
%}

classdef TEST_spheroidalMatVec < matlab.unittest.TestCase
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
        function testEffectOfRotationOnMatvec(testCase)
            params = SpheroidalParameters;
            params.isReal = true;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = false;
            params.matvec_eta = 10;
            
            p = 16;
            np = 2*p*(p+1);

            params.sigma = rand(np, 1);

            X_src = params.get_X();
            N_src = params.get_Norm_rot(p);
            X_trg = X_src + N_src;

            unrotated_res = spheroidalMatVec(params, 'DL', X_trg);

            params.thetas = pi/6;
            params.phis = pi/10;
            X_src_rot = params.get_X();
            N_src_rot = params.get_Norm_rot(p);
            X_trg_rot = X_src_rot + N_src_rot;
            rotated_res = spheroidalMatVec(params, 'DL', X_trg_rot);

            testCase.verifyLessThan(norm(unrotated_res - rotated_res) / norm(rotated_res), 1e-14, ...
                "In the case of one spheroid, rotation should have no effect on the matvec.");
        end
    end
end

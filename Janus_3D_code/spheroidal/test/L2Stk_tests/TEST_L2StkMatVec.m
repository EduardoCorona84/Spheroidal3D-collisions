classdef TEST_L2StkMatVec < matlab.unittest.TestCase
    properties
        params=SpheroidalParameters;
        p = 16;

        u0_prolate = 4/sqrt(3);
        u0_oblate = 8/sqrt(3);
        a_prolate;
        a_oblate;
    end

    methods (TestClassSetup)
        function setup(testCase)
            clc();

            testCase.a_prolate = 1/testCase.u0_prolate;
            testCase.a_oblate = 1/sqrt(1 + testCase.u0_oblate^2);
        end
    end

    methods (Test)
        function testSLPMatVec(testCase)
        end

        function testDLPMatVec(testCase)
            p = 16;
            np = 2*p*(p+1);
            params = testCase.params;
            params.matvec_eta = 10; 
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = false; 
            params.centers = [0 0 0];

            sigma_x = rand(np, 1);
            params.sigma = sigma_x; params.get_shc;
            sigma_y = rand(np, 1);
            sigma_z = rand(np, 1);

            X_self = prolate_spheroid_shape(params.p, params.u0, params.a);
            X_trg = X_self + get_norm_vecs(params.p, params.u0, params.oblate);

            [DLPmatvec_X, DLPmatvec_Y, DLPmatvec_Z] = L2StkMatVec(params, 'DLP', sigma_x, sigma_y, sigma_z, X_trg);
            target_pts = { X_trg };
            [L2StkDLP_X, L2StkDLP_Y, L2StkDLP_Z] = L2StkDLP(target_pts, params, sigma_x, sigma_y, sigma_z, 1);

            rel_errs = [
                norm(L2StkDLP_X{1} - DLPmatvec_X) ;
                norm(L2StkDLP_Y{1} - DLPmatvec_Y) ;
                norm(L2StkDLP_Z{1} - DLPmatvec_Z)
            ];

            testCase.verifyLessThan(rel_errs, 1e-10, "Matvec doesn't agree with L2StkDLP.");
        end

        function testTSLMatVec(testCase)
        end
    end
end
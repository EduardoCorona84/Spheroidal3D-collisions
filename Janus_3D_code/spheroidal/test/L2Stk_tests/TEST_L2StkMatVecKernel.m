%{
    This is to verify that the code for L2StkMatVecKernel meets the very 
    basic properties associated with their respective operators. For a more
    rigorous test, see TEST_stokes_bie.m.
%}

classdef TEST_L2StkMatVecKernel < matlab.unittest.TestCase
    properties
        params=SpheroidalParameters;
        p = 16;

        u0_prolate = 4/sqrt(3);
        u0_oblate = 4/sqrt(3);
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
        %% SLP tests
        function testSLPAction(testCase)
            %{
                The action of the SLP matrix should be the same as the
                matvec.
            %}
        end

        %% DLP tests
        function testProlateDLPRigidBodyMotion(testCase)
            %{
                Replicates the rigid body motion test for L2StkDLP. Note
                that this test is a bit different due to the formatting of
                the result we get from the matvec (see the comment in the
                L2StkMatVecKernel about the structure of the matrix).
            %}
            tol = 9e-6;
            p = 16;
            np = 2*p*(p+1);
            params = testCase.params;
            params.p = 16;
            params.matvec_eta = 10;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.thetas = 0;
            params.phis = 0;
            params.oblate = false; 
            params.centers = [0 0 0];

            %%% Matvec computation
            DP = L2StkMatVecKernel(params, 'DLP', p);

            %%%
            %%% Translational motion
            %%%
            sigma_x = 3*ones(np, 1);
            sigma_y = 2*ones(np, 1);
            sigma_z = 1*ones(np, 1);
            sig = [sigma_x ; sigma_y; sigma_z];

            DP_res_vec = DP * sig;
            DP_res_x = DP_res_vec(1:np);
            DP_res_y = DP_res_vec(np+1:2*np);
            DP_res_z = DP_res_vec(2*np+1:3*np);

            rel_errs = [
              norm(DP_res_x + 0.5*sigma_x) / norm(DP_res_x);
              norm(DP_res_y + 0.5*sigma_y) / norm(DP_res_y);
              norm(DP_res_z + 0.5*sigma_z) / norm(DP_res_z);
            ];

            testCase.verifyLessThan(rel_errs, tol, 'Rigid body motion test fails for translational motion.');

            %%%
            %%% Rotational motion
            %%%
            X_src = prolate_spheroid_shape(p, params.u0, params.a);
            omega = rand(1, 3) - 0.5;
            sigma_rot = cross(repmat(omega, np, 1), X_src);
            sigma_x = sigma_rot(:, 1);
            sigma_y = sigma_rot(:, 2);
            sigma_z = sigma_rot(:, 3);
            sig = [sigma_x ; sigma_y; sigma_z];

            DP_res_vec = DP * sig;
            DP_res_x = DP_res_vec(1:np);
            DP_res_y = DP_res_vec(np+1:2*np);
            DP_res_z = DP_res_vec(2*np+1:3*np);

            rel_errs = [
              norm(DP_res_x + 0.5*sigma_x) / norm(DP_res_x);
              norm(DP_res_y + 0.5*sigma_y) / norm(DP_res_y);
              norm(DP_res_z + 0.5*sigma_z) / norm(DP_res_z);
            ];

            testCase.verifyLessThan(rel_errs, tol, 'Rigid body motion test fails for rotational motion.');
        end

        function testOblateDLPRigidBodyMotion(testCase)
            %{
                Replicates the rigid body motion test for L2StkDLP. Note
                that this test is a bit different due to the formatting of
                the result we get from the matvec (see the comment in the
                L2StkMatVecKernel about the structure of the matrix).
            %}
            tol = 9e-6;
            p = 16;
            np = 2*p*(p+1);
            params = testCase.params;
            params.matvec_eta = 10;
            params.u0 = testCase.u0_oblate;
            params.a = testCase.a_oblate;
            params.thetas = 0;
            params.phis = 0;
            params.oblate = false; 
            params.centers = [0 0 0];

            %%% Matvec computation
            DP = L2StkMatVecKernel(params, 'DLP', p);

            %%%
            %%% Translational motion
            %%%
            sigma_x = 3*ones(np, 1);
            sigma_y = 2*ones(np, 1);
            sigma_z = 1*ones(np, 1);
            sig = [sigma_x ; sigma_y; sigma_z];

            DP_res_vec = DP * sig;
            DP_res_x = DP_res_vec(1:np);
            DP_res_y = DP_res_vec(np+1:2*np);
            DP_res_z = DP_res_vec(2*np+1:3*np);

            rel_errs = [
              norm(DP_res_x + 0.5*sigma_x) / norm(DP_res_x);
              norm(DP_res_y + 0.5*sigma_y) / norm(DP_res_y);
              norm(DP_res_z + 0.5*sigma_z) / norm(DP_res_z);
            ];

            testCase.verifyLessThan(rel_errs, tol, 'Rigid body motion test fails for translational motion.');

            %%%
            %%% Rotational motion
            %%%
            X_src = prolate_spheroid_shape(p, params.u0, params.a);
            omega = rand(1, 3) - 0.5;
            sigma_rot = cross(repmat(omega, np, 1), X_src);
            sigma_x = sigma_rot(:, 1);
            sigma_y = sigma_rot(:, 2);
            sigma_z = sigma_rot(:, 3);
            sig = [sigma_x ; sigma_y; sigma_z];

            DP_res_vec = DP * sig;
            DP_res_x = DP_res_vec(1:np);
            DP_res_y = DP_res_vec(np+1:2*np);
            DP_res_z = DP_res_vec(2*np+1:3*np);

            rel_errs = [
              norm(DP_res_x + 0.5*sigma_x) / norm(DP_res_x);
              norm(DP_res_y + 0.5*sigma_y) / norm(DP_res_y);
              norm(DP_res_z + 0.5*sigma_z) / norm(DP_res_z);
            ];

            testCase.verifyLessThan(rel_errs, tol, 'Rigid body motion test fails for rotational motion.');
        end

        %% TLP tests
        function testTLPisAdjointOfDLP(testCase)
            %{
                The traction of the Stokes SLP is the adjoint of the Stokes 
                DLP. We verify this fact.
            %}
            p = 16;
            np = 2*p*(p+1);
            params = testCase.params;
            params.matvec_eta = 10;
            params.u0 = testCase.u0_oblate;
            params.a = testCase.a_oblate;
            params.thetas = 0;
            params.phis = 0;
            params.oblate = false; 
            params.centers = [0 0 0];

            %%% Matrix computation
            DP = L2StkMatVecKernel(params, 'DLP', p);

            sigma_x = 3*ones(np, 1);
            sigma_y = 2*ones(np, 1);
            sigma_z = 1*ones(np, 1);
            sig = [sigma_x ; sigma_y; sigma_z];

            %%% Matvec computation
            DPmatvec_res = L2StkMatVec(params, 'DLP', sigma_x, sigma_y, sigma_z, X_trg)
        end

        function testTLPRigidBodyMotion(testCase)
            %{
                In rigid body motion, the traction of the SLP should be
                zero.
            %}
            tol = 9e-6;
            p = 16;
            np = 2*p*(p+1);
            params = testCase.params;
            params.p = 16;
            params.matvec_eta = 10;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.thetas = 0;
            params.phis = 0;
            params.oblate = false; 
            params.centers = [0 0 0];

            %%% Matvec computation
            TLP = L2StkMatVecKernel(params, 'TLP', p, 0);

            %%%
            %%% Translational motion
            %%%
            u = [ones(np, 1); zeros(2*np, 1)];

            TLP_res_vec = TLP * u;
            TLP_res_x = TLP_res_vec(1:np);
            TLP_res_y = TLP_res_vec(np+1:2*np);
            TLP_res_z = TLP_res_vec(2*np+1:3*np);

            errs = [
                norm(TLP_res_x);
                norm(TLP_res_y);
                norm(TLP_res_z);
            ]

            testCase.verifyLessThan(errs, tol, 'Rigid body motion test fails for translational motion.');

            %%%
            %%% Rotational motion
            %%%
            X_src = prolate_spheroid_shape(p, params.u0, params.a);
            omega = rand(1, 3) - 0.5;
            sigma_rot = cross(repmat(omega, np, 1), X_src);
            sigma_x = sigma_rot(:, 1);
            sigma_y = sigma_rot(:, 2);
            sigma_z = sigma_rot(:, 3);
            sig = [sigma_x ; sigma_y; sigma_z];

            TLP_res_vec = TLP * sig;
            TLP_res_x = TLP_res_vec(1:np);
            TLP_res_y = TLP_res_vec(np+1:2*np);
            TLP_res_z = TLP_res_vec(2*np+1:3*np);

            rel_errs = [
              norm(TLP_res_x + 0.5*sigma_x) / norm(TLP_res_x);
              norm(TLP_res_y + 0.5*sigma_y) / norm(TLP_res_y);
              norm(TLP_res_z + 0.5*sigma_z) / norm(TLP_res_z);
            ];

            testCase.verifyLessThan(rel_errs, tol, 'Rigid body motion test fails for rotational motion.');
        end
    end
end

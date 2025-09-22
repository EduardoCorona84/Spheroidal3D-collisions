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
        function testDLPAction(testCase)
            %{
                The action of the DLP matrix should be the same as the
                matvec.
            %}
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

            [u, v] = gl_grid(p);
            sigma_x = cos(u) .* sin(u).^2 + sin(u) .* cos(v);
            sigma_y = cos(u).^2 + sin(u).*cos(u).*sin(v);
            sigma_z = cos(u) .* sin(u).^2 - sin(u) .* cos(v) + exp(cos(u).^3);

            X_self = prolate_spheroid_shape(p, params.u0, params.a);

            %%% Matvec matrix computation
            np = 2*p*(p+1);
            DLP = L2StkMatVecKernel(params, 'DLP', p);
            DLP_matvecmtx_result = DLP*[sigma_x ; sigma_y; sigma_z];
            DLP_matvec_x = DLP_matvecmtx_result(1:np);
            DLP_matvec_y = DLP_matvecmtx_result(np+1:2*np);
            DLP_matvec_z = DLP_matvecmtx_result(2*np+1:end);

            %%% Matvec comparison
            [DLP_x, DLP_y, DLP_z] = L2StkMatVec(params, 'DLP', sigma_x, sigma_y, sigma_z, X_self);

            rel_errs = [
                norm(DLP_x - DLP_matvec_x) ./ norm(DLP_x);
                norm(DLP_y - DLP_matvec_y) ./ norm(DLP_y);
                norm(DLP_z - DLP_matvec_z) ./ norm(DLP_z)
            ];

            testCase.verifyLessThan(rel_errs, 1e-6, ...
                "Action of DLP matrix does not match the matvec.")
        end

        function testEffectOfRotationOnDLPMatrix(testCase)
            %{
                In this case of only one spheroid, the matrix should
                remain unaffected by any sort of rotation.
            %}
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

            X_self = prolate_spheroid_shape(p, params.u0, params.a);

            DLPunrotated = L2StkMatVecKernel(params, 'DLP', p);

            params.thetas = pi/6;
            params.phis = pi/10;
            DLProtated = L2StkMatVecKernel(params, 'DLP', p);

            testCase.verifyLessThan(norm(DLPunrotated - DLProtated), 1e-14, ...
                "Action of DLP matrix does not match the matvec for unrotated bodies.")
        end

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

        function testDLPMatrixFormationWithMultipleBodies(testCase)
            %{
                Verifies that the action.
            %}
            
        end

        %% TLP tests
        function testTLPAction(testCase)
            %{
                The action of the DLP matrix should be the same as the
                matvec.
            %}
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

            [u, v] = gl_grid(p);
            sigma_x = cos(u) .* sin(u).^2 + sin(u) .* cos(v);
            sigma_y = cos(u).^2 + sin(u).*cos(u).*sin(v);
            sigma_z = cos(u) .* sin(u).^2 - sin(u) .* cos(v) + exp(cos(u).^3);

            X_self = prolate_spheroid_shape(p, params.u0, params.a);
            nu_self = get_norm_vecs(p, params.u0, params.oblate);

            %%% Matvec matrix computation
            np = 2*p*(p+1);
            TLP = L2StkMatVecKernel(params, 'TLP', p, nu_self);
            TLP_matvecmtx_result = TLP*[sigma_x ; sigma_y; sigma_z];
            TLP_matvec_x = TLP_matvecmtx_result(1:np);
            TLP_matvec_y = TLP_matvecmtx_result(np+1:2*np);
            TLP_matvec_z = TLP_matvecmtx_result(2*np+1:end);

            %%% Matvec comparison
            [TLP_x, TLP_y, TLP_z] = L2StkMatVec(params, 'TLP', sigma_x, sigma_y, sigma_z, X_self, nu_self);

            rel_errs = [
                norm(TLP_x - TLP_matvec_x) ./ norm(TLP_x);
                norm(TLP_y - TLP_matvec_y) ./ norm(TLP_y);
                norm(TLP_z - TLP_matvec_z) ./ norm(TLP_z)
            ];

            testCase.verifyLessThan(rel_errs, 1e-10, ...
                "Action of TSL matrix does not match the action of the matvec.");

            dSmtx = kerneldS([], SurfaceSph(X_self));
            kerneldS_result = dSmtx * reshape([sigma_x, sigma_y, sigma_z].', [], 1);
            kerneldS_result = reshape(kerneldS_result,3,[]).';

            rel_dS_errs = [
                norm(TLP_matvec_x - kerneldS_result(:,1)) / norm(kerneldS_result(:,1));
                norm(TLP_matvec_y - kerneldS_result(:,2)) / norm(kerneldS_result(:,2));
                norm(TLP_matvec_z - kerneldS_result(:,3)) / norm(kerneldS_result(:,3));
            ];
            testCase.verifyLessThan(rel_dS_errs, 1e-6, ...
                'L2StkMatVecKernel does not match up with kerneldS to a reasonable tolerance.');
        end
    end
end

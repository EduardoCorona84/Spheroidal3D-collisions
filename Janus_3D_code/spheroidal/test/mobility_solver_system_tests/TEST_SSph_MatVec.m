classdef TEST_SSph_MatVec < matlab.unittest.TestCase
    properties
    end

    methods (TestClassSetup)
        function setup(testCase)
            clc();
        end
    end

    methods (Test)
        %% Self-evaluation checks
        function testMatVecMatchesL2StkSelfEvaluationCode(testCase)
            %{
            The matvec intended for spheroidal mobility code should match 
            the work done in the L2Stk files. Note that SSph_MatVec
            implements the entire BIO, and not just the principal-valued
            integral operator.
            %}
            p = 16;
            np = 2*p*(p+1);

            %% SLP
            Fparams = struct();
            Fparams.parbd.kerd = 3;
            Fparams.parbd.n3 = 1;
            Fparams.parbd.p = p;
            Fparams.parbd.Ct = [0 0 0];
            Fparams.parbd.equ_radii = 1;
            Fparams.parbd.polar_radii = 1/2;
            Fparams.parbd.eps = 1e-2;
            Fparams.parbd.mdist = 1000;
            Fparams.parbd.out = true;
            Fparams.denseMV = true;
            Fparams.tdisc = "euler";

            Fparams_updated = SpheroidalMS_initparams(Fparams);
            parbd = Fparams_updated.parbd;
            parbd.flag_pot = 'SL_Stk_3D';

            [u, v] = gl_grid(p);
            sigma_x = Ynm(3, 3, u, v);
            sigma_y = zeros(np, 1);
            sigma_z = zeros(np, 1);

            % Interleaved format
            V = reshape([sigma_x.'; sigma_y.'; sigma_z.'], [], 1);
            L = [];

            % SSph_MatVec
            Y = SSph_MatVec(V, L, parbd);
            Y_x = Y(1:3:end);
            Y_y = Y(2:3:end);
            Y_z = Y(3:3:end);

            % L2Stk
            [u0, a] = calculate_u0_and_a_from_radii(parbd.shape_type, parbd.equ_radii, parbd.polar_radii);
            sph_params = SpheroidalParameters();
            sph_params.centers = [0 0 0];
            sph_params.u0 = u0;
            sph_params.a = a;
            sph_params.oblate = strcmp('oblate', parbd.shape_type);
            sph_params.matvec_eta = 1000;
            sph_params.Rmat = eye(3);
            
            sph_params.sigma = sigma_x; % Update p.
            X_self = sph_params.get_X;
            [L2Stkmatvec_x, L2Stkmatvec_y, L2Stkmatvec_z] = L2StkMatVec(sph_params, 'SLP', sigma_x, sigma_y, sigma_z, X_self);

            % Compare
            rel_errs = [
                norm(L2Stkmatvec_x - Y_x);
                norm(L2Stkmatvec_y - Y_y);
                norm(L2Stkmatvec_z - Y_z);
            ];

            testCase.verifyLessThan(rel_errs, 1e-14, 'SSph_MatVec produces incorrect results for Stokes SLP.');

            %% DLP
            parbd.flag_pot = 'DL_Stk_3D';
            
            Y = SSph_MatVec(V, L, parbd);
            Y_x = Y(1:3:end);
            Y_y = Y(2:3:end);
            Y_z = Y(3:3:end);

            % L2Stk
            [L2Stkmatvec_x, L2Stkmatvec_y, L2Stkmatvec_z] = L2StkMatVec(sph_params, 'DLP', sigma_x, sigma_y, sigma_z, X_self);

            % Compare
            rel_errs = [
                norm(L2Stkmatvec_x - Y_x - 0.5*sigma_x);
                norm(L2Stkmatvec_y - Y_y - 0.5*sigma_y);
                norm(L2Stkmatvec_z - Y_z - 0.5*sigma_z);
            ];

            testCase.verifyLessThan(rel_errs, 1e-14, 'SSph_MatVec produces incorrect results for Stokes DLP.');

            %% TSL
            parbd.flag_pot = 'TSL_Stk_3D';
            
            Y = SSph_MatVec(V, L, parbd);
            Y_x = Y(1:3:end);
            Y_y = Y(2:3:end);
            Y_z = Y(3:3:end);

            % L2Stk
            nu_self = sph_params.get_Norm;
            [L2Stkmatvec_x, L2Stkmatvec_y, L2Stkmatvec_z] = L2StkMatVec(sph_params, 'TLP', sigma_x, sigma_y, sigma_z, X_self, nu_self);

            % Compare
            rel_errs = [
                norm(L2Stkmatvec_x - Y_x + 0.5*sigma_x);
                norm(L2Stkmatvec_y - Y_y + 0.5*sigma_y);
                norm(L2Stkmatvec_z - Y_z + 0.5*sigma_z);
            ];

            testCase.verifyLessThan(rel_errs, 1e-14, 'SSph_MatVec produces incorrect results for Stokes TSL.');
        end

        function testMatVecMatchesSphereSelfEvaluationCode(testCase)
            %{
            Same test as the one for L2Stk, except to test whether it matches previous work
            done for the sphere.
            %}
        end

        %% Evaluation with multiple bodies
    end
end

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
        function testMatVecMatchesL2StkSelfEvaluation(testCase)
            %{
            The matvec intended for spheroidal mobility code should match 
            the work done in the L2Stk files. Note that SSph_MatVec
            implements the entire BIO, and not just the principal-valued
            integral operator.
            %}
            p = 8;
            np = 2*p*(p+1);

            %% SLP
            Fparams = struct();
            Fparams.parbd.kerd = 3;
            Fparams.parbd.n3 = 1;
            Fparams.parbd.p = p;
            Fparams.parbd.Ct = [15 0 0];
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
            rotation_mtx = RotationMat([2 5 3], 0.1);
            parbd.MRot = {rotation_mtx};

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
            sph_params.centers = [3 0 0];
            sph_params.u0 = u0;
            sph_params.a = a;
            sph_params.oblate = strcmp('oblate', parbd.shape_type);
            sph_params.matvec_eta = 1000;
            sph_params.Rmat = rotation_mtx;
            % sph_params.Rmat = eye(3);
            
            sph_params.sigma = sigma_x; % Update p.
            X_self = sph_params.get_X;
            [L2Stkmatvec_x, L2Stkmatvec_y, L2Stkmatvec_z] = L2StkMatVec(sph_params, 'SLP', sigma_x, sigma_y, sigma_z, X_self);

            % Compare
            rel_errs = [
                norm(L2Stkmatvec_x - Y_x) / norm(L2Stkmatvec_x);
                norm(L2Stkmatvec_y - Y_y) / norm(L2Stkmatvec_y);
                norm(L2Stkmatvec_z - Y_z) / norm(L2Stkmatvec_z);
            ];

            testCase.verifyLessThan(rel_errs, 1e-8, 'SSph_MatVec produces incorrect results for Stokes SLP.');

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
                norm(L2Stkmatvec_x - Y_x - 0.5*sigma_x) / norm(L2Stkmatvec_x);
                norm(L2Stkmatvec_y - Y_y - 0.5*sigma_y) / norm(L2Stkmatvec_y);
                norm(L2Stkmatvec_z - Y_z - 0.5*sigma_z) / norm(L2Stkmatvec_z);
            ];

            testCase.verifyLessThan(rel_errs, 1e-8, 'SSph_MatVec produces incorrect results for Stokes DLP.');

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

            testCase.verifyLessThan(rel_errs, 1e-8, 'SSph_MatVec produces incorrect results for Stokes TSL.');
        end

        function testMatVecMatchesSphereSelfEvaluationCode(testCase)
            %{
            Same test as the one for L2Stk, except to test whether it matches previous work
            done for the sphere.
            %}
            p = 8;
            np = 2*p*(p+1);

            %% SLP
            Fparams = struct();
            Fparams.parbd.kerd = 3;
            Fparams.parbd.n3 = 1;
            Fparams.parbd.p = p;
            Fparams.parbd.Ct = [15 0 0];
            Fparams.parbd.equ_radii = 4;
            Fparams.parbd.polar_radii = 4;
            Fparams.parbd.eps = 1e-2;
            Fparams.parbd.mdist = 1000;
            Fparams.parbd.out = true;
            Fparams.denseMV = true;
            Fparams.tdisc = "euler";

            Fparams_updated = SpheroidalMS_initparams(Fparams);
            parbd = Fparams_updated.parbd;
            parbd.flag_pot = 'SL_Stk_3D';
            % rotation_mtx = RotationMat([2 5 3], 0.1);
            % parbd.MRot = {rotation_mtx};

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

            % Compare sphere self-evaluation with RBS matrices
            Sc = SurfaceSph(Fparams.parbd.equ_radii*shape_gallery(p,''));
            RBSmtx = kernelS([], Sc);
            RBSres = RBSmtx*V;
            RBSres_x = RBSres(1:3:end);
            RBSres_y = RBSres(2:3:end);
            RBSres_z = RBSres(3:3:end);

            % Compare
            rel_errs = [
                norm(RBSres_x - Y_x) / norm(RBSres_x);
                norm(RBSres_y - Y_y) / norm(RBSres_y);
                norm(RBSres_z - Y_z) / norm(RBSres_z);
            ];

            testCase.verifyLessThan(rel_errs, 1e-8, 'SSph_MatVec produces incorrect results for Stokes SLP.');

            %% DLP
            parbd.flag_pot = 'DL_Stk_3D';
            
            Y = SSph_MatVec(V, L, parbd);
            Y_x = Y(1:3:end);
            Y_y = Y(2:3:end);
            Y_z = Y(3:3:end);

            % RBS
            Sc = SurfaceSph(Fparams.parbd.equ_radii*shape_gallery(p,''));
            RBSmtx = kernelD([], Sc);
            RBSres = RBSmtx*V;
            RBSres_x = RBSres(1:3:end);
            RBSres_y = RBSres(2:3:end);
            RBSres_z = RBSres(3:3:end);

            % Compare
            rel_errs = [
                norm(RBSres_x - Y_x) / norm(RBSres_x);
                norm(RBSres_y - Y_y) / norm(RBSres_y);
                norm(RBSres_z - Y_z) / norm(RBSres_z);
            ];

            testCase.verifyLessThan(rel_errs, 1e-8, 'SSph_MatVec produces incorrect results for Stokes DLP.');

            %% TSL
            parbd.flag_pot = 'TSL_Stk_3D';
            
            Y = SSph_MatVec(V, L, parbd);
            Y_x = Y(1:3:end);
            Y_y = Y(2:3:end);
            Y_z = Y(3:3:end);

            % RBS
            Sc = SurfaceSph(Fparams.parbd.equ_radii*shape_gallery(p,''));
            RBSmtx = kerneldS([], Sc);
            RBSres = RBSmtx*V;
            RBSres_x = RBSres(1:3:end);
            RBSres_y = RBSres(2:3:end);
            RBSres_z = RBSres(3:3:end);

            % Compare
            rel_errs = [
                norm(RBSres_x - Y_x) / norm(RBSres_x);
                norm(RBSres_y - Y_y) / norm(RBSres_y);
                norm(RBSres_z - Y_z) / norm(RBSres_z);
            ];

            testCase.verifyLessThan(rel_errs, 1e-8, 'SSph_MatVec produces incorrect results for Stokes DLP.');
        end

        %% Evaluation with multiple bodies
        function testMatVecMatchesL2StkMultipleBodies(testCase)
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
            rotation_mtx = RotationMat([2 5 3], 0.1);
            parbd.MRot = {rotation_mtx}; % FAILS

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
            sph_params.Rmat = rotation_mtx;
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

            testCase.verifyLessThan(rel_errs, 1e-14, ...
                'SSph_MatVec produces incorrect results for Stokes SLP in the context of multiple bodies.');
        end

        function testMatVecMatchesOnVariedShapeTypes(testCase)
            %{
            Tests whether the matvec works for spheroids and sphere. There
            is a lot of boilerplate code here, and that is because we know
            that evaluating from sphere/spheroid. to target points works.
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
            rotation_mtx = RotationMat([2 5 3], 0.1);
            parbd.MRot = {rotation_mtx}; % FAILS

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
            sph_params.Rmat = rotation_mtx;
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

            testCase.verifyLessThan(rel_errs, 1e-14, ...
                'SSph_MatVec produces incorrect results for Stokes SLP in the context of multiple bodies.');
        end
    end
end

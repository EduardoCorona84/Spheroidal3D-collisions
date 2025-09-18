classdef TEST_L2Stk < matlab.unittest.TestCase
    properties
        params=SpheroidalParameters;
        p = 16;

        u0_prolate = 2/sqrt(3);
        u0_oblate = 2/sqrt(3);
        a_prolate;
        a_oblate;
    end

    methods (TestClassSetup)
        function setup(testCase)
            clc();

            testCase.a_prolate = 1/testCase.u0_prolate;
            testCase.a_oblate = 1/sqrt(1 + testCase.u0_oblate^2);

            addpath(genpath('../../.'))
        end
    end

    methods (Test)
        function testProlateInteriorEvaluation(testCase)
            %{
                Tests evaluation of the Stokes SLP inside a prolate
                spheroid.
            %}
            rng(42);
            %%% SETUP
            p = 16;
            np = 2*p*(p+1);
            params = testCase.params;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = 0; 
            params.centers = [0 0 0];

            % Get Cartesian coordinates from the two spheroids
            nu_trg = get_norm_vecs(p, params.u0, params.oblate);
            X_trg = prolate_spheroid_shape(p, params.u0, params.a) - 0.4*nu_trg;

            sigma_x = rand(np,1)+0.5;
            sigma_y = rand(np,1)-0.5;
            sigma_z = rand(np,1);

            target_pts = cell(1, 1);
            target_pts{1} = X_trg;
            [L2Stkx, L2Stky, L2Stkz] = L2Stk(target_pts, params, sigma_x, sigma_y, sigma_z, 1);

            %%% Kernel_Eval
            [X_src, ~] = params.get_X();
            N_src = params.get_Norm(p, 1);

            Sns = SurfaceSph(X_src);
            [~, gwt_gl] = g_grid(p + 1);
            wt_gl = pi/p * repmat(gwt_gl', 2*p, 1) ./ sin(gl_grid(p));
            wt_gl = wt_gl(:);
            W_src = Sns.geoProp.W .* wt_gl;

            pot = 'SL_Stk_3D';
            KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-12,2,400,1);
            KEparams.dim = 3;
            Xv = reshape(repmat(X_src,1,3)',3,[])';
            Wv = repmat(W_src,1,3)'; Wv=Wv(:);
            KEparams.X = Xv;
            KEparams.W2 = Wv.';
            KEparams.nor = reshape(repmat(N_src,1,3)',3,[])';
            KEparams.ci = repmat((1:3)', size(X_trg, 1), 1);
            KEparams.cj = repmat((1:3)', size(X_src, 1), 1);
            KEparams.mu = 1;

            Xtrg_ii = reshape(repmat(X_trg,1,3)',3,[])';
            dDL_mat = Kernel_Eval(Xtrg_ii, Xv, KEparams);

            % Evaluate on density
            sig = reshape([sigma_x,sigma_y,sigma_z].',[],1);
            DP_kernel_eval = dDL_mat * sig;
            DP_kernel_eval = reshape(DP_kernel_eval,3,[]).';

            %%% Compare results
            tol = 1e-6;
            rel_errs = [
                norm(L2Stkx{1} - DP_kernel_eval(:,1)) ./ norm(DP_kernel_eval(:,1));
                norm(L2Stky{1} - DP_kernel_eval(:,2)) ./ norm(DP_kernel_eval(:,2));
                norm(L2Stkz{1} - DP_kernel_eval(:,3)) ./ norm(DP_kernel_eval(:,3))
            ];
            testCase.verifyLessThan(rel_errs, tol, ...
                'Formula for Stokes SLP does not match Kernel_Eval.');
        end

        function testWithKernelS(testCase)
            p = testCase.p;
            np = 2*p*(p+1);
            params = testCase.params;
            params.matvec_eta = 10; 
            params.u0 = testCase.u0_oblate;
            params.a = testCase.a_oblate;
            params.oblate = 1; 
            params.centers = [0 0 0];

            [X_src, ~] = params.get_X();
            
            %%% Translational motion
            sigma_x = 4 * ones(np, 1);
            sigma_y = 2 * ones(np, 1);
            sigma_z = 3 * ones(np, 1);

            [L2Stkx, L2Stky, L2Stkz] = L2Stk([], params, ...
                sigma_x, sigma_y, sigma_z, 1);

            Sc = SurfaceSph(oblate_spheroid_shape(p, params.u0, params.a));
            DMat = kernelS([], Sc);
            sig = reshape([sigma_x,sigma_y,sigma_z].', [], 1);
            DP_res_vec = DMat * sig;
            DP_res_vec = reshape(DP_res_vec,3,[]).';

            norm(DP_res_vec(:,1) - L2Stkx)
        end

        %% Test vector spherical harmonics
        function testActionOfVSHOnSLP(testCase)
            %{
                Verifies the eigenvalue relationship of vector spherical
                harmonics with the Stokes SLP. Note that this emulates
                Test_Spharm_Stk.m.
            %}
            tol = 1e-6;
            p = 16;
            np = 2*p*(p+1);
            params = SpheroidalParameters(); 
            params.u0 = 10000001/sqrt(20000001); % AR = 1+1e-7
            params.a = 1/params.u0;
            params.oblate = 0; 
            params.centers = [0 0 0];
            params.sigma = ones(np, 1);
            params.isReal = false;

            X_src = params.get_X();
            Sc = SurfaceSph(X_src);

            %% Setup n and m
            n = 11; m = 3;
            [u,v]=gl_grid(p);

            %% Vnm
            V = Vnm('Vnm', Sc, n, m, u, v);
            V = reshape(V,3,[]).';

            [L2Stkx, L2Stky, L2Stkz] = L2Stk([], params, V(:,1), V(:,2), V(:,3), 1);
            L2Stk_result = [L2Stkx, L2Stky, L2Stkz];

            Vnm_eigval = n/((2*n+1)*(2*n+3));
            Vnm_rel_err = norm(L2Stk_result - Vnm_eigval*V) / norm(L2Stk_result);
            testCase.verifyLessThan(Vnm_rel_err, tol, ...
                'Spectral relationship failed to meet tolerance for Vnm.');

            %% Wnm
            W = Vnm('Wnm', Sc, n, m, u, v);
            W = reshape(W,3,[]).';

            [L2Stkx, L2Stky, L2Stkz] = L2Stk([], params, W(:,1), W(:,2), W(:,3), 1);
            L2Stk_result = [L2Stkx, L2Stky, L2Stkz];

            Wnm_eigval = (n+1)/((2*n+1)*(2*n-1));
            Wnm_rel_err = norm(L2Stk_result - Wnm_eigval*W) / norm(L2Stk_result);
            testCase.verifyLessThan(Wnm_rel_err, tol, ...
                'Spectral relationship failed to meet tolerance for Wnm.');

            %% Xnm
            X = Vnm('Xnm', Sc, n, m, u, v);
            X = reshape(X,3,[]).';

            [L2Stkx, L2Stky, L2Stkz] = L2Stk([], params, X(:,1), X(:,2), X(:,3), 1);
            L2Stk_result = [L2Stkx, L2Stky, L2Stkz];

            Xnm_eigval = 1/(2*n+1);
            Xnm_rel_err = norm(L2Stk_result - Xnm_eigval*X) / norm(L2Stk_result);
            testCase.verifyLessThan(Xnm_rel_err, tol, ...
                'Spectral relationship failed to meet tolerance for Xnm.');
        end
    end
end

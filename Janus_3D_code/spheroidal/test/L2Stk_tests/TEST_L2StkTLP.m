%{
    Test code for spheroidalDP.m.
%}


classdef TEST_L2StkTLP < matlab.unittest.TestCase
    properties
        params=SpheroidalParameters;
        p = 16;
       
        % u0_prolate = 8/(3*sqrt(7)); % AR = 8
        % u0_prolate = 4/sqrt(15); % AR = 4
        % u0_prolate = 3/(2*sqrt(2)); % AR = 3
        u0_prolate = 2/sqrt(3); % AR = 2
        % u0_prolate = 3/sqrt(5); % AR = 1.5
        % u0_prolate = 13/sqrt(69); % AR = 1.3
        % u0_prolate = 6/sqrt(11); % AR = 1.2
        % u0_prolate = 11/sqrt(21); % AR = 1.1
        % u0_prolate = 3.2796
        % u0_prolate = 21/sqrt(41); % AR = 1.05
        % u0_prolate = 21/sqrt(41) + 1e-11; % AR = 1.05
        % u0_prolate = 100001/sqrt(200001); % AR = 1+1e-5
        u0_oblate = 3/sqrt(5);
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
        %% Kernel_Eval check (off-surface)
        function testTLPProlateKernelEval(testCase)
            %{
                Should agree with the TLP implementation in Kernel_Eval
                off-surface.
            %}
            rng(42);

            p = 16;
            np = 2*p*(p+1);
            params = testCase.params;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = 0; 
            params.centers = [0 0 0];
            params.isReal = true;

            X_trg = prolate_spheroid_shape(p, params.u0, params.a);
            nu_trg = get_norm_vecs(p, params.u0, params.oblate);
            X_trg = X_trg + 5*nu_trg;
            nt = size(X_trg, 1);

            % 
            nu_trg = nu_trg + 5*rand(size(nu_trg));

            [u, v] = gl_grid(p);
            sigma_x = cos(u) .* sin(u).^2 + sin(u) .* cos(v);
            sigma_y = cos(u).^2 + sin(u).*cos(u).*sin(v);
            sigma_z = cos(u) .* sin(u).^2 - sin(u) .* cos(v) + exp(cos(u).^3);

            %%% Stokes Kernel_Eval
            % Get source geometry and weights
            X_src = prolate_spheroid_shape(p, params.u0, params.a);

            % Quadrature weights for Kernel_Eval
            Sns = SurfaceSph(X_src);
            [~, gwt_gl] = g_grid(p + 1);
            wt_gl = pi/p * repmat(gwt_gl', 2*p, 1) ./ sin(gl_grid(p));
            wt_gl = wt_gl(:);
            W_src = Sns.geoProp.W .* wt_gl;
            Wv = repmat(W_src,1,3)'; Wv=Wv(:);

            pot = 'TSL_Stk_3D';
            KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-12,2,400,1);
            KEparams.dim = 3;
            Xv = reshape(repmat(X_src,1,3)',3,[])';
            KEparams.X = Xv;
            KEparams.W2 = Wv.';
            KEparams.nor = reshape(repmat(nu_trg,1,3)',3,[])';
            KEparams.ci = repmat((1:3)', size(X_trg, 1), 1);
            KEparams.cj = repmat((1:3)', size(X_src, 1), 1);

            Xtrg_ii = reshape(repmat(X_trg,1,3)',3,[])';
            TLP_mat = Kernel_Eval(Xtrg_ii, Xv, KEparams);

            %%% Compute results
            sig = reshape([sigma_x,sigma_y,sigma_z].',[],1);
            TLP_kernel_eval = TLP_mat * sig;
            TLP_kernel_eval = reshape(TLP_kernel_eval,3,[]).';

            %%% Spectral calculation
            target_pts = cell(1, 1);
            target_pts{1} = X_trg;
            target_nu = cell(1, 1);
            target_nu{1} = nu_trg;
            [L2StkTLPx, L2StkTLPy, L2StkTLPz] = L2StkTLP(target_pts, target_nu, params, sigma_x, sigma_y, sigma_z, 1, false);

            %%% Compare results
            rel_errs = [
                norm(L2StkTLPx{1} - TLP_kernel_eval(:,1)) ./ norm(TLP_kernel_eval(:,1));
                norm(L2StkTLPy{1} - TLP_kernel_eval(:,2)) ./ norm(TLP_kernel_eval(:,2));
                norm(L2StkTLPz{1} - TLP_kernel_eval(:,3)) ./ norm(TLP_kernel_eval(:,3));
            ];
            tol = 9e-6;
            testCase.verifyLessThan(rel_errs, tol, ...
                'Formula for Stokes TLP does not match Kernel_Eval.');
        end

        function testTLPOblateKernelEval(testCase)
            %{
                Should agree with the TLP implementation in Kernel_Eval
                off-surface.
            %}
            rng(42);

            p = 16;
            np = 2*p*(p+1);
            params = testCase.params;
            params.u0 = testCase.u0_oblate;
            params.a = testCase.a_oblate;
            params.oblate = true; 
            params.centers = [0 0 0];

            X_trg = oblate_spheroid_shape(p, params.u0, params.a);
            nu_trg = get_norm_vecs(p, params.u0, params.oblate);
            X_trg = X_trg + nu_trg;
            nt = size(X_trg, 1);

            [u, v] = gl_grid(p);
            sigma_x = cos(u) .* sin(u).^2 + sin(u) .* cos(v);
            sigma_y = cos(u).^2 + sin(u).*cos(u).*sin(v);
            sigma_z = cos(u) .* sin(u).^2 - sin(u) .* cos(v) + exp(cos(u).^3);

            %%% Stokes Kernel_Eval
            % Get source geometry and weights
            X_src = oblate_spheroid_shape(p, params.u0, params.a);
            N_src = params.get_Norm(p, 1);

            % Quadrature weights for Kernel_Eval
            Sns = SurfaceSph(X_src);
            [~, gwt_gl] = g_grid(p + 1);
            wt_gl = pi/p * repmat(gwt_gl', 2*p, 1) ./ sin(gl_grid(p));
            wt_gl = wt_gl(:);
            W_src_orig = Sns.geoProp.W .* wt_gl;
            Wv = repmat(W_src_orig,1,3)'; Wv=Wv(:);

            pot = 'TSL_Stk_3D';
            KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-12,2,400,1);
            KEparams.dim = 3;
            Xv = reshape(repmat(X_src,1,3)',3,[])';
            KEparams.X = Xv;
            KEparams.W2 = Wv.';
            KEparams.nor = reshape(repmat(nu_trg,1,3)',3,[])';
            KEparams.ci = repmat((1:3)', size(X_trg, 1), 1);
            KEparams.cj = repmat((1:3)', size(X_src, 1), 1);

            Xtrg_ii = reshape(repmat(X_trg,1,3)',3,[])';
            TLP_mat = Kernel_Eval(Xtrg_ii, Xv, KEparams);

            %%% Compute results
            sig = reshape([sigma_x,sigma_y,sigma_z].',[],1);
            TLP_kernel_eval = TLP_mat * sig;
            TLP_kernel_eval = reshape(TLP_kernel_eval,3,[]).';

            %%% Spectral calculation
            target_pts = cell(1, 1);
            target_pts{1} = X_trg;
            target_nu = cell(1, 1);
            target_nu{1} = nu_trg;
            [L2StkTLPx, L2StkTLPy, L2StkTLPz] = L2StkTLP(target_pts, target_nu, params, sigma_x, sigma_y, sigma_z, 1, false);

            %%% Compare results
            rel_errs = [
                norm(L2StkTLPx{1} - TLP_kernel_eval(:,1)) ./ norm(TLP_kernel_eval(:,1));
                norm(L2StkTLPy{1} - TLP_kernel_eval(:,2)) ./ norm(TLP_kernel_eval(:,2));
                norm(L2StkTLPz{1} - TLP_kernel_eval(:,3)) ./ norm(TLP_kernel_eval(:,3));
            ];
            tol = 9e-6;
            testCase.verifyLessThan(rel_errs, tol, ...
                'Formula for Stokes TLP does not match Kernel_Eval.');
        end

        function testTSLProlateWithAnotherSpheroid(testCase)
            %{
                The point of this test is to test whether the off-diagonal
                block matrices are being done correctly.
            %}
            rng(42);
            p = 16;
            np = 2*p*(p+1);
            params = testCase.params;
            params.u0 = [4 2/sqrt(3)];
            params.a = 1./params.u0;
            params.oblate = [0 0];
            params.thetas = [0 0];
            params.phis = [0 0];
            params.centers = [0 0 0; 5 0 0];

            params.sigma = randn(np,1); % Force p=16.

            %% Random density
            % You can achieve better accuracy with an actual smooth density.
            sigma_x = randn(np,1);
            sigma_y = randn(np,1);
            sigma_z = randn(np,1);

            %% Setup points and normals
            [X_trg, X_src] = params.get_X(1);
            [Nu_trg, Nu_src] = params.get_Norm_rot(p, 1);
            Nu_trg = Nu_trg(1:544,:);

            %% Kernel_Eval evaluation
            % Quadrature weights for Kernel_Eval
            Sns = SurfaceSph(X_src);
            [~, gwt_gl] = g_grid(p + 1);
            wt_gl = pi/p * repmat(gwt_gl', 2*p, 1) ./ sin(gl_grid(p));
            wt_gl = wt_gl(:);
            W_src = Sns.geoProp.W .* wt_gl;
            Wv = repmat(W_src,1,3)'; Wv=Wv(:);

            % Kernel_Eval
            KEparams = Kernel_Eval_parameters('TSL_Stk_3D',0,1,1,1,1e-12,2,400,1);
            KEparams.dim = 3;
            Xv = reshape(repmat(X_src,1,3)',3,[])';
            KEparams.X = Xv;
            KEparams.W2 = Wv.';
            KEparams.nor = reshape(repmat(Nu_trg,1,3)',3,[])';
            KEparams.ci = repmat((1:3)', size(X_trg, 1), 1);
            KEparams.cj = repmat((1:3)', size(X_src, 1), 1);

            Xtrg_ii = reshape(repmat(X_trg,1,3)',3,[])';
            TSL_mat = Kernel_Eval(Xtrg_ii, Xv, KEparams);
            Kernel_Eval_result = TSL_mat*reshape([sigma_x,sigma_y,sigma_z].',[],1);
            Kernel_Eval_result = reshape(Kernel_Eval_result,3,[]).';

            %% S'' evaluation
            params_i = SpheroidalParameters;
            params_i.u0 = 4;
            params_i.a = 1/params_i.u0;
            params_i.oblate = 0;
            params_i.thetas = 0;
            params_i.phis = 0;
            params_i.centers = [0 0 0];
            target_pts = cell(1, 1);
            target_pts{1} = X_trg;
            target_nu = cell(1, 1);
            target_nu{1} = Nu_trg;
            [L2StkTLPx, L2StkTLPy, L2StkTLPz] = L2StkTLP(target_pts, target_nu, params_i, sigma_x, sigma_y, sigma_z, 1, false);
        
            %% Verification
            rel_TLP_errs = [
                norm(L2StkTLPx{1} - Kernel_Eval_result(:,1)) ./ norm(Kernel_Eval_result(:,1));
                norm(L2StkTLPy{1} - Kernel_Eval_result(:,2)) ./ norm(Kernel_Eval_result(:,2));
                norm(L2StkTLPz{1} - Kernel_Eval_result(:,3)) ./ norm(Kernel_Eval_result(:,3));
            ];
            testCase.verifyLessThan(rel_TLP_errs, 1e-7, ...
                'Off-surface evaluation does not work for TSL.');

            params_i.sigma = eye(np);
            mtx = TSLmatrix_at_target(params_i, X_trg, Nu_trg, p);
            mtx_result = mtx*[sigma_x; sigma_y; sigma_z];
            mtx_result_x = mtx_result(1:np);
            mtx_result_y = mtx_result(np+1:2*np);
            mtx_result_z = mtx_result(2*np+1:end);
            rel_mat_errs = [
                norm(mtx_result_x - Kernel_Eval_result(:,1)) ./ norm(Kernel_Eval_result(:,1));
                norm(mtx_result_y - Kernel_Eval_result(:,2)) ./ norm(Kernel_Eval_result(:,2));
                norm(mtx_result_z - Kernel_Eval_result(:,3)) ./ norm(Kernel_Eval_result(:,3));
            ];

            testCase.verifyLessThan(rel_mat_errs, 1e-7, ...
                'Matrix evaluation is incorrect.');

            prm = zeros(1,3*np); 
            prm(1:np) = 1:3:3*np; prm(np+1:2*np) = 2:3:3*np; prm(2*np+1:3*np)=3:3:3*np;
            testCase.verifyLessThan(norm(mtx - TSL_mat(prm, prm)), 1e-10, ...
                'Matrix from L2StkMatVecKernel does not match with Kernel_Eval off-surface.');
        end

        %% On-surface check
        function testCheckWithkerneldS(testCase)
            %{
                Compares the on-surface evaluation with the implementation
                in kerneldS.
            %}
            p = 16;
            np = 2*p*(p+1);
            params = SpheroidalParameters;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = false; 
            params.centers = [0 0 0];
            params.isReal = false;

            fine_p = p+4;
            X_self = prolate_spheroid_shape(fine_p, params.u0, params.a);
            nu_self = get_norm_vecs(fine_p, params.u0, params.oblate);

            % Smooth density with low frequencies
            [u, v] = gl_grid(p);
            % sigma_x = cos(u) .* sin(u).^2 + sin(u) .* cos(v);
            % sigma_y = cos(u).^2 + sin(u).*cos(u).*sin(v);
            % sigma_z = cos(u) .* sin(u).^2 - sin(u) .* cos(v) + exp(cos(u).^3);

            % High frequency density
            sigma_x = Ynm(14, 13, u, v);
            sigma_y = Ynm(14, 12, u, v);
            sigma_z = Ynm(15, 15, u, v);

            %%% Spectral calculation
            target_pts = cell(1, 1);
            target_pts{1} = X_self;
            target_nu = cell(1, 1);
            target_nu{1} = nu_self;
            [L2StkTLPx, L2StkTLPy, L2StkTLPz] = L2StkTLP(target_pts, target_nu, params, interpsh(sigma_x, fine_p), interpsh(sigma_y, fine_p), interpsh(sigma_z, fine_p), 1, false);
            L2StkTLPx = L2StkTLPx{1}; L2StkTLPy = L2StkTLPy{1}; L2StkTLPz = L2StkTLPz{1};
            % [L2StkTLPdivx, L2StkTLPdivy, L2StkTLPdivz] = L2StkTLP(target_pts, target_nu, params, sigma_x, sigma_y, sigma_z, 1, true);
            
            L2StkTLPx = interpsh(L2StkTLPx, p);
            L2StkTLPy = interpsh(L2StkTLPy, p);
            L2StkTLPz = interpsh(L2StkTLPz, p);

            %%% kerneldS calculation
            X_self = prolate_spheroid_shape(p, params.u0, params.a);
            Df = kerneldS([], SurfaceSph(X_self));
            kerneldS_result = Df * reshape([sigma_x, sigma_y, sigma_z].', [], 1);
            kerneldS_result = reshape(kerneldS_result,3,[]).';

            rel_errs = [
                norm(L2StkTLPx - kerneldS_result(:,1)) / norm(kerneldS_result(:,1));
                norm(L2StkTLPy - kerneldS_result(:,2)) / norm(kerneldS_result(:,2));
                norm(L2StkTLPz - kerneldS_result(:,3)) / norm(kerneldS_result(:,3));
            ];

            % rel_div_errs = [
            %     norm(L2StkTLPdivx{1} - kerneldS_result(:,1)) / norm(kerneldS_result(:,1));
            %     norm(L2StkTLPdivy{1} - kerneldS_result(:,2)) / norm(kerneldS_result(:,2));
            %     norm(L2StkTLPdivz{1} - kerneldS_result(:,3)) / norm(kerneldS_result(:,3));
            % ];

            tol = 1e-4;
            testCase.verifyLessThan(rel_errs, tol, "On-surface evaluation for the TLP failed.");
        end

        function testSurfaceDivergenceFormulaOffSurface(testCase)
            rng(42);
            p = 16;
            np = 2*p*(p+1);
            params = testCase.params;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = false; 
            params.centers = [0 0 0];

            X_self = prolate_spheroid_shape(p, params.u0, params.a);
            nu_self = get_norm_vecs(p, params.u0, params.oblate);
            X_trg = X_self + 5*nu_self;

            [u, v] = gl_grid(p);
            sigma_x = cos(u) .* sin(u).^2 + sin(u) .* cos(v);
            sigma_y = cos(u).^2 + sin(u).*cos(u).*sin(v);
            sigma_z = cos(u) .* sin(u).^2 - sin(u) .* cos(v) + exp(cos(u).^3);

            %%% Second-order derivative calculation
            target_pts = cell(1, 1);
            target_pts{1} = X_trg;
            target_nu = cell(1, 1);
            target_nu{1} = nu_self;
            [L2StkTLPx, L2StkTLPy, L2StkTLPz] = L2StkTLP(target_pts, target_nu, params, sigma_x, sigma_y, sigma_z, 1, false);

            %%% First-order derivative calculation
            [L2StkTLPx_div, L2StkTLPy_div, L2StkTLPz_div] = L2StkTLP(target_pts, target_nu, params, sigma_x, sigma_y, sigma_z, 1, true);

            %%% Comparison
            rel_errs = [
                norm(L2StkTLPx{1} - L2StkTLPx_div{1}) / norm(L2StkTLPx{1});
                norm(L2StkTLPy{1} - L2StkTLPy_div{1}) / norm(L2StkTLPy{1});
                norm(L2StkTLPz{1} - L2StkTLPz_div{1}) / norm(L2StkTLPz{1});
            ];

            testCase.verifyLessThan(rel_errs, 1e-6, "Divergence formula doesn't match.");
        end

        %% Near-surface (convergence) test
        function testTSLNearSurfaceConvergenceFromInterior(testCase)
            p = 20;
            np = 2*p*(p+1);
            params = testCase.params;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = false; 
            params.centers = [0 0 0];
            params.isReal = true;

            X_self = prolate_spheroid_shape(p, params.u0, params.a);
            nu_self = get_norm_vecs(p, params.u0, params.oblate);

            [u, v] = gl_grid(p);
            sigma_x = cos(u) .* sin(u).^2 + sin(u) .* cos(v);
            sigma_y = cos(u).^2 + sin(u).*cos(u).*sin(v);
            sigma_z = cos(u) .* sin(u).^2 - sin(u) .* cos(v) + exp(cos(u).^3);

            %%% Spectral calculation at surface
            target_pts = cell(1, 1);
            target_pts{1} = X_self;
            target_nu = cell(1, 1);
            target_nu{1} = nu_self;
            [L2StkTLPx_surf, L2StkTLPy_surf, L2StkTLPz_surf] = L2StkTLP(target_pts, target_nu, params, sigma_x, sigma_y, sigma_z, 1, false);

            lim_x = L2StkTLPx_surf{1} + 0.5*sigma_x;
            lim_y = L2StkTLPy_surf{1} + 0.5*sigma_y;
            lim_z = L2StkTLPz_surf{1} + 0.5*sigma_z;

            %%% Spectral calculation away-surface
            distances = 10.^(-2:-1:-5);
            errs_x = [];
            errs_y = [];
            errs_z = [];
            for i=1:numel(distances)
                d = distances(i);
                target_pts = cell(1, 1);
                target_pts{1} = X_self - d*nu_self;
                target_nu = cell(1, 1);
                target_nu{1} = nu_self;
                [L2StkTLPx, L2StkTLPy, L2StkTLPz] = L2StkTLP(target_pts, target_nu, params, sigma_x, sigma_y, sigma_z, 1, false);
                errs_x = [errs_x norm(L2StkTLPx{1} - lim_x) ./ norm(lim_x)];
                errs_y = [errs_y norm(L2StkTLPy{1} - lim_y) ./ norm(lim_y)];
                errs_z = [errs_z norm(L2StkTLPz{1} - lim_z) ./ norm(lim_z)];
            end
        end

        function testTSLNearSurfaceConvergenceFromExterior(testCase)
            p = 16;
            np = 2*p*(p+1);
            params = testCase.params;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = false; 
            params.centers = [0 0 0];
            params.isReal = true;

            X_self = prolate_spheroid_shape(p, params.u0, params.a);
            nu_self = get_norm_vecs(p, params.u0, params.oblate);

            [u, v] = gl_grid(p);
            sigma_x = cos(u) .* sin(u).^2 + sin(u) .* cos(v);
            sigma_y = cos(u).^2 + sin(u).*cos(u).*sin(v);
            sigma_z = cos(u) .* sin(u).^2 - sin(u) .* cos(v) + exp(cos(u).^3);

            %%% Spectral calculation at surface
            target_pts = cell(1, 1);
            target_pts{1} = X_self;
            target_nu = cell(1, 1);
            target_nu{1} = nu_self;
            [L2StkTLPx_surf, L2StkTLPy_surf, L2StkTLPz_surf] = L2StkTLP(target_pts, target_nu, params, sigma_x, sigma_y, sigma_z, 1, false);

            lim_x = L2StkTLPx_surf{1} - 0.5*sigma_x;
            lim_y = L2StkTLPy_surf{1} - 0.5*sigma_y;
            lim_z = L2StkTLPz_surf{1} - 0.5*sigma_z;

            %%% Spectral calculation away-surface
            distances = 10.^(-2:-1:-5);
            errs_x = [];
            errs_y = [];
            errs_z = [];
            for i=1:numel(distances)
                d = distances(i);
                target_pts = cell(1, 1);
                target_pts{1} = X_self + d*nu_self;
                target_nu = cell(1, 1);
                target_nu{1} = nu_self;
                [L2StkTLPx, L2StkTLPy, L2StkTLPz] = L2StkTLP(target_pts, target_nu, params, sigma_x, sigma_y, sigma_z, 1, false);
                errs_x = [errs_x norm(L2StkTLPx{1} - lim_x) ./ norm(lim_x)];
                errs_y = [errs_y norm(L2StkTLPy{1} - lim_y) ./ norm(lim_y)];
                errs_z = [errs_z norm(L2StkTLPz{1} - lim_z) ./ norm(lim_z)];
            end
        end

        %% Vector spherical harmonics test
        function testActionOfVSHOnTSL(testCase)
            %{
                Verifies the eigenvalue relationship of vector spherical
                harmonics with the Stokes TSL on-surface.

                The coefficients for the principal-valued DLP should be the 
                average of the exterior and interior coefficients.

                Note that this emulates Test_Spharm_Stk.m.
            %}
            tol = 1e-6;
            p = 20;
            np = 2*p*(p+1);
            params = SpheroidalParameters(); 
            params.u0 = 10000001/sqrt(20000001); % AR = 1+1e-7
            params.a = 1/params.u0;
            params.oblate = 0; 
            params.centers = [0 0 0];
            params.sigma = ones(np, 1);
            params.isReal = false;

            X_src = params.get_X();
            Nu_src = params.get_Norm();
            Sc = SurfaceSph(X_src);

            %% Setup n and m
            n = 15; m = 0;
            [u,v]=gl_grid(p);

            %% Vnm
            V = Vnm('Vnm', Sc, n, m, u, v);
            V = reshape(V,3,[]).';

            [L2Stkx, L2Stky, L2Stkz] = L2StkTLP({X_src}, {Nu_src}, params, V(:,1), V(:,2), V(:,3), 1, 0);
            L2Stk_result = [L2Stkx{1}, L2Stky{1}, L2Stkz{1}];

            Vnm_eigval = (3/2)/((2*n+1)*(2*n+3));
            Vnm_rel_err = norm(L2Stk_result - Vnm_eigval*V) / norm(L2Stk_result);
            testCase.verifyLessThan(Vnm_rel_err, tol, ...
                'Spectral relationship failed to meet tolerance for Vnm.');

            %% Wnm
            W = Vnm('Wnm', Sc, n, m, u, v);
            W = reshape(W,3,[]).';

            [L2Stkx, L2Stky, L2Stkz] = L2StkTLP({X_src}, {Nu_src}, params, W(:,1), W(:,2), W(:,3), 1, 0);
            L2Stk_result = [L2Stkx{1}, L2Stky{1}, L2Stkz{1}];

            Wnm_eigval = 3/(2 - 8*n^2);
            Wnm_rel_err = norm(L2Stk_result - Wnm_eigval*W) / norm(L2Stk_result);
            testCase.verifyLessThan(Wnm_rel_err, tol, ...
                'Spectral relationship failed to meet tolerance for Wnm.');

            %% Xnm
            X = Vnm('Xnm', Sc, n, m, u, v);
            X = reshape(X,3,[]).';

            [L2Stkx, L2Stky, L2Stkz] = L2StkTLP({X_src}, {Nu_src}, params, X(:,1), X(:,2), X(:,3), 1, 0);
            L2Stk_result = [L2Stkx{1}, L2Stky{1}, L2Stkz{1}];

            Xnm_eigval = -3/(4*n+2);
            Xnm_rel_err = norm(L2Stk_result - Xnm_eigval*X) / norm(L2Stk_result);
            testCase.verifyLessThan(Xnm_rel_err, tol, ...
                'Spectral relationship failed to meet tolerance for Xnm.');
        end

        function testActionOfVSHOnTSLMultipleNM(testCase)
            %{
                Verifies the eigenvalue relationship of vector spherical
                harmonics with the Stokes TSL on-surface.

                The coefficients for the principal-valued DLP should be the 
                average of the exterior and interior coefficients.

                Note that this emulates Test_Spharm_Stk.m.

                For large values of p, this may be very slow. For p=8, this
                takes approximately 3 minutes to finish on my computer. For
                p=16, it takes approximately 1 hour and 30 minutes.
            %}
            tol = 1e-6;
            p = 16;
            np = 2*p*(p+1);
            params = SpheroidalParameters(); 
            params.u0 = 1000001/sqrt(2000001); % AR = 1+1e-7
            params.a = 1/params.u0;
            params.oblate = 0; 
            params.centers = [0 0 0];
            params.sigma = ones(np, 1);
            params.isReal = false;

            X_src = params.get_X();
            Nu_src = params.get_Norm();
            Sc = SurfaceSph(X_src);

            %% Setup n and m
            n = 11; m = 10;
            [u,v]=gl_grid(p);

            sp = (p+1)^2;
            ii = (1:sp)';
            nn = floor(sqrt(ii-1));
            mm = ii - nn.^2 - nn - 1;

            %% Vnm
            Vnm_errors = zeros(size(nn));
            for i=1:sp
                fprintf("Index for Vnm %d\n", i);
                n = nn(i); m = mm(i);
                V = Vnm('Vnm', Sc, n, m, u, v);
                V = reshape(V,3,[]).';
    
                [L2Stkx, L2Stky, L2Stkz] = L2StkTLP({X_src}, {Nu_src}, params, V(:,1), V(:,2), V(:,3), 1, 0);
                L2Stk_result = [L2Stkx{1}, L2Stky{1}, L2Stkz{1}];
    
                Vnm_eigval = (3/2)/((2*n+1)*(2*n+3));
                Vnm_rel_err = norm(L2Stk_result - Vnm_eigval*V) / norm(L2Stk_result);
                Vnm_errors(i) = Vnm_rel_err;
                if Vnm_rel_err >= tol
                    fprintf("Failed to meet tolernace for Vnm with indices (n, m)=(%d, %d) with error %f\n", n, m, Vnm_rel_err);
                end
            end
            testCase.verifyLessThan(Vnm_errors, tol, ...
                    'Spectral relationship failed to meet tolerance for Vnm.');
            plot_errors_for_spherical_eigenvectors(Vnm_errors, p, sp)

            %% Wnm
            Wnm_errors = zeros(size(nn));
            for i=1:sp
                fprintf("Index for Wnm %d\n", i);
                n = nn(i); m = mm(i);
                W = Vnm('Wnm', Sc, n, m, u, v);
                W = reshape(W,3,[]).';
    
                [L2Stkx, L2Stky, L2Stkz] = L2StkTLP({X_src}, {Nu_src}, params, W(:,1), W(:,2), W(:,3), 1, 0);
                L2Stk_result = [L2Stkx{1}, L2Stky{1}, L2Stkz{1}];
    
                Wnm_eigval = 3/(2 - 8*n^2);
                Wnm_rel_err = norm(L2Stk_result - Wnm_eigval*W) / norm(L2Stk_result);
                Wnm_errors(i) = Wnm_rel_err;
                if Wnm_rel_err >= tol
                    fprintf("Failed to meet tolernace for Wnm with indices (n, m)=(%d, %d) with error %f\n", n, m, Vnm_rel_err);
                end
            end
            testCase.verifyLessThan(Wnm_errors, tol, ...
                'Spectral relationship failed to meet tolerance for Wnm.');
            plot_errors_for_spherical_eigenvectors(Wnm_errors, p, sp)

            %% Xnm
            Xnm_errors = zeros(size(nn));
            for i=1:sp
                fprintf("Index for Xnm %d\n", i);
                n = nn(i); m = mm(i);
                X = Vnm('Xnm', Sc, n, m, u, v);
                X = reshape(X,3,[]).';
    
                [L2Stkx, L2Stky, L2Stkz] = L2StkTLP({X_src}, {Nu_src}, params, X(:,1), X(:,2), X(:,3), 1, 0);
                L2Stk_result = [L2Stkx{1}, L2Stky{1}, L2Stkz{1}];
    
                Xnm_eigval = -3/(4*n+2);
                Xnm_rel_err = norm(L2Stk_result - Xnm_eigval*X) / norm(L2Stk_result);
                Xnm_errors(i) = Xnm_rel_err;
                if Xnm_rel_err >= tol
                    fprintf("Failed to meet tolernace for Xnm with indices (n, m)=(%d, %d) with error %f\n", n, m, Vnm_rel_err);
                end
            end
            testCase.verifyLessThan(Xnm_rel_err, tol, ...
                'Spectral relationship failed to meet tolerance for Xnm.');
            plot_errors_for_spherical_eigenvectors(Xnm_errors, p, sp)
        end
    end
end

%% Helper functions (mostly to debug)
function plot_errors_for_spherical_eigenvectors(errors, p, sp)
    % The maximum degree n is p
    n_max = p; 
    
    % Set up the dimensions for the matrix.
    % Columns correspond to n (from 0 to n_max).
    % Rows correspond to m (from -n_max to n_max).
    num_m_values = 2 * n_max + 1;
    num_n_values = n_max + 1;
    
    % Create a matrix to hold the errors, initializing with NaN (Not-a-Number).
    % NaN values will not be plotted, which is perfect for our triangular data.
    error_matrix = nan(num_m_values, num_n_values);
    
    % Loop through your calculated errors and place them in the correct (m, n) position.
    for i = 1:sp
        % Recalculate n and m corresponding to the i-th error value
        n = floor(sqrt(i-1));
        m = i - n^2 - n - 1;
        
        % Determine the matrix indices.
        % Column index is n+1 because MATLAB indices start at 1.
        col_idx = n + 1;
        % Row index needs an offset to handle negative m values.
        row_idx = m + n_max + 1; 
        
        % Assign the error to the correct position in the matrix.
        error_matrix(row_idx, col_idx) = errors(i);
    end
    
    % --- Plotting the result ---
    
    figure;
    % Use imagesc to plot the matrix data as an image.
    % We can use log10 of the errors for better color contrast if they vary wildly.
    imagesc(0:n_max, -n_max:n_max, log10(error_matrix));
    
    % Add labels and a title
    title('log_{10} of V_{n}^{m} relative errors');
    xlabel('n');
    ylabel('m');
    
    % Add a color bar to show the mapping of colors to error values
    h = colorbar;
    ylabel(h, 'log_{10}(relative error)');
    
    % Flip the y-axis to have m=-8 at the bottom and m=8 at the top
    set(gca, 'YDir', 'normal'); 
    % Set axis properties for clarity
    axis tight;
    xticks(0:n_max);
    yticks(-n_max:2:n_max); % Show ticks for every other m value
    grid on;
end
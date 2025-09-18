%{
    Test code for spheroidalDP.m.
%}


classdef TEST_L2StkDLP < matlab.unittest.TestCase
    properties
        params=SpheroidalParameters;
        p = 16;

        u0_prolate = 21/sqrt(41);
        u0_oblate = 3/sqrt(5);
        a_prolate;
        a_oblate;

        tol = 1e-6;
        on_surface_tol = 9*1e-4;
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
        function testLaplaceCheck(testCase)
            %{
                Verifies the spheroidalDL and spheroidalDP works.
            %}
            p = testCase.p;
            np = 2*p*(p+1);
            params = testCase.params;
            params.matvec_eta = 10; 
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = 0; 
            params.centers = [0 0 0];

            % Get Cartesian coordinates from the two spheroids
            Xtrg = prolate_spheroid_shape(p,1.2,1/1.2)+[3,3,1];
            Xself = prolate_spheroid_shape(p, params.u0, params.a);
            
            Sns = SurfaceSph(Xself);
            integrate = @(f) integrateOverS(Sns, f)';
            
            Ntrg = size(Xtrg, 1);
            Nself = size(Xself, 1);
            
            % Pre-calculate Rvec and r for all target-source pairs
            Rvec_x = Xtrg(:, 1)' - Xself(:, 1);
            Rvec_y = Xtrg(:, 2)' - Xself(:, 2);
            Rvec_z = Xtrg(:, 3)' - Xself(:, 3);
            r = sqrt(Rvec_x.^2 + Rvec_y.^2 + Rvec_z.^2);
            
            % Pre-compute the normal vectors at the source points
            norm_vecs = params.get_Norm(p, 1);
            nx_src = norm_vecs(:,1);
            ny_src = norm_vecs(:,2);
            nz_src = norm_vecs(:,3);

            sigma_x = rand(np,1)+0.5;
            sigma_y = rand(np,1)-0.5;
            sigma_z = rand(np,1);

            sigma_list = [sigma_x, sigma_y, sigma_z];
            Imat = eye(3);
            
            fprintf('\n--- Checking Laplace double layer potential and its gradient... ---\n');
            for sigma_ind = 1:3
                % First, compare the implementation of the Laplace double layer
                % potential.
                sig = sigma_list(:, sigma_ind);
                params.sigma = sig;
                params.get_shc;
            
                % Create kernel of Laplace DLP
                DL_sig_integrand = sig .* (Rvec_x .* nx_src + Rvec_y .* ny_src + Rvec_z .* nz_src) ./ (r.^3);
                DL_sig_vec(1, :, sigma_ind) = 1/(4*pi) * integrate(DL_sig_integrand);
            
                % Calculate MatVec for all target points
                DL_matvec_vec(:, 1, sigma_ind) = spheroidalMatVec(params, 'DL', Xtrg);
            
                fprintf("DL sigma_%d vs matvec max abs error: %e\n", sigma_ind, ...
                        max(norm( DL_sig_vec(1,:,sigma_ind)' - DL_matvec_vec(:, 1, sigma_ind) )) );
            
                % Now, check the gradient of the Laplace double layer potential.
                % Direct calculation of e_j \cdot \nabla(DLP) via integration
                % Kernel is grad_x ( (x-y).n_y / |x-y|^3 ) = n_y / r^3 - 3 * (x-y) * ((x-y).n_y) / r^5
                r_cubed = r.^3;
                r_fifth = r.^5;
                RdotN = (Rvec_x .* nx_src + Rvec_y .* ny_src + Rvec_z .* nz_src);
                for j = 1:3
                    if j == 1
                        grad_kernel_j = (nx_src ./ r_cubed) - 3 * Rvec_x .* RdotN ./ r_fifth;
                    elseif j == 2
                        grad_kernel_j = (ny_src ./ r_cubed) - 3 * Rvec_y .* RdotN ./ r_fifth;
                    else % j == 3
                        grad_kernel_j = (nz_src ./ r_cubed) - 3 * Rvec_z .* RdotN ./ r_fifth;
                    end
            
                    dDP_integrand = sig .* grad_kernel_j;
                    dDP_sig_vec(1, :, sigma_ind, j) = 1/(4*pi) * integrate(dDP_integrand);
            
                    % Indirect calculation of the term e_j \cdot \nabla(DLP) via
                    % spheroidalMatVec.
                    nu_stacked_vec = repmat(Imat(j,:), size(Xtrg,1), 1); % Create a normal vector for every target point
                    dDP_matvec_vec(:, 1, sigma_ind, j) = spheroidalMatVec(params, 'DP', Xtrg, nu_stacked_vec);

                    testCase.verifyTrue(max(norm( dDP_sig_vec(1, :, sigma_ind, j)' - dDP_matvec_vec(:, 1, sigma_ind, j) )) < testCase.tol, ...
                        sprintf('gradient of Laplace DP sigma_%d dx%d vs matvec max abs error above tolerance', sigma_ind, j));
                end
            end
            fprintf('------ Done... ------\n');
        end

        function testFormulaWithKernelEval(testCase)
            %{
                Note that this does not test L2StkDLP; instead, this tests
                whether or not the formula works off-surface with
                Kernel_Eval.
            %}
            rng(42);

            %%% Setup
            p = testCase.p;
            np = 2*p*(p+1);
            params = testCase.params;
            params.matvec_eta = 10; 
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = 0; 
            params.centers = [0 0 0];

            X_trg = prolate_spheroid_shape(p, params.u0, params.a);
            nu_trg = get_norm_vecs(p, params.u0, params.oblate);
            X_trg = X_trg + 1e-5*nu_trg;
            nt = size(X_trg, 1);

            sigma_x = rand(np,1)+0.5;
            sigma_y = rand(np,1)-0.5;
            sigma_z = rand(np,1);

            %%% Stokes Kernel_Eval
            % Get source geometry and weights
            X_src = prolate_spheroid_shape(p, params.u0, params.a);
            N_src = params.get_Norm(p, 1);
            normvec_x = N_src(:,1); normvec_y = N_src(:,2); normvec_z = N_src(:,3);

            % Quadrature weights for Kernel_Eval
            Sns = SurfaceSph(X_src);
            [~, gwt_gl] = g_grid(p + 1);
            wt_gl = pi/p * repmat(gwt_gl', 2*p, 1) ./ sin(gl_grid(p));
            wt_gl = wt_gl(:);
            W_src = Sns.geoProp.W .* wt_gl;
            Wv = repmat(W_src,1,3)'; Wv=Wv(:);

            pot = 'DL_Stk_3D';
            KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-12,2,400,1);
            KEparams.dim = 3;
            Xv = reshape(repmat(X_src,1,3)',3,[])';
            KEparams.X = Xv;
            KEparams.W2 = Wv.';
            KEparams.nor = reshape(repmat(N_src,1,3)',3,[])';
            KEparams.ci = repmat((1:3)', size(X_trg, 1), 1);
            KEparams.cj = repmat((1:3)', size(X_src, 1), 1);

            Xtrg_ii = reshape(repmat(X_trg,1,3)',3,[])';
            dDL_mat = Kernel_Eval(Xtrg_ii, Xv, KEparams);

            %%% Laplace Kernel_Evals
            ns = 1;
            % First, we need to evaluate all of the densities.
            y_dot_sig=zeros(size(sigma_x));
            for i=1:ns
                if ~params.oblate(i)
                    Xloc=prolate_spheroid_shape(p, params.u0, params.a);
                else
                    Xloc=oblate_spheroid_shape(p, params.u0, params.a);
                end
                y_dot_sig(:,:,i)=sigma_x(:,:,i).*Xloc(:,1)+sigma_y(:,:,i).*Xloc(:,2)+sigma_z(:,:,i).*Xloc(:,3);
            end

            % Now, we evaluate Kernel_Eval matrices (for now, let's only
            % examine the x-coordinate).
            KEparams_1=Kernel_Eval_parameters('dSL_L_3D',0,1,1,1,1e-8,2,400,1);
            KEparams_1.dim = 3; KEparams_1.mu=1;
            KEparams_1.X=X_src;
            KEparams_1.W2=W_src.';
            Nu_x=repmat([1,0,0],nt,1);
            KEparams_1.nor=Nu_x;
            SPxker=Kernel_Eval(X_trg, X_src, KEparams_1);

            KEparams_2=Kernel_Eval_parameters('dSL_L_3D',0,1,1,1,1e-8,2,400,1);
            KEparams_2.dim = 3; KEparams_2.mu=1;
            KEparams_2.X=X_src;
            KEparams_2.W2=W_src.';
            Nu_y=repmat([0,1,0],nt,1);
            KEparams_2.nor=Nu_y;
            SPyker=Kernel_Eval(X_trg, X_src, KEparams_2);

            KEparams_3=Kernel_Eval_parameters('dSL_L_3D',0,1,1,1,1e-8,2,400,1);
            KEparams_3.dim = 3; KEparams_3.mu=1;
            KEparams_3.X=X_src;
            KEparams_3.W2=W_src.';
            Nu_z=repmat([0,0,1],nt,1);
            KEparams_3.nor=Nu_z;
            SPzker=Kernel_Eval(X_trg, X_src, KEparams_3);

            KEparams_4=Kernel_Eval_parameters('dDL_L_3D',0,1,1,1,1e-8,2,400,1);
            KEparams_4.dim = 3; KEparams_4.mu=1;
            KEparams_4.X=X_src;
            KEparams_4.W2=W_src.';
            KEparams_4.nor=N_src;
            KEparams_4.targnor=Nu_x;
            DPxker=Kernel_Eval(X_trg, X_src, KEparams_4);

            KEparams_5=Kernel_Eval_parameters('dDL_L_3D',0,1,1,1,1e-8,2,400,1);
            KEparams_5.dim = 3; KEparams_5.mu=1;
            KEparams_5.X=X_src;
            KEparams_5.W2=W_src.';
            KEparams_5.nor=N_src;
            KEparams_5.targnor=Nu_y;
            DPyker=Kernel_Eval(X_trg, X_src, KEparams_5);

            KEparams_6=Kernel_Eval_parameters('dDL_L_3D',0,1,1,1,1e-8,2,400,1);
            KEparams_6.dim = 3; KEparams_6.mu=1;
            KEparams_6.X=X_src;
            KEparams_6.W2=W_src.';
            KEparams_6.nor=N_src;
            KEparams_6.targnor=Nu_z;
            DPzker=Kernel_Eval(X_trg, X_src, KEparams_6);
            
            %%% Compute results
            sig = reshape([sigma_x,sigma_y,sigma_z].',[],1);
            DP_kernel_eval = dDL_mat * sig;
            DP_kernel_eval = reshape(DP_kernel_eval,3,[]).';

            first_sum_term = X_trg(:,1).*(DPxker*sigma_x) + X_trg(:,2).*(DPxker*sigma_y) + X_trg(:,3).*(DPxker*sigma_z);
            second_sum_term = SPxker*(normvec_x.*sigma_x) + SPyker*(normvec_x.*sigma_y) + SPzker*(normvec_x.*sigma_z);
            formula_X = -first_sum_term + DPxker*y_dot_sig - second_sum_term;

            first_sum_term = X_trg(:,1).*(DPyker*sigma_x) + X_trg(:,2).*(DPyker*sigma_y) + X_trg(:,3).*(DPyker*sigma_z);
            second_sum_term = SPxker*(normvec_y.*sigma_x) + SPyker*(normvec_y.*sigma_y) + SPzker*(normvec_y.*sigma_z);
            formula_Y = -first_sum_term + DPyker*y_dot_sig - second_sum_term;

            first_sum_term = X_trg(:,1).*(DPzker*sigma_x) + X_trg(:,2).*(DPzker*sigma_y) + X_trg(:,3).*(DPzker*sigma_z);
            second_sum_term = SPxker*(normvec_z.*sigma_x) + SPyker*(normvec_z.*sigma_y) + SPzker*(normvec_z.*sigma_z);
            formula_Z = -first_sum_term + DPzker*y_dot_sig - second_sum_term;
           
            %%% Compare results
            testCase.verifyLessThan(norm(formula_X - DP_kernel_eval(:,1)) ./ norm(formula_X), testCase.tol, ...
                'Formula for Stokes DLP (x-component) does not match Kernel_Eval.');

            testCase.verifyLessThan(norm(formula_Y - DP_kernel_eval(:,2)) ./ norm(formula_Y), testCase.tol, ...
                'Formula for Stokes DLP (y-component) does not match Kernel_Eval.');

            testCase.verifyLessThan(norm(formula_Z - DP_kernel_eval(:,3)) ./ norm(formula_Z), testCase.tol, ...
                'Formula for Stokes DLP (z-component) does not match Kernel_Eval.');
        end

        function testDirectIntegration(testCase)
            %{
                Verifies that the formula works (off-surface).
            %}
            p = testCase.p;
            np = 2*p*(p+1);
            params = testCase.params;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = 0; 
            params.centers = [0 0 0];

            sigma_x = rand(np,1)+0.5;
            sigma_y = rand(np,1)-0.5;
            sigma_z = rand(np,1);

            X_self = prolate_spheroid_shape(p, params.u0, params.a);

            Sns = SurfaceSph(X_self);
            integrate = @(f) integrateOverS(Sns, f)';

            % Pre-compute the normal vectors at the source points
            norm_vecs = params.get_Norm(p, 1);
            nx_src = norm_vecs(:,1);
            ny_src = norm_vecs(:,2);
            nz_src = norm_vecs(:,3);

            nu_trg = get_norm_vecs(p, params.u0, params.oblate);
            X_trg = prolate_spheroid_shape(p, params.u0, params.a) + 2*nu_trg;

            Rvec_x = X_trg(:, 1)' - X_self(:, 1);
            Rvec_y = X_trg(:, 2)' - X_self(:, 2);
            Rvec_z = X_trg(:, 3)' - X_self(:, 3);
            r = sqrt(Rvec_x.^2 + Rvec_y.^2 + Rvec_z.^2);

            %%% Integrate through the formula off-surface.
            r_inv_5 = r .^ (-5);
            RdotN = Rvec_x .* nx_src + Rvec_y .* ny_src + Rvec_z .* nz_src; % (x-y) \cdot n(y)
            RdotSigma = Rvec_x .* sigma_x + Rvec_y .* sigma_y + Rvec_z .* sigma_z; % (x-y) \cdot \sigma(y)
        
            integrand_DLx = RdotN .* Rvec_x .* RdotSigma .* r_inv_5;
            integrand_DLy = RdotN .* Rvec_y .* RdotSigma .* r_inv_5;
            integrand_DLz = RdotN .* Rvec_z .* RdotSigma .* r_inv_5;
        
            % Integrate
            DLx = (3/(4*pi)) * integrate(integrand_DLx);
            DLy = (3/(4*pi)) * integrate(integrand_DLy);
            DLz = (3/(4*pi)) * integrate(integrand_DLz);
        
            % Now, actually calculate result from L2Stk and compare.
            target_pts = cell(1, 1);
            target_pts{1} = X_trg;
            [L2StkDLPx, L2StkDLPy, L2StkDLPz] = L2StkDLP(target_pts, params, sigma_x, sigma_y, sigma_z, 1);

            rel_errs = [
                max(norm(L2StkDLPx{1}-DLx))./norm(DLx)
                max(norm(L2StkDLPy{1}-DLy))./norm(DLy) 
                max(norm(L2StkDLPz{1}-DLz))./norm(DLz)
            ];

            testCase.verifyLessThan(rel_errs(1), testCase.tol, ...
                'Formula for Stokes DLP (x-component) does not match integration formula.');

            testCase.verifyLessThan(rel_errs(2), testCase.tol, ...
                'Formula for Stokes DLP (y-component) does not match integration formula.');

            testCase.verifyLessThan(rel_errs(3), testCase.tol, ...
                'Formula for Stokes DLP (z-component) does not match integration formula.');
        end

        function testKernelEvalProlateOffSurface(testCase)
            %{
                Compares off-surface evaluation with Kernel_Eval's
                implementation and compares error.
            %}
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
            X_trg = prolate_spheroid_shape(p, params.u0, params.a) + 2*nu_trg;

            sigma_x = rand(np,1)+0.5;
            sigma_y = rand(np,1)-0.5;
            sigma_z = rand(np,1);

            target_pts = cell(1, 1);
            target_pts{1} = X_trg;
            [L2Stkx, L2Stky, L2Stkz] = L2StkDLP(target_pts, params, sigma_x, sigma_y, sigma_z, 1);

            %%% Kernel_Eval
            [X_src, ~] = params.get_X();
            N_src = params.get_Norm(p, 1);

            Sns = SurfaceSph(X_src);
            [~, gwt_gl] = g_grid(p + 1);
            wt_gl = pi/p * repmat(gwt_gl', 2*p, 1) ./ sin(gl_grid(p));
            wt_gl = wt_gl(:);
            W_src = Sns.geoProp.W .* wt_gl;

            pot = 'DL_Stk_3D';
            KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-12,2,400,1);
            KEparams.dim = 3;
            Xv = reshape(repmat(X_src,1,3)',3,[])';
            Wv = repmat(W_src,1,3)'; Wv=Wv(:);
            KEparams.X = Xv;
            KEparams.W2 = Wv.';
            KEparams.nor = reshape(repmat(N_src,1,3)',3,[])';
            KEparams.ci = repmat((1:3)', size(X_trg, 1), 1);
            KEparams.cj = repmat((1:3)', size(X_src, 1), 1);

            Xtrg_ii = reshape(repmat(X_trg,1,3)',3,[])';
            dDL_mat = Kernel_Eval(Xtrg_ii, Xv, KEparams);

            % Evaluate on density
            sig = reshape([sigma_x,sigma_y,sigma_z].',[],1);
            DP_kernel_eval = dDL_mat * sig;
            DP_kernel_eval = reshape(DP_kernel_eval,3,[]).';

            %%% Compare results
            testCase.verifyLessThan(norm(L2Stkx{1} - DP_kernel_eval(:,1)) ./ norm(DP_kernel_eval(:,1)), testCase.tol, ...
                'Formula for Stokes DLP (x-component) does not match Kernel_Eval.');

            testCase.verifyLessThan(norm(L2Stky{1} - DP_kernel_eval(:,2)) ./ norm(DP_kernel_eval(:,2)), testCase.tol, ...
                'Formula for Stokes DLP (y-component) does not match Kernel_Eval.');

            testCase.verifyLessThan(norm(L2Stkz{1} - DP_kernel_eval(:,3)) ./ norm(DP_kernel_eval(:,3)), testCase.tol, ...
                'Formula for Stokes DLP (y-component) does not match Kernel_Eval.');
        end

        function testKernelEvalOblateOffSurface(testCase)
            %{
                Compares off-surface evaluation with Kernel_Eval's
                implementation and compares error.

                The oblate case needs a higher order for better
                convergence...
            %}
            %%% SETUP
            p = 20;
            np = 2*p*(p+1);
            params = testCase.params;
            params.u0 = testCase.u0_oblate;
            params.a = testCase.a_oblate;
            params.oblate = true; 
            params.centers = [0 0 0];

            % Get Cartesian coordinates from the two spheroids
            nu_trg = get_norm_vecs(p, params.u0, params.oblate);
            X_trg = oblate_spheroid_shape(p, params.u0, params.a) + 2*nu_trg;

            sigma_x = rand(np,1)+0.5;
            sigma_y = rand(np,1)-0.5;
            sigma_z = rand(np,1);

            target_pts = cell(1, 1);
            target_pts{1} = X_trg;
            [L2Stkx, L2Stky, L2Stkz] = L2StkDLP(target_pts, params, sigma_x, sigma_y, sigma_z, 1);

            %%% Kernel_Eval
            [X_src, ~] = params.get_X();
            N_src = params.get_Norm(p, 1);

            Sns = SurfaceSph(X_src);
            [~, gwt_gl] = g_grid(p + 1);
            wt_gl = pi/p * repmat(gwt_gl', 2*p, 1) ./ sin(gl_grid(p));
            wt_gl = wt_gl(:);
            W_src = Sns.geoProp.W .* wt_gl;

            pot = 'DL_Stk_3D';
            KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-12,2,400,1);
            KEparams.dim = 3;
            Xv = reshape(repmat(X_src,1,3)',3,[])';
            Wv = repmat(W_src,1,3)'; Wv=Wv(:);
            KEparams.X = Xv;
            KEparams.W2 = Wv.';
            KEparams.nor = reshape(repmat(N_src,1,3)',3,[])';
            KEparams.ci = repmat((1:3)', size(X_trg, 1), 1);
            KEparams.cj = repmat((1:3)', size(X_src, 1), 1);

            Xtrg_ii = reshape(repmat(X_trg,1,3)',3,[])';
            dDL_mat = Kernel_Eval(Xtrg_ii, Xv, KEparams);

            % Evaluate on density
            sig = reshape([sigma_x,sigma_y,sigma_z].',[],1);
            DP_kernel_eval = dDL_mat * sig;
            DP_kernel_eval = reshape(DP_kernel_eval,3,[]).';

            %%% Compare results
            tol = 9e-8;
            testCase.verifyLessThan(norm(L2Stkx{1} - DP_kernel_eval(:,1)) ./ norm(DP_kernel_eval(:,1)), tol, ...
                'Formula for Stokes DLP (x-component) does not match Kernel_Eval.');

            testCase.verifyLessThan(norm(L2Stky{1} - DP_kernel_eval(:,2)) ./ norm(DP_kernel_eval(:,2)), tol, ...
                'Formula for Stokes DLP (y-component) does not match Kernel_Eval.');

            testCase.verifyLessThan(norm(L2Stkz{1} - DP_kernel_eval(:,3)) ./ norm(DP_kernel_eval(:,3)), tol, ...
                'Formula for Stokes DLP (y-component) does not match Kernel_Eval.');
        end

        %% Test rotation
        function testKernelEvalProlateOffSurfaceRotatedBody(testCase)
            %%% SETUP
            p = 16;
            np = 2*p*(p+1);
            params = testCase.params;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = 0; 
            params.centers = [0 0 0];
            params.thetas = pi/6;
            params.phis = 0;

            sigma_x = rand(np,1)+0.5;
            sigma_y = rand(np,1)-0.5;
            sigma_z = rand(np,1);

            params.sigma = sigma_x; % A hack to update p.

            [X_src, ~] = params.get_X();
            N_src = params.get_Norm_rot(p);
            X_trg = X_src + 0.3*N_src;
            
            target_pts = cell(1, 1);
            target_pts{1} = X_trg;
            [L2Stkx, L2Stky, L2Stkz] = L2StkDLP(target_pts, params, sigma_x, sigma_y, sigma_z, 1);

            %%% Kernel_Eval
            Sns = SurfaceSph(X_src);
            [~, gwt_gl] = g_grid(p + 1);
            wt_gl = pi/p * repmat(gwt_gl', 2*p, 1) ./ sin(gl_grid(p));
            wt_gl = wt_gl(:);
            W_src = Sns.geoProp.W .* wt_gl;

            pot = 'DL_Stk_3D';
            KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-12,2,400,1);
            KEparams.dim = 3;
            Xv = reshape(repmat(X_src,1,3)',3,[])';
            Wv = repmat(W_src,1,3)'; Wv=Wv(:);
            KEparams.X = Xv;
            KEparams.W2 = Wv.';
            KEparams.nor = reshape(repmat(N_src,1,3)',3,[])';
            KEparams.ci = repmat((1:3)', size(X_trg, 1), 1);
            KEparams.cj = repmat((1:3)', size(X_src, 1), 1);

            Xtrg_ii = reshape(repmat(X_trg,1,3)',3,[])';
            dDL_mat = Kernel_Eval(Xtrg_ii, Xv, KEparams);

            % Evaluate on density
            sig = reshape([sigma_x,sigma_y,sigma_z].',[],1);
            DP_kernel_eval = dDL_mat * sig;
            DP_kernel_eval = reshape(DP_kernel_eval,3,[]).';

            %%% Compare results
            testCase.verifyLessThan(norm(L2Stkx{1} - DP_kernel_eval(:,1)) ./ norm(DP_kernel_eval(:,1)), testCase.tol, ...
                'Formula for Stokes DLP (x-component) does not match Kernel_Eval.');

            testCase.verifyLessThan(norm(L2Stky{1} - DP_kernel_eval(:,2)) ./ norm(DP_kernel_eval(:,2)), testCase.tol, ...
                'Formula for Stokes DLP (y-component) does not match Kernel_Eval.');

            testCase.verifyLessThan(norm(L2Stkz{1} - DP_kernel_eval(:,3)) ./ norm(DP_kernel_eval(:,3)), testCase.tol, ...
                'Formula for Stokes DLP (y-component) does not match Kernel_Eval.');
        end

        %% On-surface checks
        function testProlateRigidBodyMotion(testCase)
            %{
                For rigid body motion (translational + rotational), the 
                following jump relation
                    (-1/2 I + D)[sigma] = sigma
                should hold.
            %}
            p = 16;
            np = 2*p*(p+1);
            params = testCase.params;
            params.matvec_eta = 10; 
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = 0; 
            params.centers = [0 0 0];

            [X_src, ~] = params.get_X();
            
            %%% Translational motion
            sigma_x_trans = 4 * ones(np, 1);
            sigma_y_trans = 2 * ones(np, 1);
            sigma_z_trans = 3 * ones(np, 1);

            [L2Stkx, L2Stky, L2Stkz] = L2StkDLP([], params, ...
                sigma_x_trans, sigma_y_trans, sigma_z_trans, 1);

            Sc = SurfaceSph(prolate_spheroid_shape(p, params.u0, params.a));
            DMat = kernelD([], Sc);
            sig = reshape([sigma_x_trans,sigma_y_trans,sigma_z_trans].', [], 1);
            DP_res_vec = DMat * sig;
            DP_res_vec = reshape(DP_res_vec,3,[]).';

            testCase.verifyLessThan(norm(L2Stkx + 0.5*sigma_x_trans) / norm(L2Stkx), testCase.on_surface_tol, ...
                'Translational RBM identity failed for x-component.');
            testCase.verifyLessThan(norm(L2Stky + 0.5*sigma_y_trans) / norm(L2Stky), testCase.on_surface_tol, ...
                'Translational RBM identity failed for y-component.');
            testCase.verifyLessThan(norm(L2Stkz + 0.5*sigma_z_trans) / norm(L2Stkz), testCase.on_surface_tol, ...
                'Translational RBM identity failed for z-component.');

            %%% Rotational motion
            omega = rand(1, 3) - 0.5;
            sigma_rot = cross(repmat(omega, np, 1), X_src);
            sigma_x_rot = sigma_rot(:, 1);
            sigma_y_rot = sigma_rot(:, 2);
            sigma_z_rot = sigma_rot(:, 3);

            [L2Stkx, L2Stky, L2Stkz] = L2StkDLP([], params, ...
                sigma_x_rot, sigma_y_rot, sigma_z_rot, 1);
            
            testCase.verifyLessThan(norm(L2Stkx + 0.5*sigma_x_rot) / norm(L2Stkx), testCase.on_surface_tol, ...
                'Rotational RBM identity failed for x-component.');
            testCase.verifyLessThan(norm(L2Stky + 0.5*sigma_y_rot) / norm(L2Stky), testCase.on_surface_tol, ...
                'Rotational RBM identity failed for y-component.');
            testCase.verifyLessThan(norm(L2Stkz + 0.5*sigma_z_rot) / norm(L2Stkz), testCase.on_surface_tol, ...
                'Rotational RBM identity failed for z-component.');
        end

        function testOblateRigidBodyMotion(testCase)
            %{
                For rigid body motion (translational + rotational), the 
                following jump relation
                    (-1/2 I + D)[sigma] = sigma
                should hold.
            %}
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
            sigma_x_trans = 4 * ones(np, 1);
            sigma_y_trans = 2 * ones(np, 1);
            sigma_z_trans = 3 * ones(np, 1);

            [L2Stkx, L2Stky, L2Stkz] = L2StkDLP([], params, ...
                sigma_x_trans, sigma_y_trans, sigma_z_trans, 1);

            Sc = SurfaceSph(oblate_spheroid_shape(p, params.u0, params.a));
            DMat = kernelD([], Sc);
            sig = reshape([sigma_x_trans,sigma_y_trans,sigma_z_trans].', [], 1);
            DP_res_vec = DMat * sig;
            DP_res_vec = reshape(DP_res_vec,3,[]).';

            testCase.verifyLessThan(norm(L2Stkx + 0.5*sigma_x_trans) / norm(L2Stkx), testCase.on_surface_tol, ...
                'Translational RBM identity failed for x-component.');
            testCase.verifyLessThan(norm(L2Stky + 0.5*sigma_y_trans) / norm(L2Stky), testCase.on_surface_tol, ...
                'Translational RBM identity failed for y-component.');
            testCase.verifyLessThan(norm(L2Stkz + 0.5*sigma_z_trans) / norm(L2Stkz), testCase.on_surface_tol, ...
                'Translational RBM identity failed for z-component.');

            %%% Rotational motion
            omega = rand(1, 3) - 0.5;
            sigma_rot = cross(repmat(omega, np, 1), X_src);
            sigma_x_rot = 4*ones(np,1) + sigma_rot(:, 1);
            sigma_y_rot = sigma_rot(:, 2);
            sigma_z_rot = sigma_rot(:, 3);

            [L2Stkx, L2Stky, L2Stkz] = L2StkDLP([], params, ...
                sigma_x_rot, sigma_y_rot, sigma_z_rot, 1);
            
            testCase.verifyLessThan(norm(L2Stkx + 0.5*sigma_x_rot) / norm(L2Stkx), testCase.on_surface_tol, ...
                'Rotational RBM identity failed for x-component.');
            testCase.verifyLessThan(norm(L2Stky + 0.5*sigma_y_rot) / norm(L2Stky), testCase.on_surface_tol, ...
                'Rotational RBM identity failed for y-component.');
            testCase.verifyLessThan(norm(L2Stkz + 0.5*sigma_z_rot) / norm(L2Stkz), testCase.on_surface_tol, ...
                'Rotational RBM identity failed for z-component.');
        end

        function ProlateTestWithKernelD(testCase)
            p = 16;
            np = 2*p*(p+1);
            params = testCase.params;
            params.matvec_eta = 10; 
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = false;
            params.centers = [0 0 0];
            
            [u, v] = gl_grid(p);

            factor = 1/np;
            sigma_x = sin(factor*u).*cos(factor*u).^4;
            sigma_y = sin(factor*u).*cos(factor*u).^3;
            sigma_z = sin(factor*u).*cos(factor*u).^2;

            [L2Stkx, L2Stky, L2Stkz] = L2StkDLP([], params, ...
                sigma_x, sigma_y, sigma_z, 1);

            Sc = SurfaceSph(prolate_spheroid_shape(p, params.u0, params.a));
            DMat = kernelD([], Sc);
            sig = reshape([sigma_x,sigma_y,sigma_z].', [], 1);
            DP_res_vec = DMat * sig;
            DP_res_vec = reshape(DP_res_vec,3,[]).';

            rel_errs = [
                norm(DP_res_vec(:,1) - L2Stkx)./norm(L2Stkx) ;
                norm(DP_res_vec(:,2) - L2Stky)./norm(L2Stky) ;
                norm(DP_res_vec(:,3) - L2Stkz)./norm(L2Stkz)
            ];

            testCase.verifyLessThan(rel_errs, 1e-4, "Does not match with kernelD.");
        end

        function OblateTestWithKernelD(testCase)
            p = 16;
            np = 2*p*(p+1);
            params = testCase.params;
            params.matvec_eta = 10; 
            params.u0 = testCase.u0_oblate;
            params.a = testCase.a_oblate;
            params.oblate = 1;
            params.centers = [0 0 0];
            
            [u, v] = gl_grid(p);

            sigma_x = sin(u).^2;
            sigma_y = cos(u).^3;
            sigma_z = sin(u).*cos(2*u).^4;

            [L2Stkx, L2Stky, L2Stkz] = L2StkDLP([], params, ...
                sigma_x, sigma_y, sigma_z, 1);

            Sc = SurfaceSph(oblate_spheroid_shape(p, params.u0, params.a));
            DMat = kernelD([], Sc);
            sig = reshape([sigma_x,sigma_y,sigma_z].', [], 1);
            DP_res_vec = DMat * sig;
            DP_res_vec = reshape(DP_res_vec,3,[]).';

            rel_errs = [
                norm(DP_res_vec(:,1) - L2Stkx)./norm(L2Stkx) ;
                norm(DP_res_vec(:,2) - L2Stky)./norm(L2Stky) ;
                norm(DP_res_vec(:,3) - L2Stkz)./norm(L2Stkz)
            ];

            testCase.verifyLessThan(rel_errs, 1e-4, "Does not match with kernelD.");
        end
        
        %%% Near-surface convergence test
        function testNearSurfaceConvergenceTest(testCase)
            %{
                Analytically testing for near-surface evaluation is hard.
                Instead, we test for convergence. We test that the Stokes
                double layer potential converges to the exterior jump
                relation.
            %}
            p = 16;
            np = 2*p*(p+1);
            params = testCase.params;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = 0; 
            params.centers = [0 0 0];

            [X_src, ~] = params.get_X();
            N_src = params.get_Norm(p, 1);

            % Use a random density
            rng(42); % for reproducibility
            sigma_x = rand(np, 1) - 0.5;
            sigma_y = rand(np, 1) - 0.5;
            sigma_z = rand(np, 1) - 0.5;

            % Calculate the on-surface of the DLP
            [D_pv_x, D_pv_y, D_pv_z] = L2StkDLP([], params, ...
                sigma_x, sigma_y, sigma_z, 1);

            % The limit from the exterior is D[sigma] - 0.5*sigma
            expected_limit_x = D_pv_x - 0.5 * sigma_x;
            expected_limit_y = D_pv_y - 0.5 * sigma_y;
            expected_limit_z = D_pv_z - 0.5 * sigma_z;

            distances = 10.^(-2:-1:-3);
            errors = zeros(length(distances), 1);

            for i = 1:length(distances)
                d = distances(i);
                X_trg = X_src + d * N_src;

                target_pts = cell(1, 1);
                target_pts{1} = X_trg;

                [L2Stkx, L2Stky, L2Stkz] = L2StkDLP(target_pts, params, sigma_x, sigma_y, sigma_z, 1);
                
                total_err = norm([L2Stkx{1} - expected_limit_x; L2Stky{1} - expected_limit_y; L2Stkz{1} - expected_limit_z]);
                total_norm = norm([expected_limit_x; expected_limit_y; expected_limit_z]);
                
                errors(i) = total_err / total_norm;
            end
        end

        %% VSH test on spheres
        function testActionOfVSHOnDLP(testCase)
            %{
                Verifies the eigenvalue relationship of vector spherical
                harmonics with the Stokes DLP on-surface.

                The coefficients for the principal-valued DLP should be the 
                average of the exterior and interior coefficients.

                Note that this emulates Test_Spharm_Stk.m.
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
            n = 4; m = 3;
            [u,v]=gl_grid(p);

            %% Vnm
            V = Vnm('Vnm', Sc, n, m, u, v);
            V = reshape(V,3,[]).';

            [L2Stkx, L2Stky, L2Stkz] = L2StkDLP([], params, V(:,1), V(:,2), V(:,3), 1);
            L2Stk_result = [L2Stkx, L2Stky, L2Stkz];

            Vnm_eigval = (3/2)/((2*n+1)*(2*n+3));
            Vnm_rel_err = norm(L2Stk_result - Vnm_eigval*V) / norm(L2Stk_result);
            testCase.verifyLessThan(Vnm_rel_err, tol, ...
                'Spectral relationship failed to meet tolerance for Vnm.');

            %% Wnm
            W = Vnm('Wnm', Sc, n, m, u, v);
            W = reshape(W,3,[]).';

            [L2Stkx, L2Stky, L2Stkz] = L2StkDLP([], params, W(:,1), W(:,2), W(:,3), 1);
            L2Stk_result = [L2Stkx, L2Stky, L2Stkz];

            Wnm_eigval = 3/(2 - 8*n^2);
            Wnm_rel_err = norm(L2Stk_result - Wnm_eigval*W) / norm(L2Stk_result);
            testCase.verifyLessThan(Wnm_rel_err, tol, ...
                'Spectral relationship failed to meet tolerance for Wnm.');

            %% Xnm
            X = Vnm('Xnm', Sc, n, m, u, v);
            X = reshape(X,3,[]).';

            [L2Stkx, L2Stky, L2Stkz] = L2StkDLP([], params, X(:,1), X(:,2), X(:,3), 1);
            L2Stk_result = [L2Stkx, L2Stky, L2Stkz];

            Xnm_eigval = -3/(4*n+2);
            Xnm_rel_err = norm(L2Stk_result - Xnm_eigval*X) / norm(L2Stk_result);
            testCase.verifyLessThan(Xnm_rel_err, tol, ...
                'Spectral relationship failed to meet tolerance for Xnm.');
        end
    end
end

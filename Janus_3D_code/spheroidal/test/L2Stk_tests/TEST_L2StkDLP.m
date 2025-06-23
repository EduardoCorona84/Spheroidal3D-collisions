%{
    Test code for spheroidalDP.m.
%}


classdef TEST_L2StkDLP < matlab.unittest.TestCase
    properties
        params=SpheroidalParameters;
        p = 8;

        u0_prolate = 2/sqrt(3);
        u0_oblate = 2/sqrt(3);
        a_prolate;
        a_oblate;

        tol = 1e-6;
    end

    methods (TestClassSetup)
        function setup(testCase)
            clc();

            testCase.a_prolate = 1/testCase.u0_prolate;
            testCase.a_oblate = 1/sqrt(1 + testCase.u0_oblate^2);

            testCase.params.matvec_eta=10; 
            testCase.params.u0=1.1;
            testCase.params.a=1/1.1;
            testCase.params.oblate=0; 
            testCase.params.centers=[0 0 0];
        end
    end

    methods (Test)
        function testLaplaceCheck(testCase)
            p = testCase.p;
            np = 2*p*(p+1);
            params = testCase.params;

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
            params.p = p;

            % Get Cartesian coordinates from the two spheroids
            target_u0 = params.u0 * 3;
            X_trg = prolate_spheroid_shape(p, target_u0, params.a);
            nu_trg = get_norm_vecs(p, target_u0, params.oblate);
            nt = size(X_trg, 1);

            sigma_x = rand(np,1)+0.5;
            sigma_y = rand(np,1)-0.5;
            sigma_z = rand(np,1);

            %%% Stokes Kernel_Eval
            % Get source geometry and weights
            X_src_orig = prolate_spheroid_shape(p, params.u0, params.a);
            N_src_orig = params.get_Norm(p, 1);
            normvec_x = N_src_orig(:,1); normvec_y = N_src_orig(:,2); normvec_z = N_src_orig(:,3);

            % Quadrature weights for Kernel_Eval
            Sns = SurfaceSph(X_src_orig);
            [~, gwt_gl] = g_grid(p + 1);
            wt_gl = pi/p * repmat(gwt_gl', 2*p, 1) ./ sin(gl_grid(p));
            wt_gl = wt_gl(:);
            W_src_orig = Sns.geoProp.W .* wt_gl;
            Wv = repmat(W_src_orig,1,3)'; Wv=Wv(:);

            pot = 'DL_Stk_3D';
            KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-12,2,400,1);
            KEparams.dim = 3;
            Xv = reshape(repmat(X_src_orig,1,3)',3,[])';
            KEparams.X = Xv;
            KEparams.W2 = Wv.';
            KEparams.nor = reshape(repmat(N_src_orig,1,3)',3,[])';
            KEparams.ci = repmat((1:3)', size(X_trg, 1), 1);
            KEparams.cj = repmat((1:3)', size(X_src_orig, 1), 1);

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
            KEparams_1.X=X_src_orig;
            KEparams_1.W2=W_src_orig.';
            Nu_x=repmat([1,0,0],nt,1);
            KEparams_1.nor=Nu_x;
            SPxker=Kernel_Eval(X_trg, X_src_orig, KEparams_1);

            KEparams_2=Kernel_Eval_parameters('dSL_L_3D',0,1,1,1,1e-8,2,400,1);
            KEparams_2.dim = 3; KEparams_2.mu=1;
            KEparams_2.X=X_src_orig;
            KEparams_2.W2=W_src_orig.';
            Nu_y=repmat([0,1,0],nt,1);
            KEparams_2.nor=Nu_y;
            SPyker=Kernel_Eval(X_trg, X_src_orig, KEparams_2);

            KEparams_3=Kernel_Eval_parameters('dSL_L_3D',0,1,1,1,1e-8,2,400,1);
            KEparams_3.dim = 3; KEparams_3.mu=1;
            KEparams_3.X=X_src_orig;
            KEparams_3.W2=W_src_orig.';
            Nu_z=repmat([0,0,1],nt,1);
            KEparams_3.nor=Nu_z;
            SPzker=Kernel_Eval(X_trg, X_src_orig, KEparams_3);

            KEparams_4=Kernel_Eval_parameters('dDL_L_3D',0,1,1,1,1e-8,2,400,1);
            KEparams_4.dim = 3; KEparams_4.mu=1;
            KEparams_4.X=X_src_orig;
            KEparams_4.W2=W_src_orig.';
            KEparams_4.nor=N_src_orig;
            KEparams_4.targnor=Nu_x;
            DPxker=Kernel_Eval(X_trg, X_src_orig, KEparams_4);

            KEparams_5=Kernel_Eval_parameters('dDL_L_3D',0,1,1,1,1e-8,2,400,1);
            KEparams_5.dim = 3; KEparams_5.mu=1;
            KEparams_5.X=X_src_orig;
            KEparams_5.W2=W_src_orig.';
            KEparams_5.nor=N_src_orig;
            KEparams_5.targnor=Nu_y;
            DPyker=Kernel_Eval(X_trg, X_src_orig, KEparams_5);

            KEparams_6=Kernel_Eval_parameters('dDL_L_3D',0,1,1,1,1e-8,2,400,1);
            KEparams_6.dim = 3; KEparams_6.mu=1;
            KEparams_6.X=X_src_orig;
            KEparams_6.W2=W_src_orig.';
            KEparams_6.nor=N_src_orig;
            KEparams_6.targnor=Nu_z;
            DPzker=Kernel_Eval(X_trg, X_src_orig, KEparams_6);
            
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

            target_u0 = params.u0 * 3;
            X_trg = prolate_spheroid_shape(p, target_u0, params.a);
            nu_trg = get_norm_vecs(p, target_u0, params.oblate);

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

        function testKernelDMatrixCheck(testCase)
            %{
                Compares with kernelD.m, which produces the matvec for
                the Stokes DLP.
            %}

            %%% SETUP
            rng(42);
            p = testCase.p;
            np = 2*p*(p+1);
            params = testCase.params;

            X_self = prolate_spheroid_shape(p, params.u0, params.a);

            sigma_x = rand(np,1)+0.5;
            sigma_y = rand(np,1)-0.5;
            sigma_z = rand(np,1);

            Sc = SurfaceSph(X_self);
            DMat = kernelD([], Sc);
            [L2Stkx, L2Stky, L2Stkz] = L2StkDLP([], params, sigma_x, sigma_y, sigma_z, 1);

            sig = reshape([sigma_x,sigma_y,sigma_z].',[],1);
            DP_res_vec = DMat * sig;
            DP_res_vec = reshape(DP_res_vec,3,[]).';

            testCase.verifyLessThan(norm((L2Stkx - DP_res_vec(:,1)))./norm(DP_res_vec(:,1)), testCase.tol, ...
                'Formula for Stokes DLP (x-component) fails for on-surface evaluation.');

            testCase.verifyLessThan(norm((L2Stky - DP_res_vec(:,2)))./norm(DP_res_vec(:,2)), testCase.tol, ...
                'Formula for Stokes DLP (y-component) fails for on-surface evaluation.');

            testCase.verifyLessThan(norm((L2Stkz - DP_res_vec(:,3)))./norm(DP_res_vec(:,3)), testCase.tol, ...
                'Formula for Stokes DLP (z-component) fails for on-surface evaluation.');
        end

        function testKernelEvalOffSurface(testCase)
            %{
                Compares off-surface evaluation with Kernel_Eval's
                implementation and compares error.
            %}
            %%% SETUP
            rng(42);
            p = testCase.p;
            np = 2*p*(p+1);
            params = testCase.params;

            % Get Cartesian coordinates from the two spheroids
            target_u0 = params.u0 * 3;
            X_trg = prolate_spheroid_shape(p, target_u0, params.a);
            nu_trg = get_norm_vecs(p, target_u0, params.oblate);

            sigma_x = rand(np,1)+0.5;
            sigma_y = rand(np,1)-0.5;
            sigma_z = rand(np,1);

            target_pts = cell(1, 1);
            target_pts{1} = X_trg;
            [L2Stkx, L2Stky, L2Stkz] = L2StkDLP(target_pts, params, sigma_x, sigma_y, sigma_z, 1);

            %%% Kernel_Eval
            % Get source geometry and weights
            [X_src_orig, ~] = params.get_X(); % Cartesian coordinates of source points
            N_src_orig = params.get_Norm(p, 1);

            % Quadrature weights for Kernel_Eval
            Sns = SurfaceSph(X_src_orig);
            [~, gwt_gl] = g_grid(p + 1);
            wt_gl = pi/p * repmat(gwt_gl', 2*p, 1) ./ sin(gl_grid(p));
            wt_gl = wt_gl(:);
            W_src_orig = Sns.geoProp.W .* wt_gl;

            pot = 'DL_Stk_3D';
            KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-12,2,400,1);
            KEparams.dim = 3;
            Xv = reshape(repmat(X_src_orig,1,3)',3,[])';
            Wv = repmat(W_src_orig,1,3)'; Wv=Wv(:);
            KEparams.X = Xv;
            KEparams.W2 = Wv.';
            KEparams.nor = reshape(repmat(N_src_orig,1,3)',3,[])';
            KEparams.ci = repmat((1:3)', size(X_trg, 1), 1);
            KEparams.cj = repmat((1:3)', size(X_src_orig, 1), 1);

            Xtrg_ii = reshape(repmat(X_trg,1,3)',3,[])';
            dDL_mat = Kernel_Eval(Xtrg_ii, Xv, KEparams);

            % Evaluate on density
            sig = reshape([sigma_x,sigma_y,sigma_z].',[],1);
            DP_kernel_eval = dDL_mat * sig;
            DP_kernel_eval = reshape(DP_kernel_eval,3,[]).'

            %%% Compare results
            testCase.verifyLessThan(norm(L2Stkx{1} -DP_kernel_eval(:,1)) ./ norm(DP_kernel_eval(:,1)), testCase.tol, ...
                'Formula for Stokes DLP (x-component) does not match Kernel_Eval.');

            testCase.verifyLessThan(norm(L2Stky{1} -DP_kernel_eval(:,2)) ./ norm(DP_kernel_eval(:,2)), testCase.tol, ...
                'Formula for Stokes DLP (y-component) does not match Kernel_Eval.');

            testCase.verifyLessThan(norm(L2Stkz{1} -DP_kernel_eval(:,3)) ./ norm(DP_kernel_eval(:,3)), testCase.tol, ...
                'Formula for Stokes DLP (y-component) does not match Kernel_Eval.');
        end
    end
end

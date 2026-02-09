classdef TEST_spheroidalModifiedDLP < matlab.unittest.TestCase
%{
Tests for angular spheroidal wave function ASWFnm.
%}
    methods (Test)  
        function testExteriorExpectedEigenvalue(testCase)
            lambda = 1;
            
            p = 8;
            
            params = SpheroidalParameters;
            params.p = p;
            params.isReal = false;
            params.u0 = 2/sqrt(3);
            params.a = 1/params.u0;
            params.oblate = false;

            gamma = 1j*lambda*params.a;

            [theta, phi] = gl_grid(p);
            n = 5;
            m = 3;
            Snm = ASWFnm(n, m, theta, phi, gamma, 24);

            X_src = params.get_X;
            nu_src = params.get_Norm;
            X_trg = params.get_X + 2*nu_src;

            S = cart2spheroidal(X_trg, params.a, params.oblate);
            u_trg = S(:,1);

            % Quadrature weights for Kernel_Eval
            Sns = SurfaceSph(X_src);
            [~, gwt_gl] = g_grid(p + 1);
            wt_gl = pi/p * repmat(gwt_gl', 2*p, 1) ./ sin(gl_grid(p));
            wt_gl = wt_gl(:);
            W_src_orig = Sns.geoProp.W .* wt_gl;

            pot = 'DL_LMOD_3D';
            KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-12,2,400,1);
            KEparams.dim = 3;
            KEparams.X = X_src;
            KEparams.W2 = W_src_orig.';
            KEparams.nor = nu_src;
            KEparams.targnor = nu_src;
            KEparams.lambda = lambda;

            modDLP_mat = Kernel_Eval(X_trg, X_src, KEparams);
            modDLPres = modDLP_mat * Snm;

            [~, dR1_u0] = Rnm1(n, m, params.u0, gamma);
            [R3_u, ~] = Rnm3(n, m, u_trg, gamma);

            cnm = -1j * gamma * (params.u0^2 - 1);
            eigenvalue = cnm * dR1_u0 .* R3_u;
            diff = modDLPres - eigenvalue .* Snm;
            relErr = norm(diff) / norm(modDLPres);
            testCase.verifyLessThan(relErr, 1e-5, ...
                'Exterior eigenvalue check failed for modified DLP.');
        end

        function testOnSurfaceSpheroidModDLP(testCase)
            % Primary on-surface correctness check:
            % verify spectral result matches analytic eigenvalue action.
            % A direct matrix comparison is retained only as a loose sanity
            % check; convergence with p is tested separately below.
            lambda = 1;
            
            p = 16;
            
            params = SpheroidalParameters;
            params.p = p;
            params.isReal = false;
            params.u0 = 500;
            params.a = 1/params.u0;
            params.oblate = false;

            n = 5;
            m = 3;
            gamma = 1j * lambda * params.a;
            [theta, phi] = gl_grid(p);
            v = cos(theta);
            params.sigma = ASWFnm(n, m, v, phi, gamma, p, 0);
            params.get_shc();

            modDLPspectral = spheroidalModifiedDLP(params, lambda, []);

            [R1, dR1] = Rnm1(n, m, params.u0, gamma);
            [R3, dR3] = Rnm3(n, m, params.u0, gamma);
            cnm = 1j * gamma * (params.u0^2 - 1);
            eigenvalue = 0.5 * cnm * (R1 .* dR3 + R3 .* dR1);
            relErrAnalytic = norm(modDLPspectral - eigenvalue * params.sigma) / ...
                norm(eigenvalue * params.sigma);
            testCase.verifyLessThan(relErrAnalytic, 1e-12, ...
                'On-surface spectral result does not match analytic eigenvalue.');

            S = prolate_spheroid_shape(p, params.u0, params.a);
            [~, ~, modDLPmat] = kernelModifiedLap(SurfaceSph(S), 'DMat', lambda);
            modDLPres = modDLPmat*params.sigma;

            relErrDirect = norm(modDLPspectral - modDLPres) / norm(modDLPres);
            testCase.verifyLessThan(relErrDirect, 1e-4, ...
                'On-surface spectral vs RBS DLP failed.');
        end

        function testOnSurfaceDirectComparisonConvergesWithP(testCase)
            % Examines convergence as p increases.
            lambda = 1;
            u0 = 500;
            n = 5;
            m = 3;
            p_list = [8, 12, 16, 20];

            num_p = numel(p_list);
            rel_err_analytic = zeros(num_p, 1);
            rel_err_RBS = zeros(num_p, 1);

            for k = 1:num_p
                p = p_list(k);
                [rel_err_analytic(k), rel_err_RBS(k)] = ...
                    LOCAL_on_surface_error(p, lambda, u0, n, m);
            end

            T = table(p_list(:), rel_err_analytic, rel_err_RBS, ...
                'VariableNames', {'p', 'relErrAnalytic', 'relErrDirect'});
            disp(T);

            testCase.verifyTrue(all(rel_err_analytic < 1e-11), ...
                'Analytic on-surface analytic error large for some p.');

            testCase.verifyLessThan(rel_err_RBS(end), 1e-4, ...
                'Highest-p direct-vs-spectral error is too large.');
        end
    end
end

function [rel_err_analytic, rel_err_RBS] = LOCAL_on_surface_error(p, lambda, u0, n, m)
    params = SpheroidalParameters;
    params.p = p;
    params.isReal = false;
    params.u0 = u0;
    params.a = 1 / params.u0;
    params.oblate = false;
    
    gamma = 1j * lambda * params.a;
    [theta, phi] = gl_grid(p);
    v = cos(theta);
    params.sigma = ASWFnm(n, m, v, phi, gamma, p, 0);
    params.get_shc();
    
    modDLPspectral = spheroidalModifiedDLP(params, lambda, []);
    
    [R1, dR1] = Rnm1(n, m, params.u0, gamma);
    [R3, dR3] = Rnm3(n, m, params.u0, gamma);
    cnm = 1j * gamma * (params.u0^2 - 1);
    eigenvalue = 0.5 * cnm * (R1 .* dR3 + R3 .* dR1);
    rel_err_analytic = norm(modDLPspectral - eigenvalue * params.sigma) / ...
        norm(eigenvalue * params.sigma);
    
    S = prolate_spheroid_shape(p, params.u0, params.a);
    [~, ~, modDLPmat] = kernelModifiedLap(SurfaceSph(S), 'DMat', lambda);
    modDLPres = modDLPmat * params.sigma;
    rel_err_RBS = norm(modDLPspectral - modDLPres) / norm(modDLPres);
end

%{
    Scratchwork to develop tests for S''.
%}

classdef TEST_spheroidalgraddivSL < matlab.unittest.TestCase
    properties
        u0_prolate = 2/sqrt(3);
        u0_oblate = 8/sqrt(3);
        a_prolate;
        a_oblate;
        tol = 1e-6;
    end

    methods (TestClassSetup)
        function setup(testCase)
            clc();
            testCase.a_prolate = 1/testCase.u0_prolate;
            testCase.a_oblate = 1/sqrt(1 + testCase.u0_oblate^2);
        end
    end

    methods (Test)
        %%%
        %%% Spectral coefficient tests
        %%%
        function testProlateSpectralCoefficients(testCase)
            %{
                Examines the off-surface spectral coefficients of 
                spheroidalgraddivSL and compares it with the known formulas.
            %}
            geti = @(n,m) m+n^2+n+1; % Map (n,m) to 0 <= k <= sp
            p = 16;
            np = 2*p*(p+1);
            sp = (p+1)^2;
            
            params = SpheroidalParameters;
            params.p = p;
            params.isReal = false; % Must be false!
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = false;

            nu_trg = get_norm_vecs(p, params.u0, params.oblate);
            X_trg = prolate_spheroid_shape(p, params.u0, params.a) + 3*nu_trg;

            % Create a matrix where each column is a spherical harmonic 
            % basis function evaluated at the grid points.
            S = cart2spheroidal(X_trg, params.a, params.oblate);
            u_x = S(:,1); v_x = S(:,2); phi_x = S(:,3);
            Y = zeros(np, sp);
            ii = (1:sp)'; 
            nn = floor(sqrt(ii-1)); 
            mm = ii - nn.^2 - nn - 1;
            for k = 1:sp
                n = nn(k); m = mm(k);
                Y(:,k) = Ynm(n,m,real(acos(v_x)),phi_x);
            end

            params.sigma = Y; params.get_shc;
            shc_x = params.sigma_coefficients;
            u0 = params.u0;
            oblate = params.oblate;
            all_u0s=unique(u0);
            Gshc_x = zeros(size(shc_x));
            for i=1:length(all_u0s)
                this_u0=all_u0s(i);
            
                %Calculate basis transformation matrix G: Ynm -> Ynm / sqrt(u0^2-v^2)
                %or /sqrt(u0^2+v^2) for oblate
                this_G_pro = Gmatrix(p,this_u0,0,0);
                this_G_obl = Gmatrix(p,this_u0,0,1);
            
                % keep track of indices and shc of particles with the same u0
                this_index = find(u0==this_u0); 
                obl_index = nonzeros(this_index.*oblate(this_index));
                pro_index = nonzeros(this_index.*~oblate(this_index));
                obl_shc = reshape(shc_x(:,:,obl_index),sp,[],1);
                pro_shc = reshape(shc_x(:,:,pro_index),sp,[],1);
            
                %Apply basis transformation to spheroidal harmonic coefficients
                obl_Gshc = this_G_obl\obl_shc;
                pro_Gshc = this_G_pro\pro_shc;
                
                Gshc_x(:,:,obl_index) = reshape(obl_Gshc,sp,[],length(obl_index));
                Gshc_x(:,:,pro_index) = reshape(pro_Gshc,sp,[],length(pro_index));
            end

            G = Gmatrix(p, params.u0, 0, 0);
            zero_density = zeros(size(Y));
            [graddivSL_Ucomponent, ~, ~] = spheroidalgraddivSL(params, Y, zero_density, zero_density, X_trg);
            graddivSL_coeffs = shAna(graddivSL_Ucomponent);

            ii = (1:sp)';
            nn = floor(sqrt(ii-1));
            mm = ii - nn.^2 - nn - 1;
            bnm = params.a * factorial(nn-mm)./factorial(nn+mm) .* ((-1) .^ (mm)) .* sqrt(params.u0.^2 - 1);
            L = legendre_otc(p, params.u0, 1, 1, 1);
            gnm = L{1}; % P(u0)
            coefficient = bnm .* gnm ./ (params.a.^2); % Common coeffs
            nf = size(params.sigma_coefficients, 2);
            coefficient_mtx = repmat(coefficient,1,nf,1);
            expected_coeffs = coefficient_mtx .* Gshc_x;

            Fr_coeff = -1*(2*u_x.^6 + v_x.^2 + v_x.^4 + (u_x.^2).*(2 - 4*v_x.^2) + (u_x.^4).*(v_x.^2 - 3))./((u_x.^2 - 1).*(u_x - v_x).^3.*(u_x + v_x).^3);
            Fp_coeff = 2.* u_x ./ ((u_x - v_x) .* (u_x + v_x));
            Fpp_coeff = (u_x.^2 - 1) ./ ((u_x - v_x) .* (u_x + v_x));

            [Fr, Fp, Fpp] = solid_harmonic_prime(params.p, params.u0, u_x, params.oblate);

            errors = size(sp, 1);
            F = Fr_coeff.*Fr + Fp_coeff.*Fp + Fpp_coeff.*Fpp;
            expected_graddivSL = (F.*Y)*expected_coeffs;
            for k=1:sp
                n = nn(k); m = mm(k);
                expected_coefficient = shAna(expected_graddivSL(:,geti(n,m)));
                errors(k) = norm(expected_coefficient - graddivSL_coeffs(:,geti(n,m)));
            end

            testCase.verifyLessThan(errors, 9e-12, "On-surface spectral coefficients for the prolate case do not match the expected values to a reasonable tolernace");
        end

        function testOblateSpectralCoefficients(testCase)
            %{
                Examines the off-surface spectral coefficients of 
                spheroidalgraddivSL and compares it with the known formulas.
            %}
            geti = @(n,m) m+n^2+n+1; % Map (n,m) to 0 <= k <= sp
            p = 8;
            np = 2*p*(p+1);
            sp = (p+1)^2;
            
            params = SpheroidalParameters;
            params.p = p;
            params.isReal = false; % Must be false!
            params.u0 = testCase.u0_oblate;
            params.a = testCase.a_oblate;
            params.oblate = true;

            nu_trg = get_norm_vecs(p, params.u0, params.oblate);
            X_trg = oblate_spheroid_shape(p, params.u0, params.a) + 2*nu_trg;

            % Create a matrix where each column is a spherical harmonic 
            % basis function evaluated at the grid points.
            S = cart2spheroidal(X_trg, params.a, params.oblate);
            u_x = S(:,1); v_x = S(:,2); phi_x = S(:,3);
            Yn = zeros(np, sp); Yn1 = zeros(np, sp);
            ii = (1:sp)'; 
            nn = floor(sqrt(ii-1)); 
            mm = ii - nn.^2 - nn - 1;
            for k = 1:sp
                n = nn(k); m = mm(k);
                Yn(:,k) = Ynm(n,m,real(acos(v_x)),phi_x);
                Yn1(:,k) = Ynm(n+1,m,real(acos(v_x)),phi_x);
            end

            params.sigma = Yn; params.get_shc;
            shc_x = params.sigma_coefficients;
            u0 = params.u0;
            oblate = params.oblate;
            all_u0s=unique(u0);
            Gshc_x = zeros(size(shc_x));
            for i=1:length(all_u0s)
                this_u0=all_u0s(i);
            
                %Calculate basis transformation matrix G: Ynm -> Ynm / sqrt(u0^2-v^2)
                %or /sqrt(u0^2+v^2) for oblate
                this_G_pro = Gmatrix(p,this_u0,0,0);
                this_G_obl = Gmatrix(p,this_u0,0,1);
            
                % keep track of indices and shc of particles with the same u0
                this_index = find(u0==this_u0); 
                obl_index = nonzeros(this_index.*oblate(this_index));
                pro_index = nonzeros(this_index.*~oblate(this_index));
                obl_shc = reshape(shc_x(:,:,obl_index),sp,[],1);
                pro_shc = reshape(shc_x(:,:,pro_index),sp,[],1);
            
                %Apply basis transformation to spheroidal harmonic coefficients
                obl_Gshc = this_G_obl\obl_shc;
                pro_Gshc = this_G_pro\pro_shc;
                
                Gshc_x(:,:,obl_index) = reshape(obl_Gshc,sp,[],length(obl_index));
                Gshc_x(:,:,pro_index) = reshape(pro_Gshc,sp,[],length(pro_index));
            end

            G = Gmatrix(p, params.u0, 0, 0);
            zero_density = zeros(size(Yn));
            [~, ~, graddivSL_PHIcomponent] = spheroidalgraddivSL(params, Yn, zero_density, zero_density, X_trg);
            graddivSL_coeffs = shAna(graddivSL_PHIcomponent);

            ii = (1:sp)';
            nn = floor(sqrt(ii-1));
            mm = ii - nn.^2 - nn - 1;
            cnm = 1j .* params.a * factorial(nn-mm)./factorial(nn+mm) .* ((-1) .^ (mm)) .* sqrt(params.u0.^2 + 1);
            L = legendre_otc(p, 1j.*params.u0, 1, 2, 2);
            gnm = L{1}; % P(u0)
            coefficient = cnm .* gnm ./ (params.a.^2); % Common coeffs
            nf = size(params.sigma_coefficients, 2);
            coefficient_mtx = repmat(coefficient,1,nf,1);
            expected_coeffs = coefficient_mtx .* Gshc_x;

            n = nn'; m = mm';
            YnFr_coeff = ((-1).*(1+u_x.^2).*((-1)+v_x.^2)).^(-3/2).*(u_x.^2+v_x.^2).^(-1).*(sqrt(-1).*m.*((-1).*n.*v_x.^2.*((1+u_x.^2).*(1+(-1).*v_x.^2)).^(1/2)+2.*((-1).*(1+u_x.^2).^(-1).*((-1)+v_x.^2).^(-1)).^(1/2).*((-1)+v_x.^2).*(u_x.^2+v_x.^2)+(-1).*u_x.^2.*((-1).*((1+u_x.^2).*(1+(-1).*v_x.^2)).^(1/2)+2.*((-1).*(1+u_x.^2).^(-1).*((-1)+v_x.^2).^(-1)).^(1/2).*(u_x.^2+v_x.^2)+v_x.^2.*(n.*((1+u_x.^2).*(1+(-1).*v_x.^2)).^(1/2)+((-1).*(1+u_x.^2).*((-1)+v_x.^2)).^(1/2)+(-2).*((-1).*(1+u_x.^2).^(-1).*((-1)+v_x.^2).^(-1)).^(1/2).*(u_x.^2+v_x.^2)))).*cos(phi_x)+(n.*v_x.^2.*((1+u_x.^2).*(1+(-1).*v_x.^2)).^(1/2)+(-1).*(1+m.^2).*((-1).*(1+u_x.^2).^(-1).*((-1)+v_x.^2).^(-1)).^(1/2).*((-1)+v_x.^2).*(u_x.^2+v_x.^2)+u_x.^2.*((-1).*((1+u_x.^2).*(1+(-1).*v_x.^2)).^(1/2)+((-1).*(1+u_x.^2).^(-1).*((-1)+v_x.^2).^(-1)).^(1/2).*(u_x.^2+v_x.^2)+m.^2.*((-1).*(1+u_x.^2).^(-1).*((-1)+v_x.^2).^(-1)).^(1/2).*(u_x.^2+v_x.^2)+v_x.^2.*((-1).*(1+m.^2).*u_x.^2.*((-1).*(1+u_x.^2).^(-1).*((-1)+v_x.^2).^(-1)).^(1/2)+(-1).*(1+m.^2).*v_x.^2.*((-1).*(1+u_x.^2).^(-1).*((-1)+v_x.^2).^(-1)).^(1/2)+(1+n).*((-1).*(1+u_x.^2).*((-1)+v_x.^2)).^(1/2)))).*sin(phi_x));
            YnFp_coeff = (-1).*u_x.*(u_x.^2+v_x.^2).^(-1).*(m.*cos(phi_x)+sqrt(-1).*sin(phi_x));
            YnFpp_coeff = zeros(size(YnFp_coeff));
            ((-1)+m+(-1).*n).*v_x.*((-1).*((-1)+m+(-1).*n).^(-1).*(1+m+n).*(1+2.*n).*(3+2.*n).^(-1).*((-1)+v_x.^2).^(-2)).^(1/2).*(u_x.^2+v_x.^2).^(-1).*((sqrt(-1)*(-1)).*m.*cos(phi_x)+sin(phi_x));
            Yn1Fr_coeff = sqrt(-1).*((-1)+m+(-1).*n).*v_x.*((-1)+v_x.^2).^(-1).*((-1).*((-1)+m+(-1).*n).^(-1).*(1+m+n).*(1+2.*n).*(3+2.*n).^(-1).*(u_x.^2+v_x.^2).^(-2)).^(1/2).*(m.*cos(phi_x)+sqrt(-1).*sin(phi_x));
            Yn1Fp_coeff = zeros(size(YnFp_coeff));
            Yn1Fpp_coeff = zeros(size(YnFp_coeff));

            [Fr, Fp, Fpp] = solid_harmonic_prime(params.p, params.u0, u_x, params.oblate);

            n=1; m=2;
            solid_harmonic_Ynm = YnFr_coeff.*Fr + YnFp_coeff.*Fp + YnFpp_coeff.*Fpp;
            solid_harmonic_Yn1m = Yn1Fr_coeff.*Fr + Yn1Fp_coeff.*Fp + Yn1Fpp_coeff.*Fpp;
            expected_graddiv = solid_harmonic_Ynm * expected_coeffs + solid_harmonic_Yn1m * expected_coeffs;
            expected_coefficient = shAna(expected_graddiv(:,geti(n,m)));
            norm(expected_coefficient - graddivSL_coeffs(:,geti(n,m)))

            return;

            errors = size(sp, 1);
            F = YnFr_coeff.*Fr + YnFp_coeff.*Fp + YnFpp_coeff.*Fpp;
            expected_graddivSL = (F.*Yn)*expected_coeffs;
            for k=1:sp
                n = nn(k); m = mm(k);
                expected_coefficient = shAna(expected_graddivSL(:,geti(n,m)));
                errors(k) = norm(expected_coefficient - graddivSL_coeffs(:,geti(n,m)));
            end

            testCase.verifyLessThan(errors, 9e-12, "On-surface spectral coefficients for the prolate case do not match the expected values to a reasonable tolernace");
        end

        %%% Gradient checks
        function testCartesianProlateGradientCheck(testCase)
            %{
                Converts the resulting vector we get from
                spheroidalgraddivSL to spheroidal coordinates and then use
                a finite difference check. This is faster than the
                spheroidal implementation of the test.
            %}
            eps = 10.^(-2:-1:-5);
            n_eps = numel(eps);
            p = 8;
            np = 2*p*(p+1);
            params = SpheroidalParameters;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = 0;
            params.centers = [0 0 0];
            params.isReal = true;

            rng(42);
            sigma_x = rand(np,1)+0.5;
            sigma_y = zeros(np,1); sigma_z = sigma_y;
            densities = { sigma_x, sigma_y, sigma_z };

            %%%
            %%% EXTERIOR CHECK
            %%%

            % Setup target points away from surface
            nu_trg = get_norm_vecs(p, params.u0, params.oblate);
            X_trg = prolate_spheroid_shape(p, params.u0, params.a) + 3*nu_trg;

            % Evaluate spheroidal grad div SLP
            [graddivSL_Ucomponent, graddivSL_Vcomponent, graddivSL_PHIcomponent] = spheroidalgraddivSL(params, sigma_x, sigma_y, sigma_z, X_trg);
            
            S_trg = cart2spheroidal(X_trg, params.a, params.oblate);
            u = S_trg(:,1); v = S_trg(:,2); phi = S_trg(:,3);
            [graddivSL_Xcomponent, graddivSL_Ycomponent, graddivSL_Zcomponent] = convert_from_prolate_basis(graddivSL_Ucomponent, graddivSL_Vcomponent, graddivSL_PHIcomponent, ...
                                                                                                                u, v, phi);
            %%% Calculate divergence of Laplace SLP
            nu_x = repmat([1,0,0],size(X_trg,1),1);
            nu_y = repmat([0,1,0],size(X_trg,1),1);
            nu_z = repmat([0,0,1],size(X_trg,1),1);
            nu = { nu_x, nu_y, nu_z };
            function sum = divSLP(X)
                sum = 0;
                for j = 1:3
                    params.sigma = densities{j}; params.get_shc;
                    sum = sum + spheroidalSP(params, X, nu{j});
                end
            end

            rel_errs_x = zeros(1, n_eps);
            rel_errs_y = zeros(1, n_eps);
            rel_errs_z = zeros(1, n_eps);
            for i=1:n_eps
                epsilon = eps(i);
            
                % Calculate Cartesian gradient approximation
                findif_x = (divSLP(X_trg + epsilon*nu_x) - divSLP(X_trg - epsilon*nu_x))./(2*epsilon);
                findif_y = (divSLP(X_trg + epsilon*nu_y) - divSLP(X_trg - epsilon*nu_y))./(2*epsilon);
                findif_z = (divSLP(X_trg + epsilon*nu_z) - divSLP(X_trg - epsilon*nu_z))./(2*epsilon);

                rel_errs_x(i) = norm(findif_x - graddivSL_Xcomponent) / norm(findif_x);
                rel_errs_y(i) = norm(findif_y - graddivSL_Ycomponent) / norm(findif_y);
                rel_errs_z(i) = norm(findif_z - graddivSL_Zcomponent) / norm(findif_z);
            end

            tol = 1e-8;
            testCase.verifyLessThan(rel_errs_x(end), tol, "S'' fails finite difference check in the x-direction for the exterior case.");
            testCase.verifyLessThan(rel_errs_y(end), tol, "S'' fails finite difference check in the y-direction for the exterior case.");
            testCase.verifyLessThan(rel_errs_z(end), tol, "S'' fails finite difference check in the z-direction for the exterior case.");

            %%%
            %%% INTERIOR CHECK
            %%%

            % Setup target points away from surface
            nu_trg = get_norm_vecs(p, params.u0, params.oblate);
            X_trg = prolate_spheroid_shape(p, params.u0, params.a) - 1*nu_trg;

            % Evaluate spheroidal grad div SLP
            [graddivSL_Ucomponent, graddivSL_Vcomponent, graddivSL_PHIcomponent] = spheroidalgraddivSL(params, sigma_x, sigma_y, sigma_z, X_trg);
            
            S_trg = cart2spheroidal(X_trg, params.a, params.oblate);
            u = S_trg(:,1); v = S_trg(:,2); phi = S_trg(:,3);
            [graddivSL_Xcomponent, graddivSL_Ycomponent, graddivSL_Zcomponent] = convert_from_prolate_basis(graddivSL_Ucomponent, graddivSL_Vcomponent, graddivSL_PHIcomponent, ...
                                                                                                                u, v, phi);
            %%% Calculate divergence of Laplace SLP
            rel_errs_x = zeros(1, n_eps);
            rel_errs_y = zeros(1, n_eps);
            rel_errs_z = zeros(1, n_eps);
            for i=1:n_eps
                epsilon = eps(i);
            
                % Calculate Cartesian gradient approximation
                findif_x = (divSLP(X_trg + epsilon*nu_x) - divSLP(X_trg - epsilon*nu_x))./(2*epsilon);
                findif_y = (divSLP(X_trg + epsilon*nu_y) - divSLP(X_trg - epsilon*nu_y))./(2*epsilon);
                findif_z = (divSLP(X_trg + epsilon*nu_z) - divSLP(X_trg - epsilon*nu_z))./(2*epsilon);

                rel_errs_x(i) = norm(findif_x - graddivSL_Xcomponent) / norm(findif_x);
                rel_errs_y(i) = norm(findif_y - graddivSL_Ycomponent) / norm(findif_y);
                rel_errs_z(i) = norm(findif_z - graddivSL_Zcomponent) / norm(findif_z);
            end

            tol = 1e-8;
            testCase.verifyLessThan(rel_errs_x(end), tol, "S'' fails finite difference check in the x-direction for the interior case.");
            testCase.verifyLessThan(rel_errs_y(end), tol, "S'' fails finite difference check in the y-direction for the interior case.");
            testCase.verifyLessThan(rel_errs_z(end), tol, "S'' fails finite difference check in the z-direction for the interior case.");

        end

        function testProlateGradientCheckOffSurface(testCase)
            eps = 10.^(-2:-1:-5);
            n_eps = numel(eps);
            p = 8;
            np = 2*p*(p+1);
            params = SpheroidalParameters;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = 0;
            params.centers = [0 0 0];
            params.isReal = true;

            sigma_x = rand(np,1)+0.5;
            sigma_y = rand(np,1)-0.5; sigma_z = rand(np,1);
            % sigma_y = zeros(np,1); sigma_z = sigma_y;
            densities = { sigma_x, sigma_y, sigma_z };

            % Setup target points away from surface
            nu_trg = get_norm_vecs(p, params.u0, params.oblate);
            X_trg = prolate_spheroid_shape(p, params.u0, params.a) + 2*nu_trg;

            % Evaluate spheroidal grad div SLP
            [graddivSL_Ucomponent, graddivSL_Vcomponent, graddivSL_PHIcomponent] = spheroidalgraddivSL(params, sigma_x, sigma_y, sigma_z, X_trg);
            
            %%% Calculate divergence of Laplace SLP
            nu_x = repmat([1,0,0],size(X_trg,1),1);
            nu_y = repmat([0,1,0],size(X_trg,1),1);
            nu_z = repmat([0,0,1],size(X_trg,1),1);
            nu = { nu_x, nu_y, nu_z };
            function sum = divSLP(X)
                sum = 0;
                for j = 1:3
                    params.sigma = densities{j}; params.get_shc;
                    sum = sum + spheroidalSP(params, X, nu{j});
                end
            end

            %%% Get spheroidal coordinates and scale factor for the target
            %%% points
            S_trg = cart2spheroidal(X_trg, params.a, params.oblate);
            u = S_trg(:,1); v = S_trg(:,2); phi = S_trg(:,3);
            a = params.a;
            h_u = a * sqrt((u.^2 - v.^2)./(u.^2 - 1));
            h_v = a * sqrt((u.^2 - v.^2)./(1 - v.^2));
            h_phi = a * sqrt(u.^2 - 1) .* sqrt(1 - v.^2);

            %%% u direction
            u_rel_errs = zeros(1, n_eps);
            for i=1:n_eps
                epsilon = eps(i);
            
                X_plus_cart  = spheroidal2cart([u + epsilon, v, phi], params.a, params.oblate);
                X_minus_cart = spheroidal2cart([u - epsilon, v, phi], params.a, params.oblate);
                f_plus  = divSLP(X_plus_cart);
                f_minus = divSLP(X_minus_cart);
            
                df_du_approx = (f_plus - f_minus) / (2*epsilon);

                u_rel_errs(i) = norm(df_du_approx./h_u - graddivSL_Ucomponent) / norm(graddivSL_Ucomponent);
            end

            %%% v direction
            v_rel_errs = zeros(1, n_eps);
            for i=1:n_eps
                epsilon = eps(i);
            
                % Define perturbed points in spheroidal space and convert to Cartesian
                X_plus_cart  = spheroidal2cart([u, v + epsilon, phi], params.a, params.oblate);
                X_minus_cart = spheroidal2cart([u, v - epsilon, phi], params.a, params.oblate);
            
                % Calculate divergence at perturbed points
                f_plus  = divSLP(X_plus_cart);
                f_minus = divSLP(X_minus_cart);
            
                df_dv_approx = (f_plus - f_minus) / (2*epsilon);

                v_rel_errs(i) = norm(df_dv_approx./h_v - graddivSL_Vcomponent) / norm(graddivSL_Vcomponent);
            end

            %%% phi direction
            phi_rel_errs = zeros(1, n_eps);
            for i=1:n_eps
                epsilon = eps(i);
            
                X_plus_cart  = spheroidal2cart([u, v, phi + epsilon], params.a, params.oblate);
                X_minus_cart = spheroidal2cart([u, v, phi - epsilon], params.a, params.oblate);
                f_plus  = divSLP(X_plus_cart);
                f_minus = divSLP(X_minus_cart);

                df_dphi_approx = (f_plus - f_minus) / (2*epsilon);

                phi_rel_errs(i) = norm(df_dphi_approx./h_phi - graddivSL_PHIcomponent) / norm(graddivSL_PHIcomponent);
            end

            %%% Comparison
            testCase.verifyLessThan(u_rel_errs(end), 1e-8, "S'' fails finite difference check in the u-direction.");
            testCase.verifyLessThan(v_rel_errs(end), 1e-8, "S'' fails finite difference check in the u-direction.");
            testCase.verifyLessThan(phi_rel_errs(end), 1e-8, "S'' fails finite difference check in the u-direction.");
        end

        function testCartesianOblateGradientCheck(testCase)
            %{
                Converts the resulting vector we get from
                spheroidalgraddivSL to spheroidal coordinates and then use
                a finite difference check. This is faster than the
                spheroidal implementation of the test.
            %}
            eps = 10.^(-2:-1:-4);
            n_eps = numel(eps);
            p = 16;
            np = 2*p*(p+1);
            params = SpheroidalParameters;
            params.u0 = 8/sqrt(2);
            params.a = 1/sqrt(1 + params.u0);
            params.oblate = true;
            params.centers = [0 0 0];
            params.isReal = true;

            rng(42);
            sigma_x = rand(np,1)+0.5;
            sigma_y = rand(np,1)-0.5;
            sigma_z = rand(np,1)+1.0;
            % sigma_y = zeros(np,1); sigma_z = sigma_y;
            densities = { sigma_x, sigma_y, sigma_z };

            %%%
            %%% EXTERIOR CHECK
            %%%
            % Setup target points away from surface
            nu_trg = get_norm_vecs(p, params.u0, params.oblate);
            X_trg = oblate_spheroid_shape(p, params.u0, params.a) + 3*nu_trg;

            % Evaluate spheroidal grad div SLP
            [graddivSL_Ucomponent, graddivSL_Vcomponent, graddivSL_PHIcomponent] = spheroidalgraddivSL(params, sigma_x, sigma_y, sigma_z, X_trg);
            
            S_trg = cart2spheroidal(X_trg, params.a, params.oblate);
            u = S_trg(:,1); v = S_trg(:,2); phi = S_trg(:,3);
            [graddivSL_Xcomponent, graddivSL_Ycomponent, graddivSL_Zcomponent] = convert_from_oblate_basis(graddivSL_Ucomponent, graddivSL_Vcomponent, graddivSL_PHIcomponent, ...
                                                                                                                u, v, phi);
            %%% Calculate divergence of Laplace SLP
            nu_x = repmat([1,0,0],size(X_trg,1),1);
            nu_y = repmat([0,1,0],size(X_trg,1),1);
            nu_z = repmat([0,0,1],size(X_trg,1),1);
            nu = { nu_x, nu_y, nu_z };
            function sum = divSLP(X)
                sum = 0;
                for j = 1:3
                    params.sigma = densities{j}; params.get_shc;
                    sum = sum + spheroidalSP(params, X, nu{j});
                end
            end

            rel_errs_x = zeros(1, n_eps);
            rel_errs_y = zeros(1, n_eps);
            rel_errs_z = zeros(1, n_eps);
            for i=1:n_eps
                epsilon = eps(i);
            
                % Calculate Cartesian gradient approximation
                findif_x = (divSLP(X_trg + epsilon*nu_x) - divSLP(X_trg - epsilon*nu_x))./(2*epsilon);
                findif_y = (divSLP(X_trg + epsilon*nu_y) - divSLP(X_trg - epsilon*nu_y))./(2*epsilon);
                findif_z = (divSLP(X_trg + epsilon*nu_z) - divSLP(X_trg - epsilon*nu_z))./(2*epsilon);

                rel_errs_x(i) = norm(findif_x - graddivSL_Xcomponent) / norm(findif_x);
                rel_errs_y(i) = norm(findif_y - graddivSL_Ycomponent) / norm(findif_y);
                rel_errs_z(i) = norm(findif_z - graddivSL_Zcomponent) / norm(findif_z);
            end

            tol = 9e-6;
            testCase.verifyLessThan(rel_errs_x(end), tol, "S'' fails finite difference check in the x-direction for the exterior case.");
            testCase.verifyLessThan(rel_errs_y(end), tol, "S'' fails finite difference check in the y-direction for the exterior case.");
            testCase.verifyLessThan(rel_errs_z(end), tol, "S'' fails finite difference check in the z-direction for the exterior case.");
        end

        function testCartesianOblateGradientCheckInterior(testCase)
            eps = 10.^(-4:-1:-5);
            n_eps = numel(eps);
            p = 16;
            np = 2*p*(p+1);
            params = SpheroidalParameters;
            params.u0 = 2/sqrt(3);
            params.a = 1/sqrt(1 + params.u0);
            params.oblate = true;
            params.centers = [0 0 0];
            params.isReal = true;

            rng(42);
            sigma_x = rand(np,1)+0.5;
            sigma_y = rand(np,1)-0.5;
            sigma_z = rand(np,1)+1.0;
            densities = { sigma_x, sigma_y, sigma_z };

            function sum = divSLP(X)
                sum = 0;
                for j = 1:3
                    params.sigma = densities{j}; params.get_shc;
                    sum = sum + spheroidalSP(params, X, nu{j});
                end
            end

            % Setup target points away from surface
            nu_trg = get_norm_vecs(p, params.u0, params.oblate);
            X_trg = oblate_spheroid_shape(p, params.u0, params.a) + 3*nu_trg;

            % Evaluate spheroidal grad div SLP
            [graddivSL_Ucomponent, graddivSL_Vcomponent, graddivSL_PHIcomponent] = spheroidalgraddivSL(params, sigma_x, sigma_y, sigma_z, X_trg);
            
            S_trg = cart2spheroidal(X_trg, params.a, params.oblate);
            u = S_trg(:,1); v = S_trg(:,2); phi = S_trg(:,3);
            [graddivSL_Xcomponent, graddivSL_Ycomponent, graddivSL_Zcomponent] = convert_from_oblate_basis(graddivSL_Ucomponent, graddivSL_Vcomponent, graddivSL_PHIcomponent, ...
                                                                                                                u, v, phi);
            %%% Calculate divergence of Laplace SLP
            nu_x = repmat([1,0,0],size(X_trg,1),1);
            nu_y = repmat([0,1,0],size(X_trg,1),1);
            nu_z = repmat([0,0,1],size(X_trg,1),1);
            nu = { nu_x, nu_y, nu_z };

            rel_errs_x = zeros(1, n_eps);
            rel_errs_y = zeros(1, n_eps);
            rel_errs_z = zeros(1, n_eps);
            for i=1:n_eps
                epsilon = eps(i);
            
                % Calculate Cartesian gradient approximation
                findif_x = (divSLP(X_trg + epsilon*nu_x) - divSLP(X_trg - epsilon*nu_x))./(2*epsilon);
                findif_y = (divSLP(X_trg + epsilon*nu_y) - divSLP(X_trg - epsilon*nu_y))./(2*epsilon);
                findif_z = (divSLP(X_trg + epsilon*nu_z) - divSLP(X_trg - epsilon*nu_z))./(2*epsilon);

                rel_errs_x(i) = norm(findif_x - graddivSL_Xcomponent) / norm(findif_x);
                rel_errs_y(i) = norm(findif_y - graddivSL_Ycomponent) / norm(findif_y);
                rel_errs_z(i) = norm(findif_z - graddivSL_Zcomponent) / norm(findif_z);
            end

            tol = 9e-6;
            testCase.verifyLessThan(rel_errs_x(end), tol, "S'' fails finite difference check in the x-direction for the interior case.");
            testCase.verifyLessThan(rel_errs_y(end), tol, "S'' fails finite difference check in the y-direction for the interior case.");
            testCase.verifyLessThan(rel_errs_z(end), tol, "S'' fails finite difference check in the z-direction for the interior case.");
        end

        function testOblateGradientCheckOffSurface(testCase)
            eps = 10.^(-3:-1:-6);
            n_eps = numel(eps);
            p = 16;
            np = 2*p*(p+1);
            params = SpheroidalParameters;
            params.u0 = 2/sqrt(3);
            params.a = testCase.a_oblate;
            params.oblate = true;
            params.centers = [0 0 0];
            params.isReal = true;

            rng(42);
            sigma_x = rand(np,1)+0.5;
            sigma_y = rand(np,1)-0.5; sigma_z = rand(np,1);
            % sigma_y = zeros(np,1); sigma_z = zeros(np,1);
            densities = { sigma_x, sigma_y, sigma_z };

            % Setup target points away from surface
            nu_trg = get_norm_vecs(p, params.u0, params.oblate);
            X_trg = oblate_spheroid_shape(p, params.u0, params.a) + 2*nu_trg;

            % Evaluate spheroidal grad div SLP
            [graddivSL_Ucomponent, graddivSL_Vcomponent, graddivSL_PHIcomponent] = spheroidalgraddivSL(params, sigma_x, sigma_y, sigma_z, X_trg);
            
            %%% Calculate divergence of Laplace SLP
            nu_x = repmat([1,0,0],size(X_trg,1),1);
            nu_y = repmat([0,1,0],size(X_trg,1),1);
            nu_z = repmat([0,0,1],size(X_trg,1),1);
            nu = { nu_x, nu_y, nu_z };
            function sum = divSLP(X)
                sum = 0;
                for j = 1:3
                    params.sigma = densities{j}; params.get_shc;
                    sum = sum + spheroidalSP(params, X, nu{j});
                end
            end

            %%% Get spheroidal coordinates and scale factor for the target
            %%% points
            S_trg = cart2spheroidal(X_trg, params.a, params.oblate);
            u = S_trg(:,1); v = S_trg(:,2); phi = S_trg(:,3);
            a = params.a;
            h_u = a * sqrt((u.^2 + v.^2)./(u.^2 + 1));
            h_v = a .* sqrt((u.^2 + v.^2)./(1 - v.^2));
            h_phi = a .* sqrt(u.^2 + 1) .* sqrt(1 - v.^2);

            %%% u direction
            u_rel_errs = zeros(1, n_eps);
            for i=1:n_eps
                epsilon = eps(i);

                X_plus_cart  = spheroidal2cart([u + epsilon, v, phi], params.a, params.oblate);
                X_minus_cart = spheroidal2cart([u - epsilon, v, phi], params.a, params.oblate);
                f_plus  = divSLP(X_plus_cart);
                f_minus = divSLP(X_minus_cart);

                df_du_approx = (f_plus - f_minus) / (2*epsilon);

                u_rel_errs(i) = norm(df_du_approx./h_u - graddivSL_Ucomponent) / norm(graddivSL_Ucomponent);
            end

            %%% v direction
            v_rel_errs = zeros(1, n_eps);
            for i=1:n_eps
                epsilon = eps(i);
            
                % Define perturbed points in spheroidal space and convert to Cartesian
                X_plus_cart  = real(spheroidal2cart([u, v + epsilon, phi], params.a, params.oblate));
                X_minus_cart = real(spheroidal2cart([u, v - epsilon, phi], params.a, params.oblate));
            
                % Calculate divergence at perturbed points
                f_plus  = divSLP(X_plus_cart);
                f_minus = divSLP(X_minus_cart);
            
                df_dv_approx = (f_plus - f_minus) / (2*epsilon);

                v_rel_errs(i) = norm(df_dv_approx./h_v - graddivSL_Vcomponent) / norm(graddivSL_Vcomponent);
            end

            %%% phi direction
            phi_rel_errs = zeros(1, n_eps);
            for i=1:n_eps
                epsilon = eps(i);

                X_plus_cart  = spheroidal2cart([u, v, phi + epsilon], params.a, params.oblate);
                X_minus_cart = spheroidal2cart([u, v, phi - epsilon], params.a, params.oblate);
                f_plus  = divSLP(X_plus_cart);
                f_minus = divSLP(X_minus_cart);

                df_dphi_approx = (f_plus - f_minus) / (2*epsilon);

                phi_rel_errs(i) = norm(df_dphi_approx./h_phi - graddivSL_PHIcomponent) / norm(graddivSL_PHIcomponent);
            end

            %%% Comparison
            testCase.verifyLessThan(u_rel_errs(end), 1e-8, "S'' fails finite difference check in the u-direction.");
            testCase.verifyLessThan(v_rel_errs(end), 1e-8, "S'' fails finite difference check in the v-direction.");
            testCase.verifyLessThan(phi_rel_errs(end), 1e-8, "S'' fails finite difference check in the phi-direction.");
        end

        %%%
        %%% Convergence check
        %%%
        function testOblateExteriorConvergenceCheck(testCase)
        end

        function testOblateInteriorConvergenceCheck(testCase)
        end

        %%%
        %%% Formula check in terms of single derivative
        %%%
        function testProlateFormulaCheckInTermsOfFirstDerivative(testCase)
            %{
                Checks if the formula derived for S'' even works.
            %}
            rng(42);
            p = 16;
            np = 2*p*(p+1);
            params = SpheroidalParameters;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = false;
            params.centers = [0 0 0];
            params.isReal = true;

            [u, v] = gl_grid(p);
            sigma_x = cos(u) .* sin(u).^2 + sin(u) .* cos(v);
            sigma_y = cos(u).^2 + sin(u).*cos(u).*sin(v);
            sigma_z = cos(u) .* sin(u).^2 - sin(u) .* cos(v) + exp(cos(u).^3);
            densities = { sigma_x, sigma_y, sigma_z };

            X_self = prolate_spheroid_shape(p, params.u0, params.a);
            nu_trg = get_norm_vecs(p, params.u0, params.oblate);
            X_trg = X_self + nu_trg;

            %%% Evaluate spheroidal grad div SLP
            [graddivSL_Ucomponent, graddivSL_Vcomponent, graddivSL_PHIcomponent] = spheroidalgraddivSL(params, sigma_x, sigma_y, sigma_z, X_trg);
            S_trg = cart2spheroidal(X_trg, params.a, params.oblate);
            u = S_trg(:,1); v = S_trg(:,2); phi = S_trg(:,3);
            [graddivSL_X, graddivSL_Y, graddivSL_Z] = convert_from_prolate_basis(graddivSL_Ucomponent, graddivSL_Vcomponent, graddivSL_PHIcomponent, u, v, phi);
            graddivSL_result = graddivSL_X .* nu_trg(:,1) + graddivSL_Y .* nu_trg(:,2) + graddivSL_Z .* nu_trg(:,3);

            %%% Testing
            nu_x = repmat([1,0,0],size(X_trg,1),1);
            nu_y = repmat([0,1,0],size(X_trg,1),1);
            nu_z = repmat([0,0,1],size(X_trg,1),1);
            nu_x_sph = cartNu2spheroidal(nu_x, S_trg, params.a, params.oblate);
            nu_y_sph = cartNu2spheroidal(nu_y, S_trg, params.a, params.oblate);
            nu_z_sph = cartNu2spheroidal(nu_z, S_trg, params.a, params.oblate);
            sigma_x = sigma_x .* X_self(:,1);
            sigma_y = sigma_y .* X_self(:,1);
            sigma_z = sigma_z .* X_self(:,1);
            [graddivSL_Ucomponent, graddivSL_Vcomponent, graddivSL_PHIcomponent] = spheroidalgraddivSL(params, sigma_x, sigma_y, sigma_z, X_trg);
            graddivSL_X = graddivSL_Ucomponent .* nu_x_sph(:,1) + graddivSL_Vcomponent .* nu_x_sph(:,2) + graddivSL_PHIcomponent .* nu_x_sph(:,3);
            graddivSL_Y = graddivSL_Ucomponent .* nu_y_sph(:,1) + graddivSL_Vcomponent .* nu_y_sph(:,2) + graddivSL_PHIcomponent .* nu_y_sph(:,3);
            graddivSL_Z = graddivSL_Ucomponent .* nu_z_sph(:,1) + graddivSL_Vcomponent .* nu_z_sph(:,2) + graddivSL_PHIcomponent .* nu_z_sph(:,3);

            %%% Formula calculation
            % Calculate the surface divergence of the density.
            X_self = prolate_spheroid_shape(p, params.u0, params.a);
            S_self = cart2spheroidal(X_self, params.a, params.oblate);
            surf_div_density = calculate_surf_div(params, S_self, sigma_x, sigma_y, sigma_z); % divergence of y1 times sig

            %%% Calculate SP term
            params.sigma = surf_div_density; params.get_shc;
            % SP_term = spheroidalSP(params, X_trg, nu_trg);
            [SP_term_x, SP_term_y, SP_term_z] = spheroidalSP(params, X_trg, nu_x, nu_y, nu_z);

            %%% Calculate DP term
            % Here, the normals on the surface are equal to the target
            % normals.
            n_dot_sigma = nu_trg(:,1).*sigma_x + nu_trg(:,2).*sigma_y + nu_trg(:,3).*sigma_z; %n(y) \cdot (y1 times sigma)
            params.sigma = n_dot_sigma; params.get_shc;
            % DP_term = spheroidalDP(params, X_trg, nu_trg);
            [DP_term_x, DP_term_y, DP_term_z] = spheroidalDP(params, X_trg, nu_x, nu_y, nu_z);

            %%% Compare
            tol = 1e-6;
            rel_errs = [
                norm(SP_term_x - DP_term_x - graddivSL_X) / norm(graddivSL_X);
                norm(SP_term_y - DP_term_y - graddivSL_Y) / norm(graddivSL_Y);
                norm(SP_term_z - DP_term_z - graddivSL_Z) / norm(graddivSL_Z);
            ];
            testCase.verifyLessThan(rel_errs, tol, "Formula not satisfied for prolate case.");
        end

        function testOblateFormulaCheckInTermsOfFirstDerivative(testCase)
            %{
                Checks if the formula derived for S'' even works.
            %}
            rng(42);
            p = 16;
            np = 2*p*(p+1);
            params = SpheroidalParameters;
            params.u0 = testCase.u0_oblate;
            params.a = testCase.a_oblate;
            params.oblate = true;
            params.centers = [0 0 0];
            params.isReal = true;

            [u, v] = gl_grid(p);
            sigma_x = cos(u) .* sin(u).^2 + sin(u) .* cos(v);
            sigma_y = cos(u).^2 + sin(u).*cos(u).*sin(v);
            sigma_z = cos(u) .* sin(u).^2 - sin(u) .* cos(v) + exp(cos(u).^3);
            densities = { sigma_x, sigma_y, sigma_z };

            nu_trg = get_norm_vecs(p, params.u0, params.oblate);
            X_trg = prolate_spheroid_shape(p, params.u0, params.a) + nu_trg;

            %%% Evaluate spheroidal grad div SLP
            [graddivSL_Ucomponent, graddivSL_Vcomponent, graddivSL_PHIcomponent] = spheroidalgraddivSL(params, sigma_x, sigma_y, sigma_z, X_trg);
            S_trg = cart2spheroidal(X_trg, params.a, params.oblate);
            u = S_trg(:,1); v = S_trg(:,2); phi = S_trg(:,3);
            [graddivSL_X, graddivSL_Y, graddivSL_Z] = convert_from_oblate_basis(graddivSL_Ucomponent, graddivSL_Vcomponent, graddivSL_PHIcomponent, u, v, phi);
            graddivSL_result = graddivSL_X .* nu_trg(:,1) + graddivSL_Y .* nu_trg(:,2) + graddivSL_Z .* nu_trg(:,3);

            %%% Formula calculation
            % Calculate the surface divergence of the density.
            X_self = oblate_spheroid_shape(p, params.u0, params.a);
            S_self = cart2spheroidal(X_self, params.a, params.oblate);
            surf_div_density = calculate_surf_div(params, S_self, sigma_x, sigma_y, sigma_z);

            %%% Calculate SP term
            params.sigma = surf_div_density; params.get_shc;
            SP_term = spheroidalSP(params, X_trg, nu_trg);

            %%% Calculate DP term
            % Here, the normals on the surface are equal to the target
            % normals.
            n_dot_sigma = nu_trg(:,1).*sigma_x + nu_trg(:,2).*sigma_y + nu_trg(:,3).*sigma_z;
            params.sigma = n_dot_sigma; params.get_shc;
            DP_term = spheroidalDP(params, X_trg, nu_trg);

            %%% Compare
            tol = 1e-6;
            testCase.verifyLessThan(norm(SP_term - DP_term - graddivSL_result), tol, "Formula not satisfied for the oblate case.");
        end

        %%%
        %%% On-surface tests
        %%%
        function testOnSurfaceCheckWithFormula(testCase)
            %{
                Checks if the formula for S'' works on-surface.
            %}
            %%%
            %%% Setup for prolate case
            %%%
            rng(23);
            tol = 1e-6;
            p = 16;
            np = 2*p*(p+1);
            params = SpheroidalParameters;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = false;
            params.centers = [0 0 0];
            params.isReal = true;

            [u, v] = gl_grid(p);
            sigma_x = Ynm(14, 13, u, v);
            sigma_y = Ynm(15, 10, u, v);
            sigma_z = Ynm(14, 13, u, v);

            nu_trg = get_norm_vecs(p, params.u0, params.oblate);
            X_self = prolate_spheroid_shape(p, params.u0, params.a);
            S_self = cart2spheroidal(X_self, params.a, params.oblate);
            u = S_self(:,1); v = S_self(:,2); phi = S_self(:,3);

            %%% Evaluate spheroidal grad div SLP
            [graddivSL_Ucomponent, graddivSL_Vcomponent, graddivSL_PHIcomponent] = spheroidalgraddivSL(params, sigma_x, sigma_y, sigma_z, X_self);
            graddivSL_result = graddivSL_Ucomponent; % Normal derivative with respect to e_u

            %%% Formula calculation
            % Calculate the surface divergence of the density.
            surf_div_density = calculate_surf_div(params, S_self, sigma_x, sigma_y, sigma_z);

            %%% Calculate SP term
            params.sigma = surf_div_density; params.get_shc;
            SP_term = spheroidalSP(params, [], nu_trg);

            %%% Calculate DP term
            % Here, the normals on the surface are equal to the target
            % normals.
            n_dot_sigma = nu_trg(:,1).*sigma_x + nu_trg(:,2).*sigma_y + nu_trg(:,3).*sigma_z;
            params.sigma = n_dot_sigma; params.get_shc;
            DP_term = spheroidalDP(params, [], nu_trg);

            %%% Compare
            testCase.verifyLessThan(norm(SP_term - DP_term - graddivSL_result)./norm(graddivSL_result), tol, "Formula not satisfied for prolate case.");
        end

        function testSurfaceIsAverage(testCase)
            p = 16;
            np = 2*p*(p+1);
            params = SpheroidalParameters;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = false;
            params.centers = [0 0 0];
            params.isReal = true;

            [u, v] = gl_grid(p);
            sigma_x = sin(0.005*[1:1:np]).';
            sigma_y = cos(0.005*[1:1:np]).';
            sigma_z = sin(0.005*[1:1:np]).';
            
            X_surf = prolate_spheroid_shape(p, params.u0, params.a);
            S_self = cart2spheroidal(X_surf, params.a, params.oblate);
            nu_trg = get_norm_vecs(p, params.u0, params.oblate);

            eps = 1e-5;
            X_int = X_surf - eps * nu_trg;
            X_ext = X_surf + eps * nu_trg;

            surf_div_density = calculate_surf_div(params, S_self, sigma_x, sigma_y, sigma_z);
            
            [graddivSL_U_surf, ~, ~] = spheroidalgraddivSL(params, sigma_x, sigma_y, sigma_z, X_surf);
            [graddivSL_U_int, ~, ~] = spheroidalgraddivSL(params, sigma_x, sigma_y, sigma_z, X_int);
            [graddivSL_U_ext, ~, ~] = spheroidalgraddivSL(params, sigma_x, sigma_y, sigma_z, X_ext);
            
            average_val = (graddivSL_U_int + graddivSL_U_ext) / 2;

            tol = 1e-6;
            testCase.verifyLessThan(norm(graddivSL_U_surf - average_val) / norm(average_val), tol, ...
                'On-surface evaluation is not the average of the interior and exterior limits.');

            %%% Jump relation
            testCase.verifyLessThan(norm((graddivSL_U_ext - graddivSL_U_int) - -1*surf_div_density)./norm(surf_div_density), 1e-4, ...
                'Jump relation not satisfied according to the formula.');
        end

        function testJumpRelationConvergenceTest(testCase)
            p = 16;
            np = 2*p*(p+1);
            params = SpheroidalParameters;
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = false;
            params.centers = [0 0 0];
            params.isReal = true;

            [u, v] = gl_grid(p);
            sigma_x = 5*sin(0.005*u);
            sigma_y = sin(v).^3 .* cos(u);
            sigma_z = sin(u).^2 .* cos(2*v);

            X_surf = prolate_spheroid_shape(p, params.u0, params.a);
            S_self = cart2spheroidal(X_surf, params.a, params.oblate);
          
            nu_trg = get_norm_vecs(p, params.u0, params.oblate);
            surf_div_density = calculate_surf_div(params, S_self, p, sigma_x, sigma_y, sigma_z);

            eps = 1e-9;
            X_int = X_surf - eps * nu_trg;
            X_ext = X_surf + eps * nu_trg;

            params.sigma = surf_div_density; params.get_shc;
            SL_int = spheroidalSP(params, X_int, nu_trg);
            SL_surf = spheroidalSP(params, X_surf, nu_trg);
            SL_ext = spheroidalSP(params, X_ext, nu_trg);

            %%% Simultaneously test the jump relation.
            calculated_jump = SL_ext - SL_int;
            analytical_jump = -1*surf_div_density;
            testCase.verifyLessThan(norm(calculated_jump - analytical_jump)/norm(analytical_jump), 1e-6, "Jump condition failed.")
        end

        %%%
        %%% Misc. tests
        %%%
        function testRotationMatrix(testCase)
            eps = 10.^(-3:-1:-6);
            n_eps = numel(eps);
            p = 8;
            np = 2*p*(p+1);
            params = SpheroidalParameters;
            params.u0 = 60;
            params.a = testCase.a_oblate;
            params.oblate = true;
            params.centers = [0 0 0];
            params.isReal = true;

            % Setup target points away from surface
            nu_trg = get_norm_vecs(p, params.u0, params.oblate);
            X_trg = oblate_spheroid_shape(p, params.u0, params.a) + 3*nu_trg;

            pt_idx = 1; % Pick one point to test
            X_pt = X_trg(pt_idx, :);
            S_pt = cart2spheroidal(X_pt, params.a, params.oblate);
            u = S_pt(1); v = S_pt(2); phi = S_pt(3);
            a = params.a;
            
            h_u = a * sqrt((u^2 + v^2)/(u^2 + 1));
            h_v = a * sqrt((u^2 + v^2)/(1 - v^2));
            h_phi = a * sqrt(u^2 + 1) * sqrt(1 - v^2);
            
            eps_fd = 1e-7; % Small epsilon for finite difference
        end

        function testProlateSpectralCoefficientsOnSurface(testCase)
            %{
                Uses the orthogonality of the spheroidal harmonics to extract
                the coefficients of the Laplace double-layer potential.
            %}
            geti = @(n,m) m+n^2+n+1; % Map (n,m) to 0 <= k <= sp
            p = 16;
            np = 2*p*(p+1);
            sp = (p+1)^2;
            
            params = SpheroidalParameters;
            params.p = p;
            params.isReal = false; % Must be false!
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = false;

            % Create a matrix where each column is a spherical harmonic 
            % basis function evaluated at the grid points.
            [u,v] = gl_grid(p);
            Y = zeros(np, sp);
            ii = (1:sp)'; 
            nn = floor(sqrt(ii-1)); 
            mm = ii - nn.^2 - nn - 1;
            for k = 1:sp
                n = nn(k); m = mm(k);
                Y(:,k) = Ynm(n,m,u,v);
            end

            params.sigma = Y; params.get_shc;
            shc = params.sigma_coefficients;
            u0 = params.u0;
            oblate = params.oblate;
            all_u0s=unique(u0);
            Gshc = zeros(size(shc));
            for i=1:length(all_u0s)
                this_u0=all_u0s(i);
            
                %Calculate basis transformation matrix G: Ynm -> Ynm / sqrt(u0^2-v^2)
                %or /sqrt(u0^2+v^2) for oblate
                this_G_pro = Gmatrix(p,this_u0,0,0);
                this_G_obl = Gmatrix(p,this_u0,0,1);
            
                % keep track of indices and shc of particles with the same u0
                this_index = find(u0==this_u0); 
                obl_index = nonzeros(this_index.*oblate(this_index));
                pro_index = nonzeros(this_index.*~oblate(this_index));
                obl_shc = reshape(shc(:,:,obl_index),sp,[],1);
                pro_shc = reshape(shc(:,:,pro_index),sp,[],1);
            
                %Apply basis transformation to spheroidal harmonic coefficients
                obl_Gshc = this_G_obl\obl_shc;
                pro_Gshc = this_G_pro\pro_shc;
                
                Gshc(:,:,obl_index) = reshape(obl_Gshc,sp,[],length(obl_index));
                Gshc(:,:,pro_index) = reshape(pro_Gshc,sp,[],length(pro_index));
            end

            G = Gmatrix(p, params.u0, 0, 0);
            params.sigma = Y;
            SL_Ynm = spheroidalSL(params);
            SL_coeffs = shAna(SL_Ynm);

            ii = (1:sp)';
            nn = floor(sqrt(ii-1));
            mm = ii - nn.^2 - nn - 1;
            bnm = params.a * factorial(nn-mm)./factorial(nn+mm) .* ((-1) .^ (mm)) .* sqrt(params.u0.^2 - 1);
            L = legendre_otc(p, params.u0, 1, 1, 1);
            P = L{1}; Q = L{2};
            legendre_terms = P .* Q;
            coefficient = bnm .* legendre_terms;
            nf = size(params.sigma_coefficients, 2);
            coefficient_mtx = repmat(coefficient,1,nf,1);
            expected_coeffs = coefficient_mtx .* Gshc;

            errors = size(sp, 1);

            n=1; m=2;
            expected_SL = Y * expected_coeffs;
            expected_coefficient = shAna(expected_SL(:,geti(n,m)));
            norm(expected_coefficient - SL_coeffs(:,geti(n,m)))

            for k=1:sp
                n = nn(k); m = mm(k);
                expected_coefficient = shAna(expected_SL(:,geti(n,m)));
                errors(k) = norm(expected_coefficient - SL_coeffs(:,geti(n,m)));
            end

            testCase.verifyLessThan(errors, 9e-12, "On-surface spectral coefficients for the prolate case do not match the expected values to a reasonable tolernace");
        end

        function testProlateSpectralCoefficientsOffSurface(testCase)
            %{
                Uses the orthogonality of the spheroidal harmonics to extract
                the coefficients of the Laplace double-layer potential.
            %}
            geti = @(n,m) m+n^2+n+1; % Map (n,m) to 0 <= k <= sp
            p = 16;
            np = 2*p*(p+1);
            sp = (p+1)^2;
            
            params = SpheroidalParameters;
            params.p = p;
            params.isReal = false; % Must be false!
            params.u0 = testCase.u0_prolate;
            params.a = testCase.a_prolate;
            params.oblate = false;

            nu_trg = get_norm_vecs(p, params.u0, params.oblate);
            X_trg = prolate_spheroid_shape(p, params.u0, params.a) + 3*nu_trg;

            % Create a matrix where each column is a spherical harmonic 
            % basis function evaluated at the grid points.
            S = cart2spheroidal(X_trg, params.a, params.oblate);
            u_x = S(:,1); v_x = S(:,2); phi_x = S(:,3);
            Y = zeros(np, sp);
            ii = (1:sp)'; 
            nn = floor(sqrt(ii-1)); 
            mm = ii - nn.^2 - nn - 1;
            for k = 1:sp
                n = nn(k); m = mm(k);
                Y(:,k) = Ynm(n,m,real(acos(v_x)),phi_x);
            end

            params.sigma = Y; params.get_shc;
            shc = params.sigma_coefficients;
            u0 = params.u0;
            oblate = params.oblate;
            all_u0s=unique(u0);
            Gshc = zeros(size(shc));
            for i=1:length(all_u0s)
                this_u0=all_u0s(i);
            
                %Calculate basis transformation matrix G: Ynm -> Ynm / sqrt(u0^2-v^2)
                %or /sqrt(u0^2+v^2) for oblate
                this_G_pro = Gmatrix(p,this_u0,0,0);
                this_G_obl = Gmatrix(p,this_u0,0,1);
            
                % keep track of indices and shc of particles with the same u0
                this_index = find(u0==this_u0); 
                obl_index = nonzeros(this_index.*oblate(this_index));
                pro_index = nonzeros(this_index.*~oblate(this_index));
                obl_shc = reshape(shc(:,:,obl_index),sp,[],1);
                pro_shc = reshape(shc(:,:,pro_index),sp,[],1);
            
                %Apply basis transformation to spheroidal harmonic coefficients
                obl_Gshc = this_G_obl\obl_shc;
                pro_Gshc = this_G_pro\pro_shc;
                
                Gshc(:,:,obl_index) = reshape(obl_Gshc,sp,[],length(obl_index));
                Gshc(:,:,pro_index) = reshape(pro_Gshc,sp,[],length(pro_index));
            end

            params.sigma = Y;
            SL_Ynm = spheroidalSL(params, X_trg);
            SL_coeffs = shAna(SL_Ynm);

            S = cart2spheroidal(X_trg, params.a, params.oblate);
            u_x = S(:,1);

            ii = (1:sp)';
            nn = floor(sqrt(ii-1));
            mm = ii - nn.^2 - nn - 1;
            bnm = params.a * factorial(nn-mm)./factorial(nn+mm) .* ((-1) .^ (mm)) .* sqrt(params.u0.^2 - 1);
            L = legendre_otc(p, params.u0, 1, 1, 1);
            Pu0 = L{1};
            Qu = solid_harmonic(p, params.u0, u_x);
            legendre_terms = Pu0;
            coefficient = bnm .* legendre_terms;
            nf = size(params.sigma_coefficients, 2);
            coefficient_mtx = repmat(coefficient,1,nf,1);
            expected_coeffs = coefficient_mtx .* Gshc;

            errors = size(sp, 1);

            n=1; m=2;
            expected_SL = (Qu .* Y) * expected_coeffs;
            expected_coefficient = shAna(expected_SL(:,geti(n,m)));
            norm(expected_coefficient - SL_coeffs(:,geti(n,m)))

            for k=1:sp
                n = nn(k); m = mm(k);
                expected_coefficient = shAna(expected_SL(:,geti(n,m)));
                errors(k) = norm(expected_coefficient - SL_coeffs(:,geti(n,m)));
            end

            testCase.verifyLessThan(errors, 9e-12, "On-surface spectral coefficients for the prolate case do not match the expected values to a reasonable tolernace");
        end
    end
end

%% Helper functions
function [Fr, Fp, Fpp]=solid_harmonic_prime(p, u0, u_x, oblate)
    %{
        Solid spheroidal harmonics can be written as f_n^m(u)Y_n^m(v, phi).
        This function returns what Fr = f_n^m is (and its derivative as Fp), depending 
        on whether we are in the exterior or interior.
    %}
    Fr = ones(size(u_x,1),(p+1)^2);
    Fp = ones(size(u_x,1),(p+1)^2);
    Fpp = ones(size(u_x,1),(p+1)^2);

    if oblate
        u_x = 1j.*u_x;
    end

    if abs(u_x)-u0 < -1e-14 % Interior
        PQ=legendre_otc(p,u_x,1,2,2);
        P=PQ{1}; dP=PQ{3}; ddP = PQ{5}
        Fr=P.'; Fp=dP.'; Fpp=ddP.';
    elseif abs(u_x)-u0 > 1e-14 % Exterior
        PQ=legendre_otc(p,u_x,1,2,2);
        Q=PQ{2}; dQ=PQ{4}; ddQ=PQ{6};
        Fr=Q.'; Fp=dQ.'; Fpp=ddQ.';
    end
end
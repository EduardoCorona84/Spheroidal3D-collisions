%{
    Test code for Laplace to Stokes double layer potential with two prolate
    spheroids.

    The density, sigma, is randomized.
%}


clear;

CHECK_LAPLACE_FLAG = true;
CHECK_STOKES_FLAG = false;

params=SpheroidalParameters; 
params.matvec_eta=10; 
params.u0=1.1;
params.a=1/1.1;
params.oblate=0; 
params.centers=[0 0 0];
p=8; np=2*p*(p+1);

sigma_x = rand(np,1)+0.5;
sigma_y = rand(np,1)-0.5;
sigma_z = rand(np,1);

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

if CHECK_LAPLACE_FLAG
    sigma_list = [sigma_x, sigma_y, sigma_z]; % Size Nself x 3
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
    
            fprintf("gradient of Laplace DP sigma_%d dx%d vs matvec max abs error: %e\n", sigma_ind, j, ...
                    max(norm( dDP_sig_vec(1, :, sigma_ind, j)' - dDP_matvec_vec(:, 1, sigma_ind, j) )) );
        end
    end
    fprintf('------ Done... ------\n');
end

if CHECK_STOKES_FLAG
    % First, we calculate the Stokes double layer potential manually.
    % Precomputations
    r_inv_5 = r .^ (-5);
    RdotN = Rvec_x .* nx_src + Rvec_y .* ny_src + Rvec_z .* nz_src; % (x-y) \cdot n(y)
    RdotSigma = Rvec_x .* sigma_x + Rvec_y .* sigma_y + Rvec_z .* sigma_z; % (x-y) \cdot \sigma(y)

    integrand_DLx = RdotN .* Rvec_x .* RdotSigma .* r_inv_5;
    integrand_DLy = RdotN .* Rvec_y .* RdotSigma .* r_inv_5;
    integrand_DLz = RdotN .* Rvec_z .* RdotSigma .* r_inv_5;

    % Integrate
    SLx = -(3/(4*pi)) * integrate(integrand_SLx);
    SLy = -(3/(4*pi)) * integrate(integrand_SLy);
    SLz = -(3/(4*pi)) * integrate(integrand_SLz);

    % Now, actually calculate result from L2Stk and compare.
    target_pts = cell(1, 1);
    target_pts{1} = Xtrg;
    [L2Stkx, L2Stky, L2Stkz] = L2StkDLP(target_pts, params, sigma_x, sigma_y, sigma_z, 1);

    fprintf("\n inf error of Stokes double layer potential in x: %e, in y: %e, in z: %e\n", ...
        max(norm(L2Stkx{1}-SLx)), max(norm(L2Stky{1}-SLy)),max(norm(L2Stkz{1}-SLz)));
end
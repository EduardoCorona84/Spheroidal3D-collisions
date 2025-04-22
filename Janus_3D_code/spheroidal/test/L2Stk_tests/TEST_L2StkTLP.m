%{
    Test code for Laplace to Stokes double layer potential with two prolate
    spheroids.

    The density, sigma, is randomized.
%}


clear;

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
target_u0 = 1.2;
target_a = 1/1.2;
Xtrg = prolate_spheroid_shape(p, target_u0, target_a)+[3,3,1];
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
norm_vecs = get_norm_vecs(p, target_u0, oblate_flag=false);
nx_target = norm_vecs(:,1);
ny_target = norm_vecs(:,2);
nz_target = norm_vecs(:,3);

% First, we calculate the Stokes double layer potential manually.
% Precomputations
r_inv_5 = r .^ (-5);
RdotN = Rvec_x .* nx_target + Rvec_y .* ny_target + Rvec_z .* nz_target; % (x-y) \cdot n(y)
RdotSigma = Rvec_x .* sigma_x + Rvec_y .* sigma_y + Rvec_z .* sigma_z; % (x-y) \cdot \sigma(y)

integrand_TLPx = RdotN .* Rvec_x .* RdotSigma .* r_inv_5;
integrand_TLPy = RdotN .* Rvec_y .* RdotSigma .* r_inv_5;
integrand_TLPz = RdotN .* Rvec_z .* RdotSigma .* r_inv_5;

% Integrate
TLPx = -(3/(4*pi)) * integrate(integrand_TLPx);
TLPy = -(3/(4*pi)) * integrate(integrand_TLPy);
TLPz = -(3/(4*pi)) * integrate(integrand_TLPz);

% Now, actually calculate result from L2Stk and compare.
target_pts = cell(1, 1);
target_pts{1} = Xtrg;
[L2StkTLPx, L2StkTLPy, L2StkTLPz] = L2StkDLP(target_pts, params, sigma_x, sigma_y, sigma_z, 1);

fprintf("\n inf error of Stokes traction layer potential in x: %e, in y: %e, in z: %e\n", ...
    max(norm(L2StkTLPx{1}-TLPx)), max(norm(L2StkTLPy{1}-TLPy)),max(norm(L2StkTLPz{1}-TLPz)));
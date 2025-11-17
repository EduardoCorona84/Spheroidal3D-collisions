function eval = SpheroidalMS_L2Stk_MatVec_near(pot, u0, a, oblate, sigma_x, sigma_y, sigma_z, target_points, target_normals)
%{
A version of L2StkMatVec implemented for the mobility solver codeHaThe reason for this new file is to be handle BOTH spheres and spheroids.

Note that the target points are assumed to be already rotated to the local frame of
spheroid i, and that this function assumes only 1 source spheroid (for now?).

Inputs
    pot - (string) the following cases are supported: 
        'SL_Stk_3D'   single layer, Stokes
        'DL_Stk_3D'   double layer, Stokes
        'TSL_Stk_3D'  traction kernel of single layer, Stokes
    u0 - (double) controls the aspect ratio of spheroid
    a - (double) scale parameter
    sigma_x - (double array) surface density for the x-component
        size np × 1
    sigma_y - (double array) surface density for the y-component
    sigma_z - (double array) surface density for the z-component
    target_points - (double) target points to evaluate layer potential at
    target_normals - (double) normals at target points
        only necessary if doing TSL_Stk_3D.

Outputs
    eval - (double, 3*np × 1) Interleaved velocity (or traction) values at targets
        [vx1; vy1; vz1; vx2; vy2; vz2; ...].
%}

%% Initialize source particle in its local frame
params_i = SpheroidalParameters();
params_i.sigma = sigma_x; % Force p to update
params_i.u0 = u0;
params_i.a = a;
params_i.oblate = oblate;
params_i.centers = [0 0 0];
params_i.thetas = 0;
params_i.phis = 0;
params_i.Rmat = eye(3);
ns = 1;

if isempty(target_points)
    X_eval = {params_i.get_X};
    Nu_eval = {params_i.get_Norm};
else
    X_eval = {target_points};
    Nu_eval = {target_normals};
end

switch pot
    case 'SL_Stk_3D'
        [vx, vy, vz] = L2Stk(X_eval, params_i, sigma_x, sigma_y, sigma_z, ns);
    case 'DL_Stk_3D'
        [vx, vy, vz] = L2StkDLP(X_eval, params_i, sigma_x, sigma_y, sigma_z, ns);
    case 'TSL_Stk_3D'
        [vx, vy, vz] = L2StkTLP(X_eval, Nu_eval, params_i, sigma_x, sigma_y, sigma_z, ns, false);
    otherwise
        error('Incorrect potential passed in.');
end

if iscell(vx); vx = vx{1}; end
if iscell(vy); vy = vy{1}; end
if iscell(vz); vz = vz{1}; end

% Interleave result
eval = reshape([vx(:), vy(:), vz(:)].', [], 1);

end

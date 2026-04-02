function parbd = SpheroidalMS_set_params(opts)
% function parbd = SpheroidalMS_set_params(equ_radii,polar_radii,p,C,collision_eps,mdist,doAna,flag_pot,kerd,dense) 
%{
Construct parbd for Fparams in the mobility solver code.
Note that this pre-computes points, so this is a place for memory optimization if we're
struggling on this aspect.

Inputs
    equ_radii   - (double) n_b x 1 array of equatorial radii for spheres/spheroids
    polar_radii - (double) n_b x 1 array of polar radii for spheres/spheroids
    p           - (int) spheroidal harmonic order
    C           - (double) n_b x 3 array of centers 
    collision_eps - (double) collision buffer/tolerance
    mdist       - (double) collision buffer for body-body interactions
    doAna       - (bool) indicates whether density that is passed in needs to be transformed to harmonic space
    flag_pot    - (string) type of Stokes potential to be used for problem in mobility oslver
    kerd        - (int) dimension of kernel (should be 1 for Laplace potentials, 3 for Stokes potentials)
    dense       - (bool) whether to use FMM or not
    tsl_backend - (string, optional) 'spheroidal' (default) or 'cartesian'
    tsl_dealiasing - (bool, optional) enable dealiasing in TSL near-evaluation (default true)
    tsl_dealiasing_pad - (int, optional) padding for dealiasing

Outputs
    parbd       - (struct) struct with rigid body parameters; see spheroidal_mobility for properties.
%}

arguments
    opts.equ_radii (:,1) double 
    opts.polar_radii (:,1) double
    opts.p (1,1) double
    opts.C (:,3) double
    opts.collision_eps (1,1) double = NaN
    opts.mdist (1,1) double
    opts.doAna (1,1) logical
    opts.flag_pot (1,1) string {mustBeMember(opts.flag_pot,["SL_Stk_3D", "TSL_Stk_3D"])}
    opts.kerd (1,1) double
    opts.dense (1,1) logical
    opts.bodydist (1,1) struct
    opts.tsl_backend (1,1) string = "spheroidal"
    opts.tsl_dealiasing (1,1) logical = true
    opts.tsl_dealiasing_pad (1,1) double = 4
    opts.MRot = []
end

% Alias variables
equ_radii = opts.equ_radii;
polar_radii = opts.polar_radii;
p = opts.p;
C = opts.C;
collision_eps = opts.collision_eps;
mdist = opts.mdist;
doAna = opts.doAna;
flag_pot = opts.flag_pot;
kerd = opts.kerd;
dense = opts.dense;
bodydist = opts.bodydist;
tsl_backend = opts.tsl_backend;
tsl_dealiasing = opts.tsl_dealiasing;
tsl_dealiasing_pad = opts.tsl_dealiasing_pad;
MRot = LOCAL_normalize_rotations(opts.MRot, size(C,1));

np=2*p*(p+1); 
Nb = kerd*np; % DOF per particle   
n3 = size(C,1); % no of particles
N = kerd*np*n3; % DOF total

% Smooth quadrature weights (GL x Trapezoidal)
[~, gwt]=g_grid(p+1);
wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
wt = wt(:);

X = cell(n3,1); W=cell(n3,1); Nr=cell(n3,1);
tau = cell(n3,1);
Cg = reshape(repmat(C.',np,1),3,[]).';
Xg = zeros(np*n3,3); Wg=zeros(np*n3,1); Nrg=Xg; Xrp = Xg;

body_shape_types = calculate_shape_type(equ_radii, polar_radii);

for j=1:n3
    shape_type = body_shape_types(j);
    equatorial_radius = equ_radii(j);
    polar_radius = polar_radii(j);
    Rj = MRot{j};
    surface = LOCAL_build_axisymmetric_shape(p, shape_type, equatorial_radius, polar_radius);
    
    % Pre-compute points of shape
    X_local = [surface.cart.x, surface.cart.y, surface.cart.z];
    X{j} = X_local*Rj.';

    % Pre-compute normals
    Nr_local = [surface.geoProp.nor.x, surface.geoProp.nor.y, surface.geoProp.nor.z];
    Nr{j} = Nr_local*Rj.';

    % Pre-compute area elements
    W{j} = surface.geoProp.W .* wt;

    % Pre-compute moment of inertia matrix
    Xj = X{j};
    Wj = W{j}(:);
    stau = sum(Wj .* sum(Xj.^2, 2)); % \sum W_i \|x_i\|^2
    S = Xj.' * (Xj .* Wj);           % \sum W_i x_i x_i^T
    tau{j} = stau * eye(3) - S;

    % Assign to correct block
    indx=(1:np)+np*(j-1);
    Xrp(indx,:) = X{j};
    Xg(indx,:) = Xrp(indx,:)+Cg(indx,:);  
    Wg(indx,:) = W{j}; 
    Nrg(indx,:) = Nr{j}; 
end

%{
Let's restrict ourselves to the case of kerd = 3 (e.g. Stokes).

The point of duplicating the points here is that each geometric point has
three DOFs: the x-, y-, and z-components of the density. Duplicating here
allows for easier manipulation of this fact (hence, the purpose of ci and cj
below).
%}
Xv = reshape(repmat(Xg,1,kerd)',3,[])'; 
Nrv = reshape(repmat(Nrg,1,kerd)',3,[])';    
Wv = repmat(Wg,1,kerd)'; Wv = Wv(:);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Form and fill params struct
parbd = struct(...
    'flag_pot',flag_pot, ...
    'kh',0, ...
    'transinv',1, ...
    'sym',1, ...
    'proxy',0,...
    'keval','points', ...
    'dim',3, ...
    'mu',1, ...
    'Xrp',Xrp, ...
    'Xp',Xg, ...
    'Nrp',Nrg, ...
    'X',Xv, ...
    'nor',Nrv, ...
    'Wg',Wg, ...
    'W2',Wv.', ...
    'p',p, ...
    'np',np, ...
    'kerd',kerd, ...
    'n3',n3, ...
    'Nb',Nb, ...
    'N',N, ...
    'equ_radii',equ_radii, ...
    'polar_radii',polar_radii, ...
    'dense',dense, ...
    'C',C, ...
    'mdist',mdist, ...
    'collision_eps',collision_eps, ...
    'a',0, ...
    'doAna',doAna, ...
    'shape_type',body_shape_types, ...
    'bodydist', bodydist, ...
    'MRot', {MRot}, ...
    'tsl_backend', tsl_backend, ...
    'tsl_dealiasing', tsl_dealiasing, ...
    'tsl_dealiasing_pad', tsl_dealiasing_pad ...
); 

parbd.W = W; 
parbd.tau = tau; 

if kerd>1 % For Kernel_Eval?
    parbd.ci = repmat((1:kerd)',n3*np,1); 
    parbd.cj = parbd.ci; 
end

end

function surface = LOCAL_build_axisymmetric_shape(p, shape_type, equatorial_radius, polar_radius)
    switch shape_type
        case 'sphere'
            % The matvec will scale back to a radius 1 sphere anyways, so 
            % let's just build the desired sphere for consistency sake.
            surface = SurfaceSph(equatorial_radius*shape_gallery(p,''));
        case 'prolate'
            [u0, a] = calculate_u0_and_a_from_radii(shape_type, equatorial_radius, polar_radius);
            surface = SurfaceSph(prolate_spheroid_shape(p, u0, a));
        case 'oblate'
            [u0, a] = calculate_u0_and_a_from_radii(shape_type, equatorial_radius, polar_radius);
            surface = SurfaceSph(oblate_spheroid_shape(p, u0, a));
    end
end

function MRot = LOCAL_normalize_rotations(raw_rotations, n3)
    if isempty(raw_rotations)
        MRot = repmat({eye(3)}, n3, 1);
        return;
    end

    if iscell(raw_rotations)
        if numel(raw_rotations) ~= n3
            error('opts.MRot must contain one 3x3 matrix per body.');
        end
        MRot = reshape(raw_rotations, [], 1);
    elseif isnumeric(raw_rotations) && isequal(size(raw_rotations), [3 3 n3])
        MRot = cell(n3, 1);
        for k = 1:n3
            MRot{k} = raw_rotations(:,:,k);
        end
    elseif isnumeric(raw_rotations) && isequal(size(raw_rotations), [3 3]) && n3 == 1
        MRot = {raw_rotations};
    else
        error('opts.MRot must be empty, a 1-by-n3 or n3-by-1 cell array, or a 3x3xn3 numeric array.');
    end
end

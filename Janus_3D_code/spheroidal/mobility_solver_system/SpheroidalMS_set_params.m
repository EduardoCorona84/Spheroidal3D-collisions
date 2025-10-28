function parbd = SpheroidalMS_set_params(equ_radii,polar_radii,p,C,eps,mdist,out,doAna,flag_pot,kerd,dense) 
%{
Construct parbd for Fparams in the mobility solver code.
Note that this pre-computes points, so this is a place for memory optimization if we're
struggling on this aspect.

Inputs
    equ_radii   - (double) n_b x 1 array of equatorial radii for spheres/spheroids
    polar_radii - (double) n_b x 1 array of polar radii for spheres/spheroids
    shape_type  - (string) n_b x 1 array of strings: should be 'prolate', 'oblate', or 'sphere'.
    p           - (int) spheroidal harmonic order (bodies) 
    C           - (double) n_b x 3 array of centers 
    eps         - (double) epsilon buffer (collision dist)
    mdist       - (double) collision buffer for body-body interactions
    out         - (bool) external vs internal evaluation (set to 1) 
    doAna       - (bool) indicates whether density that is passed in needs to be transformed to harmonic space
    flag_pot    - (string) type of Stokes potential to be used for problem in mobility oslver
    kerd        - (int) dimension of kernel (should be 1 for Laplace potentials, 3 for Stokes potentials)
    dense       - (bool) whether to use FMM or not

Outputs
    parbd       - (struct) struct with rigid body parameters; see spheroidal_mobility for properties.
%}

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

shape_type = LOCAL_build_shape_types(equ_radii, polar_radii);

for j=1:n3
    shape_type = shape_type(j);
    equatorial_radius = equ_radii(j);
    polar_radius = polar_radii(j);
    surface = LOCAL_build_axisymmetric_shape(p, shape_type, equatorial_radius, polar_radius);
    
    % Pre-compute points of shape
    X{j} = [surface.cart.x, surface.cart.y, surface.cart.z];

    % Pre-compute normals
    Nr{j} = [surface.geoProp.nor.x, surface.geoProp.nor.y, surface.geoProp.nor.z];

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
    'out',out, ...
    'mdist',mdist, ...
    'eps',eps, ...
    'a',0, ...
    'doAna',doAna, ...
    'shape_type',shape_type ...
); 

parbd.W = W; 
parbd.tau = tau; 

if kerd>1 % For Kernel_Eval?
    parbd.ci = repmat((1:kerd)',n3*np,1); 
    parbd.cj = parbd.ci; 
end

end

function shape_type = LOCAL_build_shape_types(equ_radii, polar_radii)
    EQUALITY_TOL = 1e-14;
    ns = size(equ_radii, 1);
    shape_type = strings(ns, 1);
    for j=1:ns
        if abs(equ_radii - polar_radii) < EQUALITY_TOL
            shape_type(j) = 'sphere';
        elseif equ_radii > polar_radii
            shape_type(j) = 'oblate';
        else
            shape_type(j) = 'prolate';
        end
    end
end

function surface = LOCAL_build_axisymmetric_shape(p, shape_type, equatorial_radius, polar_radius)
    switch shape_type
        case 'sphere'
            error('not implemented');
        case 'prolate'
            [u0, a] = calculate_u0_and_a_from_radii(shape_type, equatorial_radius, polar_radius);
            surface = SurfaceSph(prolate_spheroid_shape(p, u0, a));
        case 'oblate'
            [u0, a] = calculate_u0_and_a_from_radii(shape_type, equatorial_radius, polar_radius);
            surface = SurfaceSph(oblate_spheroid_shape(p, u0, a));
    end
end
restoredefaultpath;
this_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(this_dir)));
addpath(repo_root, '-begin');

addpath(genpath(fullfile(repo_root, 'spheroidal')), '-begin');
addpath(genpath(fullfile(repo_root, 'support')), '-begin');
addpath(genpath(fullfile(repo_root, 'FMMLIB')), '-begin');
addpath(genpath(fullfile(repo_root, 'LCPsolvers')), '-begin');

rbs_code_path = fullfile(repo_root, 'Phoretic_Particles', 'Code', 'RBS_Code');
if exist(rbs_code_path, 'dir') == 7
    rmpath(genpath(rbs_code_path));
end
clear SSph_MatVec;

%% SETUP
scenario = '4bodytetrahedron';

% Simulation parameters
kerd = 1; % Kernel dimension
Nt = 10;
dt = 0.1;
denseMV = true; % Dense vs FMM off diagonal 
tdisc = 'euler';
mdist = 1e-1; % Collision parameter ?
collision_eps = 1e-1;

% Janus particle settings
lambda = 1;
gamma = 1;
% Use the attractive convention by default.
flip_surface_label = false;

% f(X)=\frac{1}{2} \frac{r \cos \theta}{r^2}+\frac{1}{2}=\frac{1}{2} \frac{\cos \theta}{r}+\frac{1}{2}
% Should really be 0.5*(X*y'/|X| + 1)...
boundary_label = @(X,y) 0.5*(X*y'./sqrt(sum(X.^2,2)) + 1);

% Spheroid parameters
p  = 8; % Spharm degree p

switch scenario
    case '1body'
        n3 = 1; % Number of objects
        equ_radii = [1/2];
        polar_radii = [1];
        C = [0 0 0]; % object centers
        oblate = false;
    case '2body'
        n3 = 2;
        equ_radii = [1/2 2/3];
        polar_radii = [1 1];
        C = [...
            0 1.5 0; ...
            0 -1.5 0; ...
        ];
        oblate = [false false];
    case '3bodystable'
        n3 = 3;
        equ_radii = [1/2 1/2 1/2];
        polar_radii = [1 1 1];
        % Equilateral triangle
        C = [...
            0 2*sqrt(3)/4 0; ...
            1 -2*sqrt(3)/4 0; ...
            -1 -2*sqrt(3)/4 0; ...
        ];
        oblate = [false false false];
    case '4bodytetrahedron'
        n3 = 4;
        equ_radii = 0.5*ones(1, n3);
        polar_radii = ones(1, n3);
        tetra_edge = 4;
        tetra_scale = tetra_edge/(2*sqrt(2));
        % Regular tetrahedron with centroid at the origin.
        C = tetra_scale * [ ...
             1  1  1; ...
             1 -1 -1; ...
            -1  1 -1; ...
            -1 -1  1; ...
        ];
        oblate = false(1, n3);
    otherwise
        error('Scenario not supported.');
end

% Point each spheroid's polar axis toward the origin.
body_axis_dir = -C;
init_norm = sqrt(sum(body_axis_dir.^2,2));
zero_dir = init_norm < eps; % Any spheroids sufficiently close to the origin?
if any(zero_dir)
    body_axis_dir(zero_dir,:) = repmat([0 0 1], sum(zero_dir), 1);
    init_norm = sqrt(sum(body_axis_dir.^2,2));
end
body_axis_dir = body_axis_dir./repmat(init_norm,1,3);

initial_rotation_matrices = LOCAL_align_z_axis_to_directions(body_axis_dir);
init_dir = repmat([0 0 1], n3, 1);

if flip_surface_label
    init_dir = -init_dir;
end

% Linear solver parameters
tol = 1e-4;

% Parameters for collision algorithm
bodydist = struct( ...
    'algo','moving_balls', ...
    'max_iter', 100, ...
    'tol',1e-6 ...
);

% Parameters for initial setup of bodies
parbd = struct( ...
    'n3',n3, ...
    'equ_radii',equ_radii, ...
    'polar_radii',polar_radii, ...
    'shape_type',oblate, ...
    'p',p, ...
    'Ct',C, ...
    'mdist',mdist, ...
    'collision_eps',collision_eps, ...
    'bodydist',bodydist ...
);

% Parameters for linear solvers
parslv = struct(...
    'solver','gmres', ...
    'tol',tol, ...
    'maxit',200, ...
    'rst',4, ...
    'precond_type','bkdiag', ...
    'prec',[], ...
    'colsolver','BBPGD', ...
    'coltol',1e-4, ...
    'colmaxit',100, ...
    'col_tolrel',tol, ...
    'col_tolabs',0.1*tol ...
);

Fparams = struct(...
    'parbd',parbd, ...
    'parslv',parslv, ...
    'Nt',Nt, ...
    'dt',dt, ...
    'comp',true, ...
    'type','JanusAmp', ...
    'lambda', lambda, ...
    'gamma', gamma, ...
    'init_dir', init_dir, ...
    'boundary_label', boundary_label, ...
    'denseMV',denseMV, ...
    'typeMV','SSph', ...
    'tdisc',tdisc, ...
    'initial_rotation_matrices', {initial_rotation_matrices} ...
);

% Plot parameters
Fparams.plotFlag = false;
Fparams.plotTrajectories = true;
Fparams.plotTrajMaxPoints = 200;
Fparams.plotSurfaceAlpha = 0.3;
Fparams.plotGrid = true;
Fparams.plotColor = 'sigma';
Fparams.plotColorMode = 'l2';
Fparams.plotForceVectors = true;
Fparams.plotView = [115 30];

spheroidal_mobility('JanusAmp', Fparams, []);

function rotations = LOCAL_align_z_axis_to_directions(directions)
    n3 = size(directions, 1);
    rotations = cell(1, n3);
    ez = [0 0 1];

    for body_idx = 1:n3
        target_dir = directions(body_idx,:);
        target_dir = target_dir / norm(target_dir);
        c = dot(ez, target_dir);

        if c > 1 - 1e-12
            rotations{body_idx} = eye(3);
            continue;
        end

        if c < -1 + 1e-12
            rotations{body_idx} = [1 0 0; 0 -1 0; 0 0 -1];
            continue;
        end

        v = cross(ez, target_dir);
        s = norm(v);
        vx = [ ...
            0, -v(3), v(2); ...
            v(3), 0, -v(1); ...
            -v(2), v(1), 0 ...
        ];

        rotations{body_idx} = eye(3) + vx + ((1 - c)/(s^2))*(vx*vx);
    end
end

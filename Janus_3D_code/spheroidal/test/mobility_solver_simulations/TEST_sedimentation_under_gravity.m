%{
Test of spheroids in free Stokes flow subjected to gravity.

One spheroid should not induce any rotation, but multiple spheroids might
(even if they don't collide) due to pressure disturbances.
%}

addpath(genpath('../../../spheroidal'));
addpath(genpath('../../../support'));

%% SETUP
scenario = '2body';

% Simulation parameters
kerd = 3; % Kernel dimension
dense = true; % Dense vs FMM off diagonal 
out = true; % outside vs inside sphere
Nt = 100;
dt = 0.1;
denseMV = true;
tdisc = 'euler';
p = 8;
mdist = 1e-2; % Collision parameter ?

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
            -2 0 0; ...
            2 0 0; ...
        ];
        oblate = [false false];
    otherwise
        error('Scenario not supported.');
end

% Linear solver parameters
tol = 1e-4;

% Force and torque functions (let gravity be proportional to volume)
[u0, a] = calculate_u0_and_a_from_radii([], equ_radii, polar_radii);
body_volumes = LOCAL_calculate_volume(n3, p, u0, a, oblate);
Tfun = @(t,C,q) [20*(body_volumes.').^3; zeros(2,n3) ];
% Ffun = @(t,C,q) zeros(3,n3);
% Tfun = @(t,C,q) zeros(3,n3);
Ffun = @(t,C,q) [zeros(2,n3) ; -9.81*(body_volumes.').^3];

% Parameters for initial setup of bodies
parbd = struct( ...
    'n3',n3, ...
    'equ_radii',equ_radii, ...
    'polar_radii',polar_radii, ...
    'shape_type',oblate, ...
    'p',p, ...
    'Ct',C, ...
    'mdist',mdist, ...
    'eps',eps, ...
    'out',out ...
);

% Parameters for linear solvers
parslv = struct(...
    'solver','gmres', ...
    'tol',tol, ...
    'maxit',200, ...
    'rst',4, ...
    'prtype','bkdiag', ...
    'prec',[], ...
    'colsolver','BBPGD', ... %% Collision parameters will need to be changed after update.
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
    'type','FTfun', ...
    'denseMV',denseMV, ...
    'typeMV','SSph', ...
    'tdisc',tdisc, ...
    'Tfun',Tfun, ...
    'Ffun',Ffun, ...
    'plotFlag',true ...
);

% Plot parameters
Fparams.plotFlag = true;
Fparams.plotTrajectories = true;
Fparams.plotTrajMaxPoints = 200;
Fparams.plotSurfaceAlpha = 0.3;
Fparams.plotGrid = true;
Fparams.plotColor = 'sigma';
Fparams.plotColorMode = 'l2';


spheroidal_mobility('', Fparams, []);

%% Utility functions
function volumes = LOCAL_calculate_volume(ns, p, u0, a, oblate)
    volumes = zeros(1, ns);
    for i=1:ns
        if oblate(i)
            pts = oblate_spheroid_shape(p, u0(i), a(i));
        else
            pts = prolate_spheroid_shape(p, u0(i), a(i));
        end
        s = SurfaceSph(pts);
        volumes(i) = s.volume;
    end
end
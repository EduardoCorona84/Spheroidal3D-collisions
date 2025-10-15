%{
Sedimentation test for Stokesian suspension of n^3 spherical rigid bodies 
inside a spherical shell.
 
Inputs: 
fname - (string) filename for experiment info

The following are converted from strings when necessary: 

Body centers, radii, parameters
n     - (int)    cubic lattice is n x n x n
rd    - (double) minimum radius 
Cdst  - (double) distance between spheres in lattice
p     - (int)    spherical harmonic order (bodies) 
ep    - (double) epsilon buffer (collision dist)  

Time discretization
Nt    - (int)    number of timesteps
dt    - (double) timestep length
tdisc - (string) timestepping scheme (euler,trapz,rk4)

lambda - (double) mod lap parameter 
%}
function [Fparams]=Test_ModLap_Mobility_Amphi(n,rd,Cdst,p,ep,Nt,dt,tdisc,lambda,saveLCPs)
%% Default parameters
if ~exist('p','var') || isempty(p)
    p=8; 
end
if ~exist('lambda','var') || isempty(lambda)
    lambda=0.1;
end
if ~exist('rd','var') || isempty(rd)
    rd=1;
end
if ~exist('n','var') || isempty(n)
    n=5;
end
if ~exist('Cdst','var') || isempty(Cdst)
    Cdst=2.3; 
end
if ~exist('ep','var') || isempty(ep)
    ep=.3;
end
if ~exist('Nt','var') || isempty(Nt)
    Nt=500;
end
if ~exist('dt','var') || isempty(dt)
    dt=.1;
end
if ~exist('tdisc','var') || isempty(tdisc)
    tdisc='euler';
end
if ~exist('saveLCPs','var') || isempty(saveLCPs)
    saveLCPs=true; 
end
%% misc extra parameters
tol=1e-4;
mdist=3; 
denseMV=false; 
denseforce=1;
gamma=1; 
%% Files to save results
mfilePath = mfilename('fullpath');
if contains(mfilePath,'LiveEditorEvaluationHelper')
    mfilePath = matlab.desktop.editor.getActiveFilename;
end
[basedir,~,~] = fileparts(mfilePath);
lcpDataDir = fullfile(basedir, 'data');
mkdir(lcpDataDir)
postFix = ['.n_' num2str(n) '.p_' num2str(p) '.cDist_' num2str(Cdst)];
fname=fullfile(lcpDataDir, ['amphi' postFix]);
lcpResDir = fullfile(basedir, 'LCPsolvers/data');
mkdir(lcpResDir);
LCP_file_path=fullfile(lcpResDir, ['amphi' postFix]);
%% boundary_label function
boundary_label =  @(X,y) 0.5*X*y'./sqrt(sum(X.^2,2)).^2 + 1/2;
%% Make sure all the code is on the matlabpath
% Remove addpaths if compiling in command line (mcc)
addpath(basedir);
addpath(fullfile(basedir,'support'));
addpath(genpath(fullfile(basedir, 'LCPsolvers/solvers')))
addpath(fullfile(basedir, 'FMMLIB/fmmlib3d-1.2/matlab'));
addpath(fullfile(basedir,'FMMLIB/stfmmlib3d-1.2/matlab'));
%% Create Fparams struct 
Fparams = struct('Nt',Nt,'dt',dt,'comp',1,'type','JanusAmp',...
    'lambda',lambda,'gamma',gamma,'denseMV',denseMV,...
    'typeMV','Vsh','tdisc',tdisc, ...
    'boundary_label',boundary_label,'denseforce',denseforce,...
    'saveLCPs',saveLCPs,'LCP_file_path',LCP_file_path);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fill Parameter Structs
% body parameters
n3 = n^3; 
rd=rd*ones(n3,1);
Fparams.parbd = struct('Shape','','n3',n3,'rd',rd,'diam',2*rd,'p',p,'mdist',mdist,'mxrd',rd(1),'eps',ep,'out',1);
%% Particle centers (cubic lattice in this example)
Cdst=mean(rd)*Cdst; % make center distance relative to radii
lx=0:Cdst:Cdst*(n-1);
lx = lx - mean(lx); 
[xx,yy,zz] = meshgrid(lx);
C = [xx(:) yy(:) zz(:)]; 
display(C);
% randomize centers and/or radii
try %#ok<TRYNC>
    rng('default'); 
end
while true
    C = C + 0.1*rand(size(C));
    [~,~,mindst,~] = LOCAL_check_collision_sph(C,Fparams);
    if all(mindst > ep/10 )
        break
    end
end
display(C); 
Fparams.parbd.Ct = C;

% initial particle orientations
%% NIC: set the initial direction to be towards the center 
init_dir= -C;
for i = 1:n3
    init_dir(i,:) = init_dir(i,:) / norm(init_dir(i,:));
end
Fparams.init_dir = init_dir;


% LCP solver parameters
Fparams.lcpOpts = struct(...
    'solver','proxquasinewton',...
    'max_iter',100,...
    'tol_rel',1e-6,...
    'tol_abs',1e-5, ...
    'stepSize',struct(...
        'init','uniform',...
        'kappa','uniform',...
        'eta','opt')...
);

% linear solver parameters
Fparams.parslv = struct('solver','gmres','tol',tol,'maxit',200,'rst',4,'prtype','bkdiag','prec',[],'prLCP',false); 

% low-fidelity parameters
% lofi_p = 2;
% Fparams.lofi = struct('Shape','','n3',n3,'rd',rd,'p',lofi_p,'Ct',C,'mdist',mdist,'eps',ep,'out',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Rigid Body Stokes 
RBS_mobility(fname,Fparams,true);
end
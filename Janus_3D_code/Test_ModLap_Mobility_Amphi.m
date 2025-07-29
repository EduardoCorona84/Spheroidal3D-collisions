function [Fparams]=Test_ModLap_Mobility_Amphi(n,rd,Cdst,p,ep,Nt,dt,tdisc,lambda,boundary_label,saveLCPs)
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
%% Files to save results
mfilePath = mfilename('fullpath');
if contains(mfilePath,'LiveEditorEvaluationHelper')
    mfilePath = matlab.desktop.editor.getActiveFilename;
end
[mfilePath,~,~] = fileparts(mfilePath);
resultsDir = fullfile(mfilePath, 'results');
mkdir(resultsDir)
postFix = ['.n_' num2str(n) '.p_' num2str(p) '.cDist_' num2str(Cdst)];
fname=fullfile(resultsDir, ['amphi' postFix]);
lcpResDir = fullfile(mfilePath, 'LCPsolvers/results');
mkdir(lcpResDir);
LCP_file_path=fullfile(lcpResDir, ['amphiLCPs' postFix]);
%% Make sure all the code is on the matlabpath'
%Remove addpaths if compiling in command line (mcc)
addpath ./; 
addpath ./support; 
addpath(genpath(fullfile(mfilePath, 'LCPsolvers/solvers')))
addpath ./FMMLIB/fmmlib3d-1.2/matlab/;
addpath ./FMMLIB/stfmmlib3d-1.2/matlab/;

%% Particle centers (cubic lattice in this example)
Cdst=mean(rd)*Cdst; %make center distance relative to radii
lx=0:Cdst:Cdst*(n-1);
lx = lx - mean(lx); 
[xx,yy,zz] = meshgrid(lx);
C = [xx(:) yy(:) zz(:)]; 
display(C);

% randomize centers and/or radii
try %#ok<TRYNC>
    rng('default'); 
end
C = C + 0.1*rand(size(C));
display(C); 
n3 = size(C,1); 
rd=rd*ones(n3,1);

% initial particle orientations
%% NIC: set the initial direction to be towards the center 
init_dir= -C;
for i = 1:n3
    init_dir(i,:) = init_dir(i,:) / norm(init_dir(i,:));
end
%% This is the old setting
% nrt = [zeros(n3/2,2) ones(n3/2,1)];
% init_dir=[nrt;-nrt];
%% Dont know who did this
%init_dir=rand(n3,3);
%init_dir=init_dir./repmat(sqrt(init_dir(:,1).^2+init_dir(:,2).^2+init_dir(:,3).^2),1,3);

% misc extra parameters
tol=1e-4; 
mdist=3; 
denseMV=true; 
denseforce=1;
gamma=1; 

% default boundary label 
if nargin<11
    boundary_label = @(X,y) 0.5*X*y'./sqrt(sum(X.^2,2)).^2 + 1/2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fill Parameter Structs
% body parameters
parbd = struct('Shape','','n3',n3,'rd',rd,'p',p,'Ct',C,'mdist',mdist,'eps',ep,'out',1);

% LCP solver parameters
lcpOpts = struct('solver','bbpgd','max_iter',100,'tol_rel',1e-12,'tol_abs',1e-12);

% linear solver parameters
parslv = struct('solver','gmres','tol',tol,'maxit',200,'rst',4,'prtype','bkdiag','prec',[]);  

%Create Fparams struct 
Fparams = struct('parbd',parbd,'parslv',parslv,'lcpOpts',lcpOpts,...
    'Nt',Nt,'dt',dt,'comp',1,'type','JanusAmp','lambda',lambda,'gamma',gamma,...
'denseMV',denseMV,'typeMV','Vsh','tdisc',tdisc,'init_dir',init_dir,'boundary_label',boundary_label,'denseforce',denseforce,...
'saveLCPs',saveLCPs,'LCP_file_path',LCP_file_path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Rigid Body Stokes 
RBS_mobility(fname,Fparams,[]);
end
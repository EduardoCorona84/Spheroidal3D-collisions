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
function [Fparams]=Test_ModLap_Mobility_Amphi(n,rd,Cdst,p,ep,Nt,dt,tdisc,lambda,saveLCPs,initMode, tol,mdist,denseMV,denseforce,gamma,loadIntermediate, plotFlag,polydisperseRatio)
%% Default parameters
if ~exist('p','var') || isempty(p)
    p=2; 
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
    Nt=200;
end
if ~exist('dt','var') || isempty(dt)
    dt=.5;
end
if ~exist('tdisc','var') || isempty(tdisc)
    tdisc='euler';
end
if ~exist('saveLCPs','var') || isempty(saveLCPs)
    saveLCPs=true; 
end
if ~exist('initMode','var') || isempty(initMode)
    initMode='lattice'; 
end
if ~exist('tol','var') || isempty(tol)
    tol=1e-4;
end
if ~exist('mdist','var') || isempty(mdist)
    mdist=3; 
end
if ~exist('denseMV','var') || isempty(denseMV)
    denseMV=true; 
end
if ~exist('denseforce','var') || isempty(denseforce)
    denseforce=1;
end
if ~exist('gamma','var') || isempty(gamma)
    gamma=1; 
end
if ~exist('loadIntermediate','var') || isempty(loadIntermediate)
    loadIntermediate=false; 
end
if ~exist('plotFlag','var') || isempty(plotFlag)
    plotFlag=true; 
end
if ~exist('polydisperseRatio','var') || isempty(polydisperseRatio)
    polydisperseRatio=0.2; 
end
%% boundary_label function
boundary_label =  @(X,y) 0.5*X*y'./sqrt(sum(X.^2,2)).^2 + 1/2;
%% Files to save results
mfilePath = mfilename('fullpath');
if contains(mfilePath,'LiveEditorEvaluationHelper')
    mfilePath = matlab.desktop.editor.getActiveFilename;
end
[dirname,~,~] = fileparts(mfilePath);
basedir = fullfile(dirname, '..');
lcpDataDir = fullfile(basedir, 'data');
postFix = ['.' initMode '.n_' num2str(n) '.p_' num2str(p) '.cDist_' num2str(Cdst)];
lcpResDir = fullfile(basedir, 'LCPsolvers/data');
fname = fullfile(lcpDataDir, ['amphi' postFix]); % Name of file for regular results file
LCP_file_path = fullfile(lcpResDir, ['amphi' postFix]); % Name of the LCP results file
% Make directories if they do not exist 
mkdir(lcpDataDir)
mkdir(lcpResDir);
%% Make sure all the code is on the matlabpath
% Remove addpaths if compiling in command line (mcc)
addpath(basedir);
addpath(fullfile(basedir,'test'));
addpath(fullfile(basedir,'support'));
addpath(genpath(fullfile(basedir, 'LCPsolvers/solvers')))
addpath(fullfile(basedir, 'FMMLIB/fmmlib3d-1.2/matlab'));
addpath(fullfile(basedir,'FMMLIB/stfmmlib3d-1.2/matlab'));
%% Create Fparams struct 
Fparams = struct('Nt',Nt,'dt',dt,'comp',1,'type','JanusAmp',...
    'lambda',lambda,'gamma',gamma,'denseMV',denseMV,...
    'typeMV','Vsh','tdisc',tdisc, ...
    'boundary_label',boundary_label,'denseforce',denseforce,...
    'saveLCPs',saveLCPs,'LCP_file_path',LCP_file_path, ...
    'loadIntermediate',loadIntermediate);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initialize configuration
switch initMode
    case 'lattice'
        [C, rd, init_dir] = init_lattice(n, Cdst, rd, polydisperseRatio);
    case 'vesicle'
        [C, init_dir] = init_vesicle(n, Cdst);
    case 'special'
        r_ = 10;
        [C_, rd, init_dir] = init_lattice(n, Cdst, rd, polydisperseRatio);
        C = [C_ + repmat([r_ r_ r_], n^3,1);
             C_ + repmat([-r_ r_ r_],n^3,1);
             C_ + repmat([r_ -r_ r_],n^3,1);
             C_ + repmat([r_ r_ -r_],n^3,1);
             C_ + repmat([-r_ -r_ r_],n^3,1);
             C_ + repmat([-r_ r_ -r_],n^3,1);
             C_ + repmat([r_ -r_ -r_],n^3,1);
             C_ + repmat([-r_ -r_ -r_],n^3,1);
             ];
        rd = repmat(rd, 8,1);
        init_dir = repmat(init_dir, 8,1);
    otherwise
        error([initMod ' not a recognized initMode'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fill Parameter Structs
% body parameters
n3 = size(C,1); 
Fparams.plotFlag = plotFlag;
Fparams.parbd = struct('Shape','','n3',n3,'rd',rd,'diam',2*rd,'p',p,'mdist',mdist,'mxrd',rd(1),'eps',ep,'out',1);
Fparams.parbd.Ct = C;
Fparams.init_dir = init_dir;
% LCP solver parameters
Fparams.lcpOpts = defaultLCPOpts(struct(...
    'solver','proxquasinewton', ...
    'max_iter',1000, ...
    'kkt_rel',1e-8, ...
    'kkt_abs',1e-8, ...
    'warmStart',true));

% linear solver parameters
Fparams.parslv = struct('solver','gmres','tol',tol,'maxit',200,'rst',4,'prtype','bkdiag','prec',[],'prLCP',false); 

% low-fidelity parameters
% lofi_p = 2;
% Fparams.lofi = struct('Shape','','n3',n3,'rd',rd,'p',lofi_p,'Ct',C,'mdist',mdist,'eps',ep,'out',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Rigid Body Stokes 
RBS_mobility(fname,Fparams);
end
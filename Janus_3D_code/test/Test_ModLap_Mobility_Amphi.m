%{
Sedimentation test for Stokesian suspension of n^3 spherical rigid bodies 
inside a spherical shell.
 
Inputs: 
fname - (string) filename for experiment info

The following are converted from strings when necessary: 

Body centers, radii, parameters
n     - (int)    cubic lattice is n x n x n
meanRadius    - (double) mean radius 
Cdst  - (double) distance between spheres in lattice
p     - (int)    spherical harmonic order (bodies) 
ep    - (double) epsilon buffer (collision dist)  
polydisperseRatio - (double) maximum purtubation relative to the meanRadius

LCP params
lcpSlvr - (string) {'bbpgd', 'l-bfgs-b', 'p-l-bfgs', 'proxquasinewton', 'bifi'}
lcpTol - (double) kkt conditions 
lcpMaxIter - (int) maximum iterations
lcpWarmStart - (bool) warm start from previous solution

Time discretization
Nt    - (int)    number of timesteps
dt    - (double) timestep length
tdisc - (string) timestepping scheme {'euler','trapz','rk4'}

parbd/parslv Params
gamma - (double)
lambda - (double) mod lap parameter 
mdist - (int)
denseMV - (bool) compute stokes mobility matrix densely only recommended for small systems p<= 4 n3<=20
denseforce - (bool) compute the modified laplace forces with dense matrices
gmresTol - (double) gmres relErr tolerence

Nuisance Params
saveLCPs - (bool) save the LCPs
initMode - (string) orientation of the initial config of particles {'lattice', 'vesicle', 'special'}
loadIntermediate - (bool) load results from previous run and continue from that point
plotFlag - (bool) plotting during simulation only recomended for small systems on local machine
seed - (int) random seed for reproducibility
%}
function [Fparams]=Test_ModLap_Mobility_Amphi(...
    n,meanRadius,Cdst,p,ep,polydisperseRatio, ...
	lcpSlvr,lcpTol,lcpMaxIter,lcpWarmStart, lcpPLo, lcpTolLo, ...
	Nt,dt,tdisc, ...
    gamma,lambda,mdist,denseMV,denseforce,gmresTol,...
    saveLCPs,initMode,loadIntermediate,plotFlag,seed)
%% Body default parameters
if ~exist('p','var') || isempty(p)
    p=2; 
end
if ~exist('meanRadius','var') || isempty(meanRadius)
    meanRadius=1;
end
if ~exist('n','var') || isempty(n)
    n=2;
end
if ~exist('Cdst','var') || isempty(Cdst)
    Cdst=3; 
end
if ~exist('ep','var') || isempty(ep)
    ep=.3;
end
if ~exist('polydisperseRatio','var') || isempty(polydisperseRatio)
    polydisperseRatio=0.0; 
end
%% LCP default params
if ~exist('lcpSlvr','var') || isempty(lcpSlvr)
    lcpSlvr='bbpgd'; 
end
if ~exist('lcpTol','var') || isempty(lcpTol)
    lcpTol=1e-6; 
end
if ~exist('lcpMaxIter','var') || isempty(lcpMaxIter)
    lcpMaxIter=100; 
end
if ~exist('lcpWarmStart','var') || isempty(lcpWarmStart)
    lcpWarmStart=true; 
end
if ~exist('lcpPLo','var') || isempty(lcpPLo)
    lcpPLo=6; 
end
if ~exist('lcpTolLo','var') || isempty(lcpTolLo)
    lcpTolLo=1e-6; 
end
%% Time disc default params
if ~exist('Nt','var') || isempty(Nt)
    Nt=200;
end
if ~exist('dt','var') || isempty(dt)
    dt=.1;
end
if ~exist('tdisc','var') || isempty(tdisc)
    tdisc='euler';
end
%% pardb/parslv default params
if ~exist('gamma','var') || isempty(gamma)
    gamma=1; 
end
if ~exist('lambda','var') || isempty(lambda)
    lambda=0.1;
end
if ~exist('mdist','var') || isempty(mdist)
    mdist=3; 
end
if ~exist('denseMV','var') || isempty(denseMV)
    denseMV=true; 
end
if ~exist('denseforce','var') || isempty(denseforce)
    denseforce=true;
end
if ~exist('gmresTol','var') || isempty(gmresTol)
    gmresTol=1e-6;
end
%% Nuisance default params
if ~exist('saveLCPs','var') || isempty(saveLCPs)
    saveLCPs=false; 
end
if ~exist('initMode','var') || isempty(initMode)
    initMode='lattice'; 
end
if ~exist('loadIntermediate','var') || isempty(loadIntermediate)
    loadIntermediate=false; 
end
if ~exist('plotFlag','var') || isempty(plotFlag)
    plotFlag=false; 
end
if ~exist('seed','var') || isempty(seed)
    seed=1; 
end
%% For repeatable behavior
rng(seed)
%% Files to save results
mfilePath = mfilename('fullpath');
if contains(mfilePath,'LiveEditorEvaluationHelper')
    mfilePath = matlab.desktop.editor.getActiveFilename;
end
[dirname,~,~] = fileparts(mfilePath);
basedir = fullfile(dirname, '..');
resDir = fullfile(basedir, 'resultsForRecord');
postFix = ['.' initMode '.n_' num2str(n) '.p_' num2str(p) ...
    '.cDist_' num2str(Cdst) '.lcpSlvr_' lcpSlvr ...
    '.polyDisperseRatio_' num2str(polydisperseRatio)];
fname = fullfile(resDir, ['amphi' postFix 'verification.new']); % Name of file for regular results file
LCP_file_path = fullfile(resDir, ['amphi.lcp' postFix]); % Name of the LCP results file
% Make directories if they do not exist 
mkdir(resDir)
%% Make sure all the code is on the matlabpath
% Remove addpaths if compiling in command line (mcc)
addpath(basedir);
addpath(fullfile(basedir,'test'));
addpath(fullfile(basedir,'support'));
addpath(genpath(fullfile(basedir, 'LCPsolvers/solvers')))
addpath(fullfile(basedir, 'FMMLIB/fmmlib3d-1.2/matlab'));
addpath(fullfile(basedir,'FMMLIB/stfmmlib3d-1.2/matlab'));

%% initialize configuration
switch initMode
    case 'lattice'
        [C, rd, init_dir] = init_lattice(n, Cdst, meanRadius, polydisperseRatio, ep);
    case 'vesicle'
        [C, init_dir] = init_vesicle(n, Cdst);
    case 'special'
        r_ = 10;
        [C_, rd, init_dir] = init_lattice(n, Cdst, meanRadius, polydisperseRatio);
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
%% Create Fparams struct 
% high level params
boundary_label =  @(X,y) 0.5*X*y'./sqrt(sum(X.^2,2)).^2 + 1/2;
Fparams = struct('Nt',Nt,'dt',dt,'comp',1,'type','JanusAmp',...
    'lambda',lambda,'gamma',gamma,'denseMV',denseMV,...
    'typeMV','Vsh','tdisc',tdisc, ...
    'boundary_label',boundary_label,'denseforce',denseforce,...
    'saveLCPs',saveLCPs,'LCP_file_path',LCP_file_path, ...
    'loadIntermediate',loadIntermediate);
Fparams.plotFlag = plotFlag;
Fparams.init_dir = init_dir;
% parbd - matVec params
n3 = size(C,1); 
Fparams.parbd = struct('Shape','','n3',n3,'rd',rd,'diam',2*rd,'p',p,...
    'mdist',mdist,'mxrd',max(rd),'eps',ep,'out',1,'dense',denseMV);
Fparams.parbd.Ct = C;
% parslv - linear solver parameters
Fparams.parslv = struct('solver','gmres','tol',gmresTol,'maxit',50,...
    'rst',4,'prtype','bkdiag','prec',[],'prLCP',false); 
% LCP solver parameters
Fparams.lcpOpts = defaultLCPOpts(struct(...
    'solver',lcpSlvr, ...
    'max_iter',lcpMaxIter, ...
    'kkt_rel',lcpTol, ...
    'kkt_abs',lcpTol, ...
    'warmStart',lcpWarmStart));
Fparams.lcpOpts.low.p=lcpPLo;
Fparams.lcpOpts.low.gmresTol=lcpTolLo;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Rigid Body Stokes 
disp('=======================================')
disp('Hyper Parameters')
fprintf('- n         = %d\n', n)
fprintf('- p         = %d\n', p)
fprintf('- gmresTol  = %.2g\n', gmresTol)
fprintf('- denseMV   = %d\n', denseMV)
fprintf('- polyRatio = %.2f\n', polydisperseRatio)
disp('=======================================')
RBS_mobility(fname,Fparams);
end

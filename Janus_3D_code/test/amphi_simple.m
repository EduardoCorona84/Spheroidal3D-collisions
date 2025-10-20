%% Files to save results
mfilePath = mfilename('fullpath');
if contains(mfilePath,'LiveEditorEvaluationHelper')
    mfilePath = matlab.desktop.editor.getActiveFilename;
end
[dirname,~,~] = fileparts(mfilePath);
basedir = fullfile(dirname, '..')
lcpDataDir = fullfile(basedir, 'data');
mkdir(lcpDataDir)
postFix = 'simple';
fname=fullfile(lcpDataDir, ['amphi' postFix]);
lcpResDir = fullfile(basedir, 'LCPsolvers/data');
mkdir(lcpResDir);
LCP_file_path=fullfile(lcpResDir, ['amphi' postFix]);
%% Make sure all the code is on the matlabpath
addpath(basedir);
addpath(fullfile(basedir,'support'));
addpath(genpath(fullfile(basedir, 'LCPsolvers/solvers')))
addpath(fullfile(basedir, 'FMMLIB/fmmlib3d-1.2/matlab'));
addpath(fullfile(basedir,'FMMLIB/stfmmlib3d-1.2/matlab'));
%% Set hyperparams
p=4; 
lambda=0.1;
rd=1;
n=4;
Cdst=2.5; % NIC: change back to 4
ep=.3;
Nt=500;
dt=.1;
tdisc='euler';
saveLCPs=true; 
tol=1e-4;
mdist=3; 
denseMV=true; 
denseforce=1;
gamma=1; 
loadIntermediate = true;
plotFlag = true;
%% Build the simplest configuration just 2 particles on a line headed straight for each other.
C = [-2,0,0; 
    2,0,0;
    0,-2,0; 
    0,2,0];

init_dir= -C;
n3 = size(C,1);
rd = ones(n3,1);
for i = 1:n3
    init_dir(i,:) = init_dir(i,:) / norm(init_dir(i,:));
    % init_dir(i,:) = init_dir(i,:) / 10;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fill Parameter Structs
% boundary_label function
boundary_label =  @(X,y) 0.5*X*y'./sqrt(sum(X.^2,2)).^2 + 1/2;
% body parameters
parbd = struct('Shape','','n3',n3,'rd',rd,'p',p,'Ct',C,'mdist',mdist,'eps',ep,'out',1);
% LCP solver parameters
lcpOpts = struct(...
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
parslv = struct('solver','gmres','tol',tol,'maxit',200,'rst',4,'prtype','bkdiag','prec',[],'prLCP',false); 
% Create Fparams struct 
Fparams = struct('parbd',parbd,'parslv',parslv,'lcpOpts',lcpOpts,...
    'Nt',Nt,'dt',dt,'comp',1,'type','JanusAmp','lambda',lambda,'gamma',gamma,...
    'denseMV',denseMV,'typeMV','Vsh','tdisc',tdisc,'init_dir',init_dir, ...
    'boundary_label',boundary_label,'denseforce',denseforce,...
    'saveLCPs',saveLCPs,'LCP_file_path',LCP_file_path, 'plotFlag', plotFlag, ...
    'loadIntermediate',loadIntermediate);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Rigid Body Stokes 
RBS_mobility(fname,Fparams);
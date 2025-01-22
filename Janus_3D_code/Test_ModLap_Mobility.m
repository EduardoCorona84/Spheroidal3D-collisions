function [Fparams]=Test_ModLap_Mobility(fname,n,rd,Cdst,p,ep,Nt,dt,tdisc,lambda,e_north,e_south)
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

%Remove addpaths if compiling in command line (mcc)
addpath ./; 
addpath ./support; 
addpath ./LCPsolvers/Num4LCP_MatLab; 
addpath ./LCPsolvers/Num4LCP_MatLab/ext;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Particle centers (cubic lattice in this example)
rng('default');
Cdst=mean(rd)*Cdst;
lx=0:Cdst:Cdst*(n-1);
lx = lx - mean(lx); 
[xx,yy,zz] = meshgrid(lx);
C = [xx(:) yy(:) zz(:)]; 
n3 = size(C,1); 

%random initial directions 
init_dir=rand(n3,3);
init_dir=init_dir./repmat(sqrt(init_dir(:,1).^2+init_dir(:,2).^2+init_dir(:,3).^2),1,3);

% initialize point charges (?) 
q=10*rand(n3,1)-5*ones(n3,1); %random array (n3x1) between (-5,5)
chrnrm = (rd/4).*rand(n3,1); %random array (n3x1) between (0,r/4)
chrdir = rand(n3,3); %random array (n3x3) in (0,1)

% normalize charge directions
chrdir = chrdir./repmat(sqrt(chrdir(:,1).^2+chrdir(:,2).^2+chrdir(:,3).^2),1,3);
% charge strength * direction
point_charges = repmat(chrnrm,1,3).*chrdir; 

% radii array and misc params
rd = rd*ones(n3,1); 
tol=1e-4; 
mdist=3; 
denseMV=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% body parameters
parbd = struct('Shape','','n3',n3,'rd',rd,'p',p,'Ct',C,'mdist',mdist,'eps',ep,'out',1);

% linear solver parameters
parslv = struct('solver','gmres','tol',tol,'maxit',200,'rst',4,'prtype','bkdiag','prec',[],...
    'colsolver','BBPGD','coltol',1e-4,'colmaxit',100,'col_tolrel',tol,'col_tolabs',0.1*tol);  

%Create Fparams struct 
Fparams = struct('parbd',parbd,'parslv',parslv,...
    'Nt',Nt,'dt',dt,'comp',1,'type','JanusEHD','lambda',lambda,...
'denseMV',denseMV,'typeMV','Vsh','tdisc',tdisc,'init_dir',init_dir,'epsilon_north',e_north,'epsilon_south',e_south);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fparams.point_charges=point_charges;
Fparams.init_charges=point_charges; 
Fparams.q=q;
%Run Rigid Body Stokes 
RBS_mobility(fname,Fparams,[]);

end
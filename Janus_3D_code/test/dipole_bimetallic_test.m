
%{
Sedimentation test for Stokesian suspension of n^3 spherical rigid bodies 
inside a spherical shell.

Here we place the charges in the interior of the spheres, positive and
negative at differing orientations.
 
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



rng('default');


p=4;
lambda=0.1;
e_north=1;
e_south=1;


%parameters for Fparams
rd=[1];
n=2;
Cdst=4;
fname='2_dipoles';
ep=.3;
Nt=100;
dt=0.1;
tdisc='euler';

%charge density on surface of sphere (scaled so that F,T~O(1))
charge_label = @(X,y) 25000*(X*y'./sqrt(sum(X.^2,2)).^2);


%center of the spheres
C=[-2 0 0 ; 2 0 0];
C=C+0.1*rand(size(C));

n3 = size(C,1); 
if size(rd,1) == 1
    rd=rd*ones(n3,1);
end

%sets direction of the charge north poles
init_dir=2*(rand(size(C))-0.5*ones(size(C)));

%direction of the epsilon north poles
eps_dir=[1 0 0; -1 0 0];

init_dir=init_dir./repmat(sqrt(init_dir(:,1).^2+init_dir(:,2).^2+init_dir(:,3).^2),1,3);
eps_dir=eps_dir./repmat(sqrt(eps_dir(:,1).^2+eps_dir(:,2).^2+eps_dir(:,3).^2),1,3);


tol=1e-4; mdist=3; denseMV=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% body parameters
parbd = struct('Shape','','n3',n3,'rd',rd,'p',p,'Ct',C,'mdist',mdist,'eps',ep,'out',1);

% linear solver parameters
parslv = struct('solver','gmres','tol',tol,'maxit',200,'rst',4,'prtype','bkdiag','prec',[],...
    'colsolver','BBPGD','coltol',1e-4,'colmaxit',100,'col_tolrel',tol,'col_tolabs',0.1*tol);  

%Create Fparams struct 
Fparams = struct('parbd',parbd,'parslv',parslv,...
    'Nt',Nt,'dt',dt,'comp',1,'type','JanusEHD','lambda',lambda,...
'denseMV',denseMV,'typeMV','Vsh','tdisc',tdisc,'init_dir',init_dir,'epsilon_north',e_north,'epsilon_south',e_south,'charge_label',charge_label,'eps_dir',eps_dir);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fparams.point_charges=point_charges;
%Fparams.init_charges=point_charges; 
%Fparams.q=q;
%Run Rigid Body Stokes 
RBS_mobility(fname,Fparams,[]);

end
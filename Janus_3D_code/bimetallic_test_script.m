

%Remove addpaths if compiling in command line (mcc)
addpath ./; 
addpath ./support; 
addpath ./LCPsolvers/Num4LCP_MatLab; 
addpath ./LCPsolvers/Num4LCP_MatLab/ext;
addpath ./FMMLIB/stfmmlib3d-1.2/matlab;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Particle centers (cubic lattice in this example)

rng('default');


clear;clc;
p=4;
lambda=0.01;
e_north=1.2;
e_south=1.2;

rd=[1,2]';
n=2;
Cdst=4;
fname='testing_big_metal_oi';
ep=.4;
Nt=1600;
dt=.1;
tdisc='euler';


Efield=100*[0 1 0];

C= [-2 0 0; 3 0 0];
%C = 1.5*[-3 -3 0; -1 -1 0; 1 1 0; 3 3 0];
lx=0:Cdst:Cdst*(n-1);
lx = lx - mean(lx); 
[xx,yy,zz] = meshgrid(lx);
 %C = [xx(:) yy(:) zz(:)]; 

% C=C([1:2:6],:);
%C=0.75*[-3 0 0; 3 0 0; 0 2 1; 0 -2 1];
%C=zeros(5,3);
%C(:,1)=(-6:3:6)';
%C=3*[1 0 0; -0.866 0.5 0; -0.866 -0.5 0];
%C=C([1:2],:);
%C=[2.5 0 0; -2.5 0 0];
%C=[-3 -3 0; 0 0 0; 3 3 0];

gamma = 10;
C=C+0.01*rand(size(C));

n3 = size(C,1); 
if size(rd,1) == 1
    rd=rd*ones(n3,1);
end

%init_dir=rand(n3,3);
%init_dir=rand(n3,3);


shrink = 0.75;
charges = zeros(2,3,n3);
%charge_dir=[1 0 0; -1 0 0];
%charge_dir= -C
charge_dir = repmat([1 1 0],n3,1);
%charge_dir = [1 0 0; 1 0 0];
%charge_dir = [1 0 0; 0 1 0; -1 0 0];

%charge_dir=2*(rand(n3,3)-0.5*ones(n3,3));
%charge_dir=rand(n3,3);
%charge_dir=2*(rand(n3,3)-0.5);

%charge_dir=C;
init_dir=charge_dir;
np = 2*p*(p+1);
chargelabel = @(X,n) atan(X*n')/(pi/4);
charge_init = Setparams(p,[0 0 0],1,'SL_LMOD_3D','',1,1,4,0,lambda);
charge_vec=zeros(np*n3,1);


charge_dir=charge_dir./repmat(sqrt(charge_dir(:,1).^2+charge_dir(:,2).^2+charge_dir(:,3).^2),1,3);
for j=1:n3
    charges(1,:,j)=shrink*charge_dir(j,:);
    charges(2,:,j)=-shrink*charge_dir(j,:);
    
end
   q=[ones(1,n3) ; -ones(1,n3)]; 

tol=1e-4; mdist=3; denseMV=0;
denseforce = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% body parameters
parbd = struct('Shape','','n3',n3,'rd',rd,'p',p,'Ct',C,'mdist',mdist,'eps',ep,'out',1);

% linear solver parameters
parslv = struct('solver','gmres','tol',tol,'maxit',200,'rst',4,'prtype','bkdiag','prec',[],...
    'colsolver','BBPGD','coltol',1e-2,'colmaxit',100,'col_tolrel',tol,'col_tolabs',0.1*tol);  

%Create Fparams struct 
Fparams = struct('parbd',parbd,'parslv',parslv,...
    'Nt',Nt,'dt',dt,'comp',1,'type','JanusEHD','lambda',lambda,'gamma',gamma,'charge_label',chargelabel,...
'denseMV',denseMV,'typeMV','Vsh','tdisc',tdisc,'init_dir',init_dir,...
'epsilon_north',e_north,'epsilon_south',e_south,'Efield',Efield,'denseforce',denseforce);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fparams.point_charges=charges;
Fparams.init_charges=charges; 
Fparams.q=q;
%Run Rigid Body Stokes 
RBS_mobility(fname,Fparams,[]);
%potential_cross_section;


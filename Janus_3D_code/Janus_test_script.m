clear;clc;
p=4;
lambda=1;
fname='tet_def_big_boi'
n=2;
Cdst=3;
ep=.3;
Nt=1000;
dt=.1;
tdisc='euler';




boundary_label = @(X,y) (0.5*X*y'./sqrt(sum(X.^2,2)).^2 + 1/2);

addpath ./; 
addpath ./support; 
addpath ./LCPsolvers/Num4LCP_MatLab; 
addpath ./LCPsolvers/Num4LCP_MatLab/ext;
addpath ./FMMLIB/stfmmlib3d-1.2/matlab;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Particle centers (cubic lattice in this example)
lx=0:Cdst:Cdst*(n-1);
lx = lx - mean(lx); 

%[xx yy,zz]= meshgrid(lx,lx,[-2, 2]');

[xx,yy] = meshgrid(lx);
%C = [xx(:) yy(:) zeros(size(xx(:)))]; 
%C = C([1 2 4 6],:);
C=3*[1 0 0; -0.5, 0.833,  0; -0.5 -0.833 0; 0 0 1];
% C = 3*[1 0 0; -1 0 0];
gamma = 1;
%C=[0 0 0];
%C=(2/3)*[-3 -3 0; -3 3 0; 3 -3 0; 3 3 0];
%C=[-3 -3 0; -3 3 0];

%C = C + 0.1*rand(size(C));
display(C); 
n3 = size(C,1); 
rd=[2 2 2 2];
rng('default');
%ID = [0.5*rand(n3,1), rand(n3,1)-0.5, zeros(n3,1)];
ID = - C;
  
 
      




%ID = [0.5*rand(n3,1), rand(n3,1)-0.5, zeros(n3,1)];
%{
init_dir=zeros(n3,3);
for ii=1:n3
    init_dir(ii,:) = Mt{sum(tt>0),ii} * ID(ii,:)'
end
%}

init_dir= -C;
%init_dir=[0 1 0; 0 -1 0];
%init_dir=[repmat([0 0 1],n3/2,1); repmat([0,0,-1],n3/2,1)];
init_dir=init_dir./repmat(sqrt(init_dir(:,1).^2+init_dir(:,2).^2+init_dir(:,3).^2),1,3);
tol=1e-4; mdist=3; denseMV=0; denseforce=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% body parameters
parbd = struct('Shape','','n3',n3,'rd',rd,'p',p,'Ct',C,'mdist',mdist,'eps',ep,'out',1);

% linear solver parameters
parslv = struct('solver','gmres','tol',tol,'maxit',200,'rst',4,'prtype','bkdiag','prec',[],...
    'colsolver','BBPGD','coltol',1e-2,'colmaxit',100,'col_tolrel',tol,'col_tolabs',0.1*tol);  

%Create Fparams struct 
Fparams = struct('parbd',parbd,'parslv',parslv,...
    'Nt',Nt,'dt',dt,'comp',1,'type','JanusAmp','lambda',lambda,'gamma',gamma,...
'denseMV',denseMV,'typeMV','Vsh','tdisc',tdisc,'init_dir',init_dir,'boundary_label',boundary_label,'denseforce',denseforce);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Run Rigid Body Stokes 

RBS_mobility(fname,Fparams,[]);




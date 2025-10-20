function [Fparams]=Test_ModLap_Mobility_bimetallic_charges(fname,n,rd,Cdst,p,ep,Nt,dt,tdisc,lambda,e_north,e_south)
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
addpath ./FMMLIB/fmmlib3d-1.2/matlab/;
addpath ./FMMLIB/stfmmlib3d-1.2/matlab/;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shell parameters 
psh = 8; %higher p for spherical shell (field is propagated through it) 
epsh = 0.1; %collision epsilon (relative to shell radius) 
mdsh = 0.75; %particles closer than (1-mdsh)*shrd are considered near
shrd = 8; %shell radius

%Particle centers (cubic lattice in this example)

rng('default');
Cdst = mean(rd)*Cdst; 
lx=0:Cdst:Cdst*(n-1);
lx = lx - mean(lx); 
[xx,yy,zz] = meshgrid(lx);
C = [xx(:) yy(:) zz(:)];

% particle radius
n3 = size(C,1); 
if size(rd,1) == 1
    rd=rd*ones(n3,1);
end

% initial orientations (random)
init_dir=rand(n3,3);
init_dir=init_dir./repmat(sqrt(init_dir(:,1).^2+init_dir(:,2).^2+init_dir(:,3).^2),1,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate internal point charges and charge strength
resolution = 4;
shrink_factor=0.8;
[Xcharge,Ycharge,Zcharge]=sphere(resolution);
Xcharge=reshape(Xcharge,[],1);
Ycharge=reshape(Ycharge,[],1);
Zcharge=reshape(Zcharge,[],1);

point_charges=zeros(size(Xcharge,1),3,n3);
for i=1:size(C,1)
point_charges(:,:,i)=shrink_factor*[Xcharge, Ycharge, Zcharge];  
end
q=zeros(size(Xcharge,1),1);
for i=1:size(q,1)
    if Xcharge(i)>=0
        q(i)=1;
    else 
        q(i)=-1;
    end
end
q=repmat(0.5*q,1,n3);

% Other misc parameters
tol=1e-4; 
mdist=3; 
denseMV=1;
gamma=1000; 
denseforce=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% body (particle) parameters
parbd = struct('Shape','','n3',n3,'rd',rd,'p',p,'Ct',C,'mdist',mdist,'eps',ep,'out',1);

%shell parameters
% modify shell radius to contain all bodies
shrd1 = shrd; 
NC = sqrt(sum(C.^2,2))+rd;
shrd = max(shrd,max(NC)+2*ep);
if shrd~=shrd1
   fprintf('\n Shell radius too small. New radius is %2.2e',shrd);  
else
   fprintf('\n Shell radius acceptable, shrd = %2.2e',shrd); 
end

% shell parameters
parsh = struct('psh',psh,'shrd',shrd,'mdist',mdsh,'eps',epsh,'out',0);

% External electric field
E0 = rand(3,1); E0=E0./norm(E0); 
Sc = cell(psh,1); 
Sc{psh} = SurfaceSph(shrd*shape_gallery(psh,'')); 
Xs = reshape(Sc{psh}.cart.to_array,[],3);
EdX = Xs*E0; %E(x) = E0^TX
Esh = shAna(EdX); %spharm coeffs for E(x)

% Get Eig for (0.5*I + D_lambda)
DEig = Get_Spectra_Vec(psh,lambda,1,'DMat',0);
DEig=DEig(:); 
sp = (psh+1)^2; 
ii=(1:sp)';
nn=floor(sqrt(ii-1));
DEigv = DEig(nn+1);  

% Solve for sh coeffs of mu
mush = Esh./DEigv; 

% For constant E0, we know mush is non-zero only for n=1. 
mush(1)=0; mush(5:end)=0;
boundarydensity=real(shSyn(mush));
shell_int = RBS_set_params(psh,[0 0 0],shrd,'DL_LMOD_3D','',Sc,1,denseMV,0,mdsh,epsh,0);
shell_int.boundarydensity=boundarydensity; 
shell_int.lambda=lambda; shell_int.Vh=mush; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% linear solver parameters
parslv = struct('solver','gmres','tol',tol,'maxit',200,'rst',4,'prtype','bkdiag','prec',[],...
    'colsolver','BBPGD','coltol',1e-4,'colmaxit',100,'col_tolrel',tol,'col_tolabs',0.1*tol);  

%Create Fparams struct 
Fparams = struct('parbd',parbd,'parsh',parsh,'shell_int',shell_int,'parslv',parslv,...
    'Nt',Nt,'dt',dt,'comp',1,'type','JanusEHD','lambda',lambda,'gamma',gamma,'denseforce',denseforce,...
'denseMV',denseMV,'typeMV','Vsh','tdisc',tdisc,'init_dir',init_dir,'epsilon_north',e_north,'epsilon_south',e_south);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fparams.point_charges=point_charges;
Fparams.init_charges=point_charges; 
Fparams.q=q;
%Run Rigid Body Stokes 
RBS_mobility(fname,Fparams,[]);

end
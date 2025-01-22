function Test_Shell_Mobility(fname,n,rd,Cdst,p,ep,Nt,dt,tdisc,psh,epsh,mdsh,shrd)
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

Shell parameters (if nargin is <=9 there is no shell)
psh   - (int)    spherical harmonic order (shell) 
epsh  - (double) epsilon buffer for shell (collision dist) 
mdsh  - (double) near-eval for relative distance to shell center (<1, ~0.75 recommended)
shrd  - (double) shell radius
%}

%Remove addpaths if compiling in command line (mcc)
addpath ./support/; 
addpath ./LCPsolvers/Num4LCP_MatLab; 
addpath ./LCPsolvers/Num4LCP_MatLab/ext;

% variables are converted from strings (to convert from command line call)
%n = str2num(n); p = str2num(p); ep = str2num(ep); Cdst = str2num(Cdst); 
%Nt = str2num(Nt); dt = str2num(dt); rd = str2num(rd); 

if nargin>9
    shell = true; 
    %psh = str2num(psh); epsh = str2num(epsh); mdsh = str2num(mdsh); shrd = str2num(shrd);
else
    shell=false; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Particle centers (cubic lattice in this example)
lx=0:Cdst:Cdst*(n-1);
lx = lx - mean(lx); 
[xx,yy,zz] = meshgrid(lx);
C = [xx(:) yy(:) zz(:)]; 
n3 = size(C,1); 

%random sphere radii in [rd,rd + rng]
rng = 0; 
rd = rd*ones(n3,1) + rng*rand(n3,1); 
NC = sqrt(sum(C.^2,2))+rd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Force and torque functions (weight proportional to volume)
Tfun = @(t,C) zeros(3,n3); 
Ffun = @(t,C) [zeros(2,n3) ; -10*(rd.').^3];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol=1e-4; mdist=3; 
% body parameters
parbd = struct('Shape','','n3',n3,'rd',rd,'p',p,'Ct',C,'mdist',mdist,'eps',ep,'out',1);

if shell
% set shell radius (should contain all bodies)
shrd1 = shrd; 
shrd = max(shrd,max(NC));
if shrd~=shrd1
   fprintf('\n Shell radius too small. New radius is %d',shrd);  
end

% shell parameters
parsh = struct('psh',psh,'shrd',shrd,'mdist',mdsh,'eps',epsh,'out',0);
end

% linear solver parameters
parslv = struct('solver','gmres','tol',tol,'maxit',200,'rst',4,'prtype','bkdiag','prec',[],...
    'colsolver','BBPGD','coltol',1e-4,'colmaxit',100,'col_tolrel',tol,'col_tolabs',0.1*tol); 

denseMV=1; 

%Create Fparams struct 
if shell
    Fparams = struct('parbd',parbd,'parsh',parsh,'parslv',parslv,...
    'Nt',Nt,'dt',dt,'comp',1,'type','FTfun','denseMV',denseMV,...
    'typeMV','Vsh','tdisc',tdisc,'Tfun',Tfun,'Ffun',Ffun);
else
    Fparams = struct('parbd',parbd,'parslv',parslv,...
    'Nt',Nt,'dt',dt,'comp',1,'type','FTfun','denseMV',denseMV,...
    'typeMV','Vsh','tdisc',tdisc,'Tfun',Tfun,'Ffun',Ffun);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Run Rigid Body Stokes 
RBS_mobility(fname,Fparams,[]);

end
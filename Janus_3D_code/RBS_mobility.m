function RBS_mobility(fname,Fparams,init)
%{
General purpose test code for the mobility problem for Stokesian rigid body 
suspensions 
 
Inputs: 
fname - (string) filename for experiment info

Fparams - (struct) parameter struct for rigid body simulation, with fields:

type     - (string) mobility problem type (FTfun,Purcell, 3sphx, 3sphGolx, MHD)
typeMV   - (string) 'Vsh' for spherical harmonics, 'Rbs' for rotation based singular quad
denseMV  - (bool) dense vs FMM for far-field
Nt       - (int)    number of timesteps
dt       - (double) timestep length
timedisc - (string) time discretization (euler, trapz,rk4)
comp     - (bool) compute intermediate quantities FT aznd VW
lambda   - (double) parameter for modified laplace case
e_north/e_south - (double) Janus particle relative permittivity 
parbd - (struct) struct with rigid body parameters:
    Shape - (string) rigid body shape, '' is default for sphere
    n3    - (int) number of rigid bodies n_b
    rd   - (double) radius (monodisperse) or array of n_b radii (polydisperse) of bodies
    p     - (int)    spherical harmonic order (bodies) 
    Ct    - (double) n_b x 3 array of centers 
    eps   - (double) epsilon buffer (collision dist)
    mdist - (double) collision buffer for body-body interactions
    out   - (bool) external vs internal evaluation (set to 1) 

parslv - (struct) linear solver parameters such as 
     prec     - (string) preconditioner type, '' for unprec, 'bkdiag'
     (block diagonal), 'TT' (tensor train)
     solver   - gmres, pcg, bicg, etc. 
     tol (tolerance), maxit (maximum iterations), rst (restart), etc.
       
parsh - (struct) optional shell geometry parameters: 
    psh - (int)    spherical harmonic order
    shrd - (double) shell radius
    mdsh - (double) min relative distance to origin for near-sing (<1)
    epsh - (double) collision buffer with geometry
    out   - (bool) external vs internal evaluation (set to 0)

Depending on type, extra parameters might be required. 
The default 'FTfun' (force and torque prescription) requires functions 
Ffun,Tfun = @(t,C) with output of size 3 x n_b. 

init    - (string) optional filename to resume a simulation from last
recorded timestep
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global timings;
 
%(0.1) (optional) Load data in init, initialize output arrays
Nt = Fparams.Nt; n3=Fparams.parbd.n3; 
tt = zeros(Nt+1,1); 
sigma = cell(Nt+1,1); mu=sigma; U=sigma; VW=U; Xt=U; Ct=Xt; FT=mu; psi_Lap = mu; 
Mt = cell(Nt+1,n3);
dt0 = Fparams.dt; 

if isempty(init)
    initfl = false; 
    
    % Evolution
    t=0; tt(1)=0;  
    Fparams.parslv.prev=[]; %initialize preconditioner params
 
    for k=1:n3
    Mt{1,k}=eye(3);   
    end
    
    Mt0 = Mt(1,:);
else
   initfl = true;  
    
   % Load previous file 
   load(init);
   
   lid = sum(tt>0);
   % Recover initial data at time t0 = tt(lid): 
   T0 = tt(lid); Xt0 = Xt{lid}; Ct0 = Ct{lid}; Mt0 = Mt(lid,:);
   VW0 = VW{lid}; 
   
   Fparams.parslv.prec=[]; 
   Fparams.parslv.prev=[]; 
  % Fparams.parslv.prtype=prtype;    
   %Re-initialize arrays (after load)
   tt = zeros(Nt+1,1); 
   sigma = cell(Nt+1,1); mu=sigma; U=sigma; VW=U; Xt=U; Ct=Xt; FT=mu;
   Mt = cell(Nt+1,n3);
    
    t=T0; tt(1)=T0; 
    Xt{1}=Xt0; Ct{1}=Ct0; Mt(1,:) = Mt0; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(0.2) Initialize timings and Fparams struct before simulation

zN = zeros(Nt,1); 
timings = struct('setup_surf',0,'setup_kernel',0,'incoming',zN,...
    'velocities',struct('solve',zN,'apply',zN,'vw',zN,'col',zN,'shell',zN,'total',zN),...
    'advance',zN,...
    'operator',struct('surf',zN,'diag',zN,'offd',zN,'total',zN),'total',zN);

tic;
Fparams = RBS_Initialize_params(Fparams);
timings.setup_surf = toc;

timedisc = Fparams.tdisc; Sc = Fparams.parbd.Sc; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%(0.3) Initialize Kernels (MatVecs) and Nullspace info
tic; 
Xrp = Fparams.parbd.Xrp; C0 = Fparams.parbd.C; 
np = Fparams.parbd.np; nrmW = zeros(n3,1);  
if initfl
for k=1:n3
   nrmW(k) = norm(VW0(4:6,k)); 
   xind = (1:np)+np*(k-1);       
   Xrp(xind,:) = Xrp(xind,:)*Mt0{k}';
end 
end

Kernels=[]; 
[Kernels,Nullsp,Fparams,timings] = RBS_Update_Operators(Xrp,C0,Mt0,nrmW,Kernels,Fparams,timings,0); 

timings.setup_kernel=toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (0.4) Initialize collision info  
[colevent,collist,~,~] = LOCAL_check_collision_sph(C0,Fparams);

if ~strcmp(Fparams.parbd.Shape,'') 
   Sc2 = SurfaceSph(rad*shape_gallery(2*p,Fparams.parbd.Shape)); 
   % Model surface pts
   X2 = repmat(reshape(Sc2.cart.to_array,[],3),n3,1);    
   np2 = 2*(2*p)*(2*p+1);    
else  
   X2=[]; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotation matrix (Rodrigues rotation formula)
MRot = @(wh,t) RotationMat(wh,t);

Xt{1}=Xrp; Ct{1}=Fparams.parbd.C;
Fparams.saveLCPs = true; % TODO: change this to false
for i=1:Nt
    if i == Nt 
        Fparams.endFlag = true; 
    else 
        Fparams.endFlag = false; 
    end
dt=dt0; 

if strcmp(timedisc,'euler')
fprintf('\n ---------------------------------------------------------- \n')
fprintf('\n (1) Explicit euler step \n')
[Xt{i+1},Mt(i+1,:),Ct{i+1},U{i},FT{i},sigma{i},mu{i},VW{i},Kernels,Nullsp,...
    Fparams,colevent,collist,dt,psi_Lap{i},Energy(i)] = ...
    LOCAL_euler_step(Xt{i},Xt{1},X2,Mt(i,:),Ct{i},Kernels,Nullsp,Fparams,colevent,collist,t,dt,i);

elseif strcmp(timedisc,'trapz')
% (1) Predictor step: 
fprintf('\n ---------------------------------------------------------- \n')
fprintf('\n (1) Trapezoidal, predictor step \n')
[Xt1,Mt1,Ct1,U1,FT1,sigma1,mu1,VW1,Kernels,Nullsp,...
    Fparams,colevent,collist,dt] = ...
    LOCAL_euler_step(Xt{i},Xt{1},X2,Mt(i,:),Ct{i},Kernels,Nullsp,Fparams,colevent,collist,t,dt,i);

% (2) Corrector step: 
fprintf('\n ---------------------------------------------------------- \n')
fprintf('\n (2) Trapezoidal, corrector step \n')
% Get incoming force distribution: 
tic; 
[FT2,sigma2,VW2,psi_Lap,Energy] = LOCAL_get_incoming_Fc(Fparams,t+dt,dt,Kernels,Nullsp,Xt1,Sc); 
fprintf('\n Time to compute incoming force: %e',toc);
timings.incoming(i) = timings.incoming(i)+toc;  
fprintf('\n Fluid Solve at time %.2f',t)
tic; 
[sigma2,mu2,U2,VW2] = LOCAL_compute_velocities(sigma2,VW2,Ct1,Kernels,Nullsp,Fparams,colevent,collist,i,dt); 
timings.velocities.total(i) = timings.velocities.total(i) + timings.velocities.solve(i) + timings.velocities.apply(i) + timings.velocities.vw(i) + timings.velocities.col(i); 
fprintf('\n Time to compute velocities / fluid solve: %e',timings.velocities.total(i)); 

% Correct VW{i} as average of VW0 and VW1 (and associated quantities)
VW{i}    = 0.5*(VW1+VW2); 
sigma{i} = 0.5*(sigma1+sigma2); 
mu{i}    = 0.5*(mu1+mu2);
U{i}     = 0.5*(U1+U2); 
FT{i}    = 0.5*(FT1+FT2); 

[Xt{i+1},Mt(i+1,:),Ct{i+1},Kernels,Nullsp,Fparams,colevent,collist,dt] = ...
    LOCAL_advance_step(VW{i},mu{i},sigma{i},Xt{i},Xt{1},X2,Mt(i,:),Ct{i},Kernels,Fparams,dt,i);

elseif strcmp(timedisc,'rk4')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Explicit Runge-Kutta 4th order method
% (1) First step, f1 = f(t0,u0)
fprintf('\n ---------------------------------------------------------- \n')
fprintf('\n (1) Runge-Kutta 4th order, t=t0, x=x(t0) \n')
dt41 = 0.5*dt; 
t41 = t;

[Xt41,Mt41,Ct41,U41,FT41,sigma41,mu41,VW41,Ker41,Null41,...
    Fpar41,cev41,clst41,dt41] = ...
    LOCAL_euler_step(Xt{i},Xt{1},X2,Mt(i,:),Ct{i},Kernels,Nullsp,Fparams,colevent,collist,t41,dt41,i);

% (2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n ---------------------------------------------------------- \n')
fprintf('\n (2) Runge-Kutta 4th order, t_{1/2}=t0+(dt/2), x=x(t_0)+(dt/2)*v1 \n')
t42 = t+dt41; 
dt42 = 0.5*dt; 

% Get incoming force distribution: 
tic; 
[FT42,sigma42,VW42,psi_Lap, Energy] = LOCAL_get_incoming_Fc(Fpar41,t42,dt42,Ker41,Null41,Xt41,Sc); 
fprintf('\n Time to compute incoming force: %e',toc);
timings.incoming(i) = timings.incoming(i)+toc;  
fprintf('\n Fluid Solve at time %.2f',t)
tic; 
[sigma42,mu42,U42,VW42] = LOCAL_compute_velocities(sigma42,VW42,Ct41,Ker41,Null41,Fpar41,cev41,clst41,i,dt42); 
timings.velocities.total(i) = timings.velocities.total(i) + timings.velocities.solve(i) + timings.velocities.apply(i) + timings.velocities.vw(i) + timings.velocities.col(i); 
fprintf('\n Time to compute velocities / fluid solve: %e',timings.velocities.total(i));

% Advance Xt42 ~ X(t) + (dt/2)*V42 
[Xt42,Mt42,Ct42,Ker42,Null42,Fpar42,cev42,clst42,dt42] = ...
    LOCAL_advance_step(VW42,mu42,sigma42,Xt{i},Xt{1},X2,Mt(i,:),Ct{i},Kernels,Fparams,dt42,i);

% (3) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n ---------------------------------------------------------- \n')
fprintf('\n (3) Runge-Kutta 4th order, t_{1/2}=t0+0.5*dt, x=x(t_{1/2}) \n')
t43 = t+dt42; 
dt43 = dt; 

% Get incoming force distribution: 
tic; 
[FT43,sigma43,VW43] = LOCAL_get_incoming_Fc(Fpar42,t43,dt43,Ker42,Null42,Xt42,Sc); 
fprintf('\n Time to compute incoming force: %e',toc);
timings.incoming(i) = timings.incoming(i)+toc;  
fprintf('\n Fluid Solve at time %.2f',t)
tic; 
[sigma43,mu43,U43,VW43] = LOCAL_compute_velocities(sigma43,VW43,Ct42,Ker42,Null42,Fpar42,cev42,clst42,i,dt43); 
timings.velocities.total(i) = timings.velocities.total(i) + timings.velocities.solve(i) + timings.velocities.apply(i) + timings.velocities.vw(i) + timings.velocities.col(i); 
fprintf('\n Time to compute velocities / fluid solve: %e',timings.velocities.total(i));

% Advance Xt43 ~ X(t) + (dt)*V43 
[Xt43,Mt43,Ct43,Ker43,Null43,Fpar43,cev43,clst43,dt43] = ...
    LOCAL_advance_step(VW43,mu43,sigma43,Xt{i},Xt{1},X2,Mt(i,:),Ct{i},Kernels,Fparams,dt43,i);

% (4) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n ---------------------------------------------------------- \n')
fprintf('\n (4) Runge-Kutta 4th order, t_{1}=t0+dt, x=x(t_{1}) \n')
t44 = t+dt43; 
dt44 = dt; 
% Get incoming force distribution: 
tic; 
[FT44,sigma44,VW44] = LOCAL_get_incoming_Fc(Fpar43,t44,dt44,Ker43,Null43,Xt43,Sc); 
timings.incoming(i) = timings.incoming(i)+toc;
fprintf('\n Fluid Solve at time %.2f \n',t)
[sigma44,mu44,U44,VW44] = LOCAL_compute_velocities(sigma44,VW44,Ct43,Ker43,Null43,Fpar43,cev43,clst43,i,dt44); 

timings.velocities.total(i) = timings.velocities.total(i) + timings.velocities.solve(i) + timings.velocities.apply(i) + timings.velocities.vw(i) + timings.velocities.col(i); 
fprintf('\n Time to compute velocities / fluid solve %.2f \n',timings.velocities.total(i));

% Average velocities and related quantities using RK4 weights
VW{i}    = (1/6)*(VW41+2*VW42+2*VW43+VW44); 
sigma{i} = (1/6)*(sigma41+2*sigma42+2*sigma43+sigma44); 
mu{i}    = (1/6)*(mu41+2*mu42+2*mu43+mu44); 
U{i}     = (1/6)*(U41+2*U42+2*U43+U44); 
FT{i}    = (1/6)*(FT41+2*FT42+2*FT43+FT44);

% Advance Xtp ~ X(t) + (dt/6)*(V41 + 2*V42 + 2*V43 + V44) 
[Xt{i+1},Mt(i+1,:),Ct{i+1},Kernels,Nullsp,Fparams,colevent,collist,dt] = ...
    LOCAL_advance_step(VW{i},mu{i},sigma{i},Xt{i},Xt{1},X2,Mt(i,:),Ct{i},Kernels,Fparams,dt,i);

end

timings.total(i) = timings.velocities.total(i) + timings.operator.total(i) + ...
    timings.advance(i) + timings.incoming(i); 
fprintf('\n Total computing time for timestep %d : %e ',i,timings.total(i))
fprintf('\n -------------------------------------------------------------\n'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = t+dt;
tt(i+1)=t;
fprintf('\n dt: %2.2f ',dt)
if mod(i,2)==1                                                                                                                                              
save(fname,'-v7.3','tt','Xt','Mt','Ct','FT','sigma','mu','U','VW','psi_Lap','Energy');                                                                                         
else                                                                                                                                                        
save([fname '2'],'-v7.3','tt','Xt','Mt','Ct','FT','sigma','mu','U','VW','psi_Lap','Energy');                                                                                   
end  
save([fname '_profile'],'timings'); 
end


end

function [Xtp,Mtp,Ctp,U,FT,sigma,mu,VW,Kernels,Nullsp,Fparams,colevent,collist,dt,psi_Lap,Energy] = LOCAL_euler_step(Xt,X0,X2,Mt,Ct,Kernels,Nullsp,Fparams,colevent,collist,t,dt,it)
global timings; 
Sc = Fparams.parbd.Sc; 
diam = Fparams.parbd.diam; 
mxrd = Fparams.parbd.mxrd; 
Shape = Fparams.parbd.Shape; 
np = Fparams.parbd.np; 
n3 = Fparams.parbd.n3; 

MRot = @(wh,t) RotationMat(wh,t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Get incoming force distribution: 
tic; 
[FT,sigma,VW,psi_Lap,Energy] = LOCAL_get_incoming_Fc(Fparams,t,dt,Kernels,Nullsp,Xt,Sc); 
Ct
fprintf('\n Time to compute incoming force: %e ',toc)
timings.incoming(it) = timings.incoming(it) + toc; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n Fluid Solve at time %.2f ',t)
tic; 
[sigma,mu,U,VW] = LOCAL_compute_velocities(sigma,VW,Ct,Kernels,Nullsp,Fparams,colevent,collist,it,dt); 
timings.velocities.total(it) = timings.velocities.total(it) + timings.velocities.solve(it) + timings.velocities.apply(it) ...
    + timings.velocities.vw(it) + timings.velocities.col(it); 
fprintf('\n Time to compute velocities / fluid solve: %e',timings.velocities.total(it)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Advance center Ct 
tic; 
Ctp = LOCAL_advance_center(Ct,dt,VW,Fparams); 
fprintf('\n Time to advance centers C(t): %e',toc); 
timings.advance(it) = timings.advance(it) + toc; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for collision after moving centers
tic; 
[colevent,collist,dt,Ctp] ...
= LOCAL_collision_info(Fparams,X2,Ct,Ctp,VW,MRot,Mt,dt);
fprintf('\n Time for collision detection: %e',toc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Advance rotation matrix Mt and X
tic; 
[Mtp,Xtp,nrmW] = LOCAL_advance_rotation(MRot,VW,Mt,Xt,X0,dt,np,n3); 
fprintf('\n Time to advance R(t) and X(t): %e',toc)
timings.advance(it) = timings.advance(it) + toc; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Update Sc and operators
fprintf('\n Surface and operator update')
if isfield(Fparams,'parsh')
tic; 
Fparams = LOCAL_compute_shell_velocity(mu+sigma,Fparams);
fprintf('\n Time to compute boundary correction: %e',toc)
timings.velocities.shell(it) = toc;
timings.velocities.total(it) = timings.velocities.total(it) + timings.velocities.shell(it); 
end

[Kernels,Nullsp,Fparams,timings] ...
   = RBS_Update_Operators(Xtp,Ctp,Mtp,nrmW,Kernels,Fparams,timings,it);
timings.operator.total(it) = timings.operator.total(it) + timings.operator.surf(it) + timings.operator.diag(it) + timings.operator.offd(it);
fprintf('\n Time to update surface and operators: %e',timings.operator.total(it));  

end

function [Xtp,Mtp,Ctp,Kernels,Nullsp,Fparams,colevent,collist,dt] = LOCAL_advance_step(VW,mu,sigma,Xt,X0,X2,Mt,Ct,Kernels,Fparams,dt,it)
global timings;  
 
np = Fparams.parbd.np; 
n3 = Fparams.parbd.n3; 

MRot = @(wh,t) RotationMat(wh,t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Advance center Ct 
tic; 
Ctp = LOCAL_advance_center(Ct,dt,VW,Fparams); 
fprintf('\n Time to advance centers C(t): %e',toc); 
timings.advance(it) = timings.advance(it) + toc; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for collision after moving centers
tic; 
[colevent,collist,dt,Ctp] ...
= LOCAL_collision_info(Fparams,X2,Ct,Ctp,VW,MRot,Mt,dt);
fprintf('\n Time for collision detection: %e',toc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Advance rotation matrix Mt and X
tic; 
[Mtp,Xtp,nrmW] = LOCAL_advance_rotation(MRot,VW,Mt,Xt,X0,dt,np,n3); 
fprintf('\n Time to advance R(t) and X(t): %e',toc)
timings.advance(it) = timings.advance(it) + toc; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Update Sc and operators
fprintf('\n Surface and operator update')
if isfield(Fparams,'parsh')
tic; 
Fparams = LOCAL_compute_shell_velocity(mu+sigma,Fparams);
fprintf('\n Time to compute boundary correction: %e',toc)
timings.velocities.shell(it) = toc;
timings.velocities.total(it) = timings.velocities.total(it) + timings.velocities.shell(it); 
end
 
[Kernels,Nullsp,Fparams,timings] ...
   = RBS_Update_Operators(Xtp,Ctp,Mtp,nrmW,Kernels,Fparams,timings,it);
timings.operator.total(it) = timings.operator.total(it) + timings.operator.surf(it) + timings.operator.diag(it) + timings.operator.offd(it);
fprintf('\n Time to update surface and operators: %e',timings.operator.total(it));

end

function [FT,fM,VW,psi_Lap,Energy] = LOCAL_get_incoming_Fc(Fparams,t,dt,Kernels,Nullsp,Xt,Sc)
%global acc rst maxit lprec; 
Energy = 0;
parslv = Fparams.parslv; 
lprec=[]; 
acc = parslv.tol; rst = parslv.rst; maxit = parslv.maxit;

%W = Fparams.parbd.W; 
n3=Fparams.parbd.n3; p = Fparams.parbd.p; np = Fparams.parbd.np; 
VW=[]; C = Fparams.parbd.C; 
Bk = Nullsp.B; Ck = Nullsp.C; Lk = Nullsp.L; 

switch Fparams.type
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Modified Laplace, Electrostatic case (Janus particles) 
    case 'JanusEHD'
        %obtain kernels for S,D,S',D' in Laplace, Modified Laplace case
        q=Fparams.q;
        SLD=Kernels.SLD; SLMODD=Kernels.SLMODD;
        DLD=Kernels.DLD; DLMODD=Kernels.DLMODD;
        dSLD=Kernels.dSLD;dSLMODD=Kernels.dSLMODD;
        dDLD=Kernels.dDLD;dDLMODD=Kernels.dDLMODD;

        Epsilon=Fparams.Epsilon;
        EpsilonMat=repmat(Epsilon,1,size(Epsilon,1));
       
        if isa(SLD,'function_handle')
           K1 = @(V) DLMODD(V) - DLD(V) + V;
           K2 = @(V) SLMODD(V) - SLD(V);
           K3 = @(V) dDLMODD(V) - Epsilon.*(dDLD(V));
           K4 = @(V) dSLMODD(V) - Epsilon.*dSLD(V) - ((1+Epsilon)/2).*V;
            
            
           M = @(V)  [ K1(V(1:end/2)) + K2(V(end/2+1:end)) ; K3(V(1:end/2)) + K4(V(end/2+1:end))];
        
        else
           K1=DLMODD-DLD+eye(size(DLD));
           K2=SLMODD-SLD;
           K3=dDLMODD-EpsilonMat.*dDLD;
           K4=dSLMODD-EpsilonMat.*dSLD - diag((1+Epsilon)/2);
           M=[K1 K2; K3 K4];
        
        
        end

        % compute rhs
        params=Fparams.parmod;
        charges=Fparams.point_charges;
        charges=permute(charges,[1 3 2]);
        charges=reshape(charges,[],size(Fparams.point_charges,2),1);
        res=size(Fparams.point_charges,1);
        %make centers better
        c1=reshape(repmat(C(:,1)',res,1),[],1);
        c2=reshape(repmat(C(:,2)',res,1),[],1);
        c3=reshape(repmat(C(:,3)',res,1),[],1);
        Cup=[c1 c2 c3];
        
        
        params.W2=reshape(reshape(q,[],1),[],1)';
        params.flag_pot='SL_L_3D';
        Q1=sum(Kernel_Eval(params.X,charges + Cup,params),2);
        Q2=zeros(size(Q1));
        params.flag_pot='dSL_L_3D';
       % Q2=sum(Kernel_Eval(params.X,charges + Cup,params),2);
        
        %Electric Field Condition
        Epotential=zeros(size(Q1));

   %     Epotential = params.X*Fparams.Efield';
   %     dEpotential = params.nor*Fparams.Efield';

        shell_int = Fparams.shell_int;
        % Compute D[mu_inf] at target points
        Epotential = VSh_Mod_MatVec_RB_trg(shell_int.boundarydensity,params.X,params.nor,shell_int);
        %shell_int.flag_pot='SL_LMOD_3D';
        %Epotential = Epotential + VSh_Mod_MatVec_RB_trg(shell_int.boundarydensity,params.X,params.nor,shell_int);
        dEpotential = Q2;
        
      for j=1:n3
            indx=(1:np)+np*(j-1); 
            %indv=(1:3*np)+3*np*(j-1); 
            sj = min(j,size(Sc,2));   
           % H_i{j} = Pot2Field(phi(indx), phi_n_i(indx),Sc{p,sj});
           S = Sc{p,sj};
           G = S.geoProp.Grad(Q1(indx));
           Q2(indx)=G.dot(vec3d(params.nor(indx,:)));
           Epot= S.geoProp.Grad(Epotential(indx));
           dEpotential(indx) = Epot.dot(vec3d(params.nor(indx,:)));
      end
%}      
           
      
        
   
 %       Q2=sum(Kernel_Eval(params.X,charges + Cup,params),2);
        Q=[(Q1 - Epotential);Epsilon.*(Q2 - dEpotential)];
   %              Q=[(Q1);Epsilon.*(Q2)];
  %solves for potential
        [mupsi,~,~,I]=gmres(M, Q,100,1e-6);
        display('gmres iters for force:');
        display(I);
        Kond = cond(M);
        mu=mupsi(1:end/2); psi=mupsi(end/2+1:end);
        
        if isa(SLD,'function_handle')
        
        phi=SLD(psi)+DLD(mu)+Q1;% + Epotential;
        phi_n_e=-0.5*psi+ dSLMODD(psi) + dDLMODD*(mu);
        phi_n_i=0.5*psi + dSLD(psi) + dDLD*(mu) + Q2;
        
        else
         phi=SLD*(psi)+DLD*(mu)+Q1;% + Epotential;
        phi_n_e=-0.5*psi+ dSLMODD*(psi) + dDLMODD*(mu);
        phi_n_i=0.5*psi + dSLD*(psi) + dDLD*(mu) + Q2;
        end
        Pot2Field = @(phi, phi_n,S) -1*S.geoProp.Grad(phi) -1*vec3d([phi_n; phi_n; phi_n]).*S.geoProp.nor;  %  -Grad phi - phi_n n
        maxwellSnor = @(E,S) times(dot(E, S.geoProp.nor), E) - times(dot(E,E), S.geoProp.nor)/2; % n \cdot (E \oprod E - 1/2 |E|^2 I)

        
        fM = zeros(3*np*n3,1); 
        H_i = cell(n3,1); H_e=H_i; 
        for j=1:n3
            indx=(1:np)+np*(j-1); 
            indv=(1:3*np)+3*np*(j-1); 
            sj = min(j,size(Sc,2));   
            H_i{j} = Pot2Field(phi(indx), phi_n_i(indx),Sc{p,sj});
            H_e{j} = Pot2Field(phi(indx), phi_n_e(indx),Sc{p,sj});
            %Maxwell stress . normal (traction)
            ftmp = maxwellSnor(H_e{j},Sc{p,sj}) - maxwellSnor(H_i{j},Sc{p,sj}); 
            ftmp = real(reshape(ftmp.to_array,[],3))'; 
            fM(indv) = ftmp(:); 
        end
        psi_Lap=[mupsi; I(1,2) ; Kond];
        fprintf('\n Electrostatic forces and torques on Janus particles \n'); 
        FT = real(Ck*fM); 
        display(reshape(FT,6,n3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    case 'JanusAmp'
        SLMODD=Kernels.SLMODD; dSLMODD=Kernels.dSLMODD;
        DLMODD=Kernels.DLMODD; dDLMODD=Kernels.dDLMODD;
        flabel=Fparams.SurfaceLabel;
        
        %solves for density, psi & uses them to compute normal derivative
        if isa(SLMODD,'function_handle')
        K = @(V) SLMODD(V) + DLMODD(V);
        P = make_projection(Fparams.parmod.p,1);
        Proj = @(V) reshape(P*reshape(V,[],n3),[],1);
        K_full_rank = @(V) Proj(K(Proj(V))) + V - Proj(V);    
        [psi,~,~,I]=gmres(K_full_rank, Proj(flabel),100,1e-6);
        fprintf('\n Janus Amph S+D BIE solve error = %e, %e',norm(K(psi)-Proj(flabel))/norm(flabel), norm(K_full_rank(psi)-Proj(flabel))/norm(flabel)); 
        
        
        phi = K(psi); 
        phi_n_e=dSLMODD(psi) + dDLMODD(psi);
        
        
        else
        K=SLMODD+DLMODD;
        P=make_projection(Fparams.parmod.p,n3);
        K_full_rank=P*K*P+(eye(np*n3)-P);
        [psi,~,~,I]=gmres(K_full_rank, P*(flabel),100,1e-6);
        %display(cond(K_full_rank));
         Kond = cond(K_full_rank); 

        fprintf('\n Janus Amph S+D BIE solve error = %e, %e',norm(K*psi-P*flabel)/norm(flabel), norm(K_full_rank*psi-P*flabel)/norm(flabel)); 
        phi = K*psi; 
        phi_n_e=dSLMODD*psi + dDLMODD*psi;
        
        end
    
        
        %computes energy
        W = repmat(Fparams.parmod.W,n3,1);
        Energy = real(- W' * (phi.*phi_n_e));
        
        
        phisquared=phi.*phi;
        gradphi_mag_squared=zeros(n3*np,1);
%       phi_n_e_2=gradphi_mag_squared;
%        gradphi = zeros(n3*np,3);
        %equations
        Pot2Field = @(phi, phi_n,S) S.geoProp.Grad(phi) + vec3d([phi_n; phi_n; phi_n]).*S.geoProp.nor;  %  Grad phi + phi_n n
        % n \cdot (E \oprod E - 1/2 |E|^2 I) - sqrt(lambda)*u^2
      %  maxwellSnor = @(E,S,phi2) -(2/sqrt(Fparams.lambda))*(times(dot(E, S.geoProp.nor), E) - times(dot(E,E), S.geoProp.nor)/2) + sqrt(Fparams.lambda)*times(phi2,S.geoProp.nor);  
         maxwellSnor = @(E,S,phi2)  (Fparams.gamma) *((-1)*(times(dot(E, S.geoProp.nor), E) - times(dot(E,E), S.geoProp.nor)/2) +     (Fparams.lambda)*times(phi2,S.geoProp.nor));
        fM = zeros(3*np*n3,1); 
        H_i = cell(n3,1); H_e=H_i; 
 
        for j=1:n3
            indx=(1:np)+np*(j-1); 
            indv=(1:3*np)+3*np*(j-1); 
            sj = min(j,size(Sc,2));   

           S = Sc{p,sj};
           G = S.geoProp.Grad(phi(indx));


           gradphi_mag_squared(indx) = G.dot(G);
           Weight= Fparams.parmod.W;

           H_e{j} = Pot2Field(phi(indx), phi_n_e(indx),Sc{p,sj});
            %Maxwell stress . normal (traction)
            ftmp = maxwellSnor(H_e{j},Sc{p,sj},phisquared(indx)); 
            ftmp = real(reshape(ftmp.to_array,[],3))'; 
            fM(indv) = ftmp(:); 
        end

        fprintf('\n Forces and Torques on Amphiphillic Janus particles \n'); 
        FT = real(Ck*fM); 
        psi_Lap=[psi; phi_n_e; I(1,2) ; Kond];
        display(reshape(FT,6,n3))
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    % Magnetic potential solve
    case 'MHD'
        KLD = Kernels.KLD; 
        SLD = Kernels.SLD; 
        fprintf('\n Magnetic potential Solve at time %.2f \n',t)
        % Magnetic solve
        % Build rhs
        %np = 2*p*(p+1); 
        rhs = zeros(np*n3,1); 
        if size(Sc,2)>1
            for j=1:n3
                Nr = reshape(Sc{p,j}.geoProp.nor.to_array,[],3);
                indx=(1:np)+np*(j-1); 
                rhs(indx) = Fparams.eta*Nr*Fparams.H0; 
            end
        else
            Nr = reshape(Sc{p}.geoProp.nor.to_array,[],3); 
            rhs = Fparams.eta*repmat(Nr,n3,1)*Fparams.H0; 
        end
        
        % Solve
        qM = Lslv(KLD,rhs,parslv); 

        fprintf('\n Magnetic potential solve res = %1.4g \n',norm(Lapp(KLD,qM)-rhs)); 

        % Compute Maxwell stress, forces and torques
        % dphi/dn at Gamma
        phi_n_e = Fparams.mur/(1-Fparams.mur)*qM; 
        phi_n_i = 1/(1-Fparams.mur)*qM; 
        % phi at Gamma (Continuous)
        phi = -Xt*Fparams.H0 + Lapp(SLD,qM); 

        % Formulas 
        Pot2Field = @(phi, phi_n,S) -1*S.geoProp.Grad(phi) -1*vec3d([phi_n; phi_n; phi_n]).*S.geoProp.nor;  %  -Grad phi - phi_n n
        maxwellSnor = @(E,S) times(dot(E, S.geoProp.nor), E) - times(dot(E,E), S.geoProp.nor)/2; % n \cdot (E \oprod E - 1/2 |E|^2 I)

        % Magnetic Field (exterior and interior)
        fM = zeros(3*np*n3,1); 
        H_i = cell(n3,1); H_e=H_i; 
        for j=1:n3
            indx=(1:np)+np*(j-1); 
            indv=(1:3*np)+3*np*(j-1); 
            sj = min(j,size(Sc,2));   
            H_i{j} = Pot2Field(phi(indx), phi_n_i(indx),Sc{p,sj});
            H_e{j} = Pot2Field(phi(indx), phi_n_e(indx),Sc{p,sj});
            %Maxwell stress . normal (traction)
            ftmp = maxwellSnor(H_e{j},Sc{p,sj}) - maxwellSnor(H_i{j},Sc{p,sj}); 
            ftmp = real(reshape(ftmp.to_array,[],3))'; 
            fM(indv) = ftmp(:); 
        end

        fprintf('\n Magnetic forces and torques \n'); 
        FT = real(Ck*fM); 
        display(reshape(FT,6,n3))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mobility matrix method for 3-bead swimmer (x-axis)
    case '3Sphx'
TD = Kernels.TD; SD = Kernels.SD; 
vxind = 6*(0:n3-1)+1;
Bx = Bk(vxind,:)';
MM = -Lapp(TD,Bx) + Lk*Bx;
MM = (1/sum(W))*Ck(vxind,:)*(Lapp(SD,(Bx+Lslv(TD,MM,parslv)))); 
F1 = MM\ones(3,1); 

omega=Fparams.omega; ds=Fparams.d; phs=Fparams.phi;  
dL1 = -ds*omega*sin(omega*t+phs); 
dL2 = -ds*omega*sin(omega*t);     
dLV = [0;dL1;dL1+dL2]; 

FLV = MM\dLV; 

v1 = -sum(FLV)/sum(F1); 
Fx = v1*F1+FLV; 
FT = zeros(6*n3,1); 
FT(vxind)=Fx; 
display([v1;v1+dL1;v1+dL1+dL2]')

fM = Bk'*FT;  

VW = zeros(6,n3); 
VW(1,:) = [v1;v1+dL1;v1+dL1+dL2]'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mobility matrix method for 3-bead swimmer (x-axis). Original model by
% Golestanian, with prescribed (linear) strokes.
    case '3SphGolx'
TD = Kernels.TD; SD = Kernels.SD;       
vxind = 6*(0:n3-1)+1;
Bx = Bk(vxind,:)';
MM = -Lapp(TD,Bx) + Lk*Bx;
MM = (1/sum(W))*Ck(vxind,:)*(Lapp(SD,(Bx+Lslv(TD,MM,parslv)))); 
F1 = MM\ones(3,1); 

del=Fparams.d; 
L=Fparams.l; 
%omega = Fparams.omega; 
del = del/L; 
m4 = @(t) mod(floor((t+0.5*dt)/del),4);
m2 = @(t) mod(floor((t+0.5*dt)/del),2);
%meps1 = @(t) m4(t)==1 | m4(t)==2;
%meps2 = @(t) m4(t)==2 | m4(t)==3;
%L1 = @(t) L(1-eps*meps1(t)+(m4(t)-1).*(1-m2(t)).*(t-eps*floor(t/eps)));
%L2 = @(t) L(1-eps*meps2(t)+L*(m4(t)-2).*m2(t).*(t-eps*floor(t/eps)));

dL1 = L*(m4(t)-1).*(1-m2(t)); 
dL2 = L*(m4(t)-2).*m2(t);      
dLV = [0;dL1;dL1+dL2]; 

FLV = MM\dLV; 

v1 = -sum(FLV)/sum(F1); 
Fx = v1*F1+FLV; 
FT = zeros(6*n3,1); 
FT(vxind)=Fx; 
display([v1;v1+dL1;v1+dL1+dL2]')

fM = Bk'*FT;  

VW = zeros(6,n3); 
VW(1,:) = [v1;v1+dL1;v1+dL1+dL2]'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mobility matrix method for Purcell rotator (xy-axis).
    case 'Purcellxy'
TD = Kernels.TD; SD = Kernels.SD;     

%Parameters
tau = Fparams.tau; omega=Fparams.omega; ds=Fparams.d; phs=Fparams.phi; 
% Current position given by 3 angles theta_i(t) and circle radius R
tht = Fparams.tht; R = Fparams.R; 
dB1 = -ds*omega*sin(omega*t+phs); 
dB2 = -ds*omega*sin(omega*t); 
Ret = @(th) R*[-sin(th);cos(th)]; 

% Indices
xyind = reshape(repmat(6*(0:n3-1),3,1),1,[])+repmat([1 2 6],1,n3);
% velocities in xy plane, wz (rotation in xy plane)
%vxyind = reshape(repmat(6*(0:n3-1),2,1),1,[])+repmat([1 2],1,n3);
%wzind = reshape(6*(0:n3-1),1,[])+repmat(6,1,n3);

% Build Mobility matrix (9 x 9)
Bx = Bk(xyind,:)';
MM = -Lapp(TD,Bx)+Lk*Bx;
MM = Ck(xyind,:)*Lapp(SD,(Bx+Lslv(TD,MM,parslv))); 
vid = [1 2 4 5 7 8]; wid = [3 6 9];
MM(vid,:) = (1/sum(W))*MM(vid,:); 
MM(wid,:) = (1/tau(3,3)).*MM(wid,:);    

%V = dth1*d1+dBV = MF
d1  = [Ret(tht(1));1;Ret(tht(2));1;Ret(tht(3));1]; 
dBV = [zeros(3,1); dB1*[Ret(tht(2));1] ; (dB1+dB2)*[Ret(tht(3));1] ]; 

% dth1*F1 + dFV = F
F1 = MM\d1; dFV = MM\dBV; 
dth1 = -sum(dFV(wid))/sum(F1(wid)); 

% Compute forces/torques and density
Fxy = dth1*F1 + dFV; 
FT = zeros(6*n3,1); 
FT(xyind)=Fxy;  
fM = Bk'*FT; 

% display angular and translational velocities
fprintf('\n Angular velocities:\n'); 
wz = [dth1;dth1+dB1;dth1+dB1+dB2]'; 
display(wz) 
fprintf('\n Translational velocities:\n'); 
vxy = R*[dth1*Ret(tht(1)) (dth1+dB1)*Ret(tht(2)) (dth1+dB1+dB2)*Ret(tht(3))]; 
display(vxy) 

VW = zeros(6,n3); 
VW(1:2,:) = vxy; VW(6,:) = wz;  

   case 'FTfun' 
   Ffun = Fparams.Ffun; 
   Tfun = Fparams.Tfun; 
   Ct = Fparams.parbd.C; 
   
   Force = Ffun(t,Ct);  
   Torque = Tfun(t,Ct);
   FT = [Force;Torque]; FT = FT(:); 
   fM = Bk'*FT;   
   
   fprintf('\n Prescribed forces and torques \n'); 
   %display(reshape(FT,6,n3))
end

end

function y = Lapp(A,x)
if isnumeric(A)
    y=A*x; 
else
    y=real(A(x)); 
end
end

function x = Lslv(A,b,parslv)
%global prec;
prec = parslv.prec; 

if nargin<6
   if ~isempty(prec) 
       pr=prec;  
   else
       pr=[]; 
   end
end

if isnumeric(A)
    x=A\b; 
else
    x = zeros(size(b)); 
    %b = parslv.Pr(b); 
    %A = @(x) parslv.Pr(A(x)) + 0.5*(x-parslv.Pr(x)); 
    
    for i=1:size(b,2)
    if nargin==2
        [x(:,i),~,rs,it]=gmres(A,b(:,i),1,1e-6,200,pr);
    else
        [x(:,i),~,rs,it]=gmres(A,b(:,i),parslv.rst,parslv.tol,parslv.maxit,pr);
    end
       fprintf('\n gmres %d iters=%d, res=%1.4g \n',i,prod(it),rs); 
    end
end
end

function den = LOCAL_CenterDistance(C)

[Y_g1,  X_g1  ] = meshgrid(C(:,1), C(:,1));
[Y_g2,  X_g2  ] = meshgrid(C(:,2), C(:,2));
[Y_g3,  X_g3  ] = meshgrid(C(:,3), C(:,3));
d1 = (X_g1 - Y_g1); d2 = (X_g2 - Y_g2); d3 = (X_g3 - Y_g3); 
den = sqrt(d1.^2 + d2.^2 + d3.^2);
den = 10000*(den==0)+den; 

end

function den = LOCAL_Distance(X,Y)

[Y_g1,  X_g1  ] = meshgrid(Y(:,1), X(:,1));
[Y_g2,  X_g2  ] = meshgrid(Y(:,2), X(:,2));
[Y_g3,  X_g3  ] = meshgrid(Y(:,3), X(:,3));
d1 = (X_g1 - Y_g1); d2 = (X_g2 - Y_g2); d3 = (X_g3 - Y_g3); 
den = sqrt(d1.^2 + d2.^2 + d3.^2);
den = 10000*(den==0)+den; 

end

function [f,g] = objGrad(x,A,b)

Ax = A(x);
f = 1/2 * dot(x, Ax) + dot(x, b);
if nargout == 1
    g = [];
    return 
end
g = Ax + b;
end

function [F_c,mu_c,rho_c] = LOCAL_Compute_Contact_LCP(collist,Kernels,Nullsp,Fparams,Ct,VW,dt)

persistent A_list b_list save_iter

if ~isfield(Fparams, 'saveLCPs') 
    saveLCPs = false;
else 
    saveLCPs = Fparams.saveLCPs;
end
if ~isfield(Fparams, 'endFlag') 
    endFlag = false; 
else 
    endFlag = Fparams.endFlag;
end
if saveLCPs && (isempty(A_list) || isempty(b_list) || isempty(save_iter))
    A_list = {};
    b_list = {};
    save_iter = 1;
end


parslv = Fparams.parslv; 
rd = Fparams.parbd.rd; 
n3 = length(rd); 
diam = Fparams.parbd.diam; %diam(i,j) = r_i + r_j
mxrd = Fparams.parbd.mxrd; %max(r_i,r_j)
eps = Fparams.parbd.eps;
Nb = Fparams.parbd.Nb; 
tol = parslv.tol; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Setup (build A and b)
 
shflg = isfield(Fparams,'parsh');
ip = collist(:,1); jp = collist(:,2); 
if shflg
    sheps = Fparams.parsh.eps;
    idsh = jp > n3; 
    ipsh = ip(idsh); 
    ip = ip(~idsh); jp = jp(~idsh); 
    numFS = length(ipsh); 
else
    numFS = 0; 
end

TD = Kernels.TD; SD = Kernels.SD; 
Bk = Nullsp.B; Ck = Nullsp.C; Lk = Nullsp.L; 
numF = length(ip); 

% Compute vectors and normal vectors for pairs
R = Ct(ip,:)-Ct(jp,:);       %Ci - Cj numF x 3
NR = sqrt(sum(R.*R,2));      %|Ci-Cj| numF x 1 
Rhat = repmat(1./NR,1,3).*R; %eij = (Ci - Cj)/|Ci-Cj|

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Build A 

F = zeros(6*n3,numF+numFS); 

for k=1:numF
    indi = (1:3)+6*(ip(k)-1);
    indj = (1:3)+6*(jp(k)-1);
    F(indi,k) = Rhat(k,:); 
    F(indj,k) = -Rhat(k,:); 
end

for k=numF+1:numF+numFS
   indi = (1:3)+6*(ipsh(k-numF)-1);
   F(indi,k) = -Ct(ipsh(k-numF),:)./norm(Ct(ipsh(k-numF),:));
end

% A is built using the dense or matfree mobility matrix. Can be accelerated
% by employing only self interaction (block-diagonal) for TD and SD. 

% Different criteria can be added as needed
% TODO nic change this back 
denseMV = Fparams.denseMV; % || numF+numFS > 1; 
bkdiag=false; 

if denseMV
    parslv.tol = tol; 
    Bf = (Bk.')*F; 
    %(3) (-0.5I-K)*rho_c
    MNS = -Lapp(TD,Bf)+Lk*Bf;
    %(4) 3x3 MNS=VNS*S(mu_c+rho_c)  
    MuNS = Lslv(TD,MNS,parslv); 

    % Setup LCP x perp A*x + b (dense build of Amat = F^T M F)
    Amat = real(F.'*(Ck*Lapp(SD,MuNS+Bf))); 
    A = @(x) Amat*x;
else % matfree
    parslv.tol = parslv.coltol; 
    Bf = @(x) (Bk.')*(F*x);
    if bkdiag % this is just preconditioner on the solve 
        S0 = @(x) reshape(Kernels.SSD0*(repmat(rd.',Nb,size(x,2)).*reshape(x,Nb,n3*size(x,2))),[],size(x,2)); %#ok<UNRCH>
        IT0 = @(x) reshape(Kernels.ITSSD0*reshape(x,Nb,n3*size(x,2)),[],size(x,2));
        A = @(x) real(F.'*(Ck*(S0(-IT0(Lapp(TD,Bf(x))+Lk*Bf(x))+Bf(x)))));
    else 
        A = @(x) real(F.'*(Ck*Lapp(SD,Lslv(TD,-Lapp(TD,Bf(x))+Lk*Bf(x),parslv)+Bf(x))));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build constant vector b: 

% Compute (1/dt)*phi
phib = zeros(numF+numFS,1); 
if numF>0
    % phi(i,j) = |C_i-C_j|-(r_i+r_i)-eps*max(r_i,r_j)
    phib(1:numF) = (1/dt)*(NR-diam(ip+n3*(jp-1))-eps*mxrd(ip+n3*(jp-1))); 
end
if numFS>0
    NC = sqrt(sum(Ct(ipsh,:).*Ct(ipsh,:),2)); 
    % phi(i,shell) = (R-r_i) - sheps*rdsh - |C_i|
    distSh = (Fparams.parsh.rd - rd(ipsh) - sheps*Fparams.parsh.rd) - NC; 
    phib(numF+1:numF+numFS) = (1/dt)*distSh;  
end

%b_k = (1/dt)*phi_k + F.'V_k
bvec = phib + real((F.')*VW(:));
%TODO: add options for restitution / elastic collisions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LCP solve
%LCP params
max_iter=parslv.colmaxit; 
tol_rel=parslv.col_tolrel; %1e-6; 
tol_abs=parslv.col_tolabs; %1e-9; 
profile=1;

x0 = zeros(size(bvec));

opts = struct( ...
    'max_iter', max_iter, ...
    'tol_rel', tol_rel, ...
    'tol_abs', tol_abs ...
);
fg = @(x) objGrad(x, A, bvec);
switch parslv.colsolver  
    case 'BBPGD'
        [lam, info] = projectedGradientDescent(fg, x0, opts);    
    case 'L-BFGS-B'
        [lam, info] = L_BFGS_B(fg, x0, opts); 
    case 'P-L-BFGS'
        [lam, info] = projectedQuasiNewton(fg, x0, opts);    
    case 'proxQuasiNewton'
        [lam, info] = proxQuasiNewton(fg, x0, opts);
    otherwise
        [lam, info] = projectedGradientDescent(fg, x0, opts); 
end

if saveLCPs
    if length(A_list) >= 100 || endFlag
        mfilePath = mfilename('fullpath');
        if contains(mfilePath,'LiveEditorEvaluationHelper')
            mfilePath = matlab.desktop.editor.getActiveFilename;
        end
        [mfilePath,~,~] = fileparts(mfilePath);
        save([mfilePath '/LCPSolvers/test/3x3x3_amphi_data_' num2str(save_iter) '.mat'], ...
            'A_list', 'b_list')
        A_list = {};
        b_list = {};
        save_iter = save_iter + 1;
        if endFlag 
            save_iter = 1;
        end
    end
    
    A_list{end+1} = Amat; 
    b_list{end+1} = bvec;
end
fprintf(['\n minmap ' parslv.colsolver ' LCP solution error = %e, iters = %d \n'], info.kkt, info.iter);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contact forces and modified densities

if norm(lam)>0
    % Contact forces / torques
    F_c = F*lam; 
    % Obtain rho and mu densities
    [mu_c,rho_c] = Lapp_ctmat(TD,SD,Lk,Bk.',[],F_c,parslv);
else
   mu_c=[]; rho_c=[]; F_c=[];  
end

end

function [Mf,Bf] = Lapp_ctmat(TD,SD,Lk,BMR,VNS,f,parslv)

Bf = BMR*f; 
%(3) (-0.5I-K)*rho_c
MNS = -Lapp(TD,Bf)+Lk*Bf;
%(4) 3x3 MNS=VNS*S(mu_c+rho_c)  
MuNS = Lslv(TD,MNS,parslv); 

if ~isempty(VNS)
    % Measure difference between velocities at contact pts, given a force
    % density ( ui-uj = MNS*F_c ) 
    Mf = VNS*Lapp(SD,MuNS+Bf); 
else
    Mf = MuNS; 
end

end

function [sigma,mu,U,VW] = LOCAL_compute_velocities(sigma,VW,Ct,Kernels,Nullsp,Fparams,col,collist,i,dt)
%global rst maxit acc timings rd; 
global timings; 
parslv = Fparams.parslv; 
rd = Fparams.parbd.rd; tau = Fparams.parbd.tau; W = Fparams.parbd.W; 
n3 = Fparams.parbd.n3; 
rdt = repmat((rd.').^(-4),3,1); 
rdw = repmat((rd.').^(-2),3,1);  
shflg = isfield(Fparams,'parsh'); 
vind = reshape(repmat(6*(0:n3-1),3,1),1,[])+repmat((1:3),1,n3);
wind = reshape(repmat(6*(0:n3-1),3,1),1,[])+repmat((4:6),1,n3);

% If inside shell, add traction from boundary correction to rhs
Fparams.Tshell = []; 
if shflg
    if ~isempty(Fparams.parsh.shellden)
        fprintf('\n Evaluation of traction correction field from shell:')
        Fparams.parsh.flag_pot = 'TSL_Stk_3D'; 
        Xtrg = Fparams.parbd.Xp; Nrtrg = Fparams.parbd.Nrp; 
        tic; 
        Fparams.Tshell = real(VSh_MatVec_RB_trg(Fparams.parsh.shellden,Xtrg,Nrtrg,Fparams.parsh));
        fprintf('\n Time for Traction computation from boundary: %e',toc); 
        
        parshS = Fparams.parsh; parshS.flag_pot = 'SL_Stk_3D'; 
        Ush = real(VSh_MatVec_RB_trg(parshS.shellden,Xtrg,Nrtrg,parshS)); 
    end
end

if isempty(VW) || Fparams.comp
% Fluid Solve
% (1) U_inc=S[sigma] (particular solution given forces and torques)
% Sigma is the incoming traction distribution 

% (2) U_sc=S[mu] ("scattered" field with zero forces and torques)
% RHS -(aI+K)*sigma
tic; 
B = Nullsp.L*sigma-Lapp(Kernels.TD,sigma);
timings.velocities.apply(i) = 0.5*toc;

if ~isempty(Fparams.Tshell)
B = B - Fparams.Tshell; 
end

% Solve Fredholm eq TD*mu = B
tic; 
mu = Lslv(Kernels.TD,B,parslv);
fprintf('\n Time for solve: %e',toc); 
timings.velocities.solve(i) = toc;  
% U = U_inc + U_sc
tic; 
U = Lapp(Kernels.SD,(mu+sigma)); 
fprintf('\n Time for apply (of S) to compute U: %e',toc);  
timings.velocities.apply(i) = timings.velocities.apply(i) + 0.5*toc;

if ~isempty(Fparams.Tshell)
   U = U + Ush; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute (V,W) and advance Ct, Xt and Mt
tic; 
CU = Nullsp.C*U; 
IU = CU(vind); 
WxI = CU(wind); 
    
% U(S_k) = V_k + W_k x (X(S_k) - C_k)
VW = zeros(6,n3); 
if ~iscell(W)
    VW(1:3,:) = (1/sum(W))*(rdw.*reshape(IU,3,n3)); 
    VW(4:6,:) = tau\(rdt.*reshape(WxI,3,n3)); 
else
    IUv = reshape(IU,3,n3); WxIv = reshape(WxI,3,n3); 
    for j=1:n3
        VW(1:3,j) = (1/sum(W{j}))*IUv(:,j); 
        VW(4:6,j) = tau{j}\WxIv(:,j);
    end
end
fprintf('\n Time to compute V and W: %e',toc); 
timings.velocities.vw(i) = toc; 

errbs = norm(U-Nullsp.D'*VW(:))/norm(U);
fprintf('\n Rigid body velocity error: %e',errbs)

elseif col || ~isempty(Fparams.Tshell) 
    tic; 
    B = Nullsp.L*sigma-Lapp(Kernels.TD,sigma);
    timings.velocities.apply(i) = toc;
    
    % If inside shell, add traction from boundary correction to rhs
    if ~isempty(Fparams.Tshell)
        B = B - Fparams.Tshell; 
    end
    
    % Solve Fredholm eq TD*mu = B
    tic; 
    mu = Lslv(Kernels.TD,B,parslv); 
    timings.velocities.solve(i) = toc;
    
    if ~col
       U   = Lapp(Kernels.SD,(mu+sigma)) + Ush; 
   
       CU  = Nullsp.C*U; 
       IU  = CU(vind); 
       WxI = CU(wind); 
      
       if ~iscell(W)
            VW(1:3,:) = (1/sum(W))*(rdw.*reshape(IU,3,n3)); 
            VW(4:6,:) = tau\(rdt.*reshape(WxI,3,n3)); 
       else
           IUv = reshape(IU,3,n3); WxIv = reshape(WxI,3,n3); 
           for j=1:n3
               VW(1:3,j) = (1/sum(W{j}))*IUv(:,j); 
               VW(4:6,j) = tau{j}\WxIv(:,j);
           end
       end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collision event
if col
   tic; 
    
   fprintf('\n-------------------------------------------------');
   fprintf('\n Computing contact force at collision sites: \n')
   display(collist(:,1:2)')
   fprintf('-------------------------------------------------\n');
   
   % Compute contact force and force distribution updates
   [F_c,mu_c,rho_c] = LOCAL_Compute_Contact_LCP(collist,Kernels,Nullsp,Fparams,Ct,VW,dt);
   
   F_c = reshape(F_c,6,[]);
   display(F_c(1:3,:));
   
   if ~isempty(mu_c)
       %Update sigma, mu, U and VW
       mu   = mu    + mu_c;
       sigma = sigma + rho_c;
       U   = Lapp(Kernels.SD,(mu+sigma));

       if ~isempty(Fparams.Tshell)
           U = U + Ush;
       end

       CU  = Nullsp.C * U;
       IU  = CU(vind);
       WxI = CU(wind);

       if ~iscell(W)
           VW(1:3,:) = (1/sum(W))*(rdw.*reshape(IU,3,n3));
           VW(4:6,:) = tau\(rdt.*reshape(WxI,3,n3));
       else
           IUv = reshape(IU,3,n3); WxIv = reshape(WxI,3,n3);
           for j=1:n3
               VW(1:3,j) = (1/sum(W{j}))*IUv(:,j);
               VW(4:6,j) = tau{j}\WxIv(:,j);
           end
       end
   end
   
   fprintf('\n Time for contact force correction: %e',toc);  
   timings.velocities.col(i) = toc;
      
   % Check rigid body velocity
   errbs = norm(U-Nullsp.D'*VW(:))/norm(U);
   fprintf('\n Rigid body velocity error after collision correction: %e',errbs)
end

VW = real(VW); 

end

function params = LOCAL_compute_shell_velocity(mu,params)

shf = @(q) VshAna([q(1:3:end);q(2:3:end);q(3:3:end)],'VW'); 

fprintf('\n Evaluation of flow on shell and boundary correction:')
Xtrg = params.parsh.Xp; Nrtrg = params.parsh.Nrp; 
Usqr = real(VSh_MatVec_RB_trg(mu,Xtrg,Nrtrg,params.parbd));
Uh = shf(Usqr); 
sigmah = -params.parsh.eigI.*Uh; 
sigma  = real(VshSyn(sigmah,'VW'));
sigma = reshape(reshape(sigma,[],3).',[],1);

% Set params for shell density apply
params.parsh.U = Usqr; 
params.parsh.Uh = Uh; 
params.parsh.shellden = sigma; 
params.parsh.Vh = sigmah; 

params.parsh.flag_pot = 'SL_Stk_3D'; params.parsh.a = 0; 

Uinf = Usqr + VSh_MatVec_RB2(sigma,[],params.parsh);  
fprintf('\n Check flow at boundary = 0: ||Uinf||_2 = %0.5g , ||Uinf||_inf = %0.5g',norm(Uinf),max(abs(Uinf))) 
if max(abs(Uinf))>1e-3
   plot(log10(abs(shf(Uinf)))); pause(0.0001);   
end

end

function Ctp = LOCAL_advance_center(Ct,dt,VW,Fparams)
%global type; 

if strcmp(Fparams.type,'Purcellxy')
    th = Fparams.tht + dt*VW(6,:);
    th = th - 2*pi*(floor(th./(2*pi))); 
    Fparams.tht = th; 
    Ctp = zeros(3,3); 
    Ctp(:,1:2) = [Rer(Fparams.tht(1));Rer(Fparams.tht(2));Rer(Fparams.tht(3))];  
else
    Ctp = Ct + dt*VW(1:3,:)';
end

end

function M = RotationMat(wh,t)

nwh = norm(wh); 
t = nwh*t; 
wh = wh./nwh; 

M = [1-(wh(2)^2+wh(3)^2)*(1-cos(t)),wh(2)*wh(1)*(1-cos(t))-wh(3)*sin(t),wh(1)*wh(3)*(1-cos(t))+wh(2)*sin(t);...
wh(1)*wh(2)*(1-cos(t))+wh(3)*sin(t),1-(wh(1)^2+wh(3)^2)*(1-cos(t)),wh(2)*wh(3)*(1-cos(t))-wh(1)*sin(t);...
wh(1)*wh(3)*(1-cos(t))-wh(2)*sin(t),wh(2)*wh(3)*(1-cos(t))+wh(1)*sin(t),1-(wh(2)^2+wh(1)^2)*(1-cos(t))];

end

function [Mtp,Xtp,nrmW] = LOCAL_advance_rotation(MRot,VW,Mt,Xt,X0,dt,np,n3)

nrmW = zeros(n3,1); Mtp=Mt; Xtp=Xt;  
for k=1:n3
    xind = (1:np)+np*(k-1); 
    nrmW(k) = norm(VW(4:6,k)); 
    if nrmW(k)>1e-10
        
        Mtp{k} = MRot(VW(4:6,k),dt)*Mt{k};          
        Xtp(xind,:) = X0(xind,:)*Mtp{k}';
    else   
        Mtp{k} = Mt{k}; 
        Xtp(xind,:) = Xt(xind,:);    
    end   
end
end

function [colevent,collist,dt,Ctp] = LOCAL_collision_info(Fparams,X2,Ct,Ctp,VW,MRot,Mt,dt)
%global diam rd mxrd sheps;  
n3 = size(Ct,1); Shape=Fparams.parbd.Shape; 
eps = Fparams.parbd.eps;  

%Check for collision between spheres (or sphere envelopes)
[colevent,collist,mindst,mindstsh] = LOCAL_check_collision_sph(Ctp,Fparams);

% Finer collision detection for non-spheres
if ~strcmp(Shape,'')
    np2 = size(X2,1)/n3;
    [colevent,collist,mindst,Xip,Xjp]=LOCAL_check_collision(collist,eps,Ctp,X2,MRot,Mt,VW,dt,np2); 
end

if colevent
    
    % If mindst<<eps, or <0, we need to adjust timestep
    bis=0; maxbis=3; shell = isfield(Fparams,'parsh');
    
    cond = mindst < 0.1*eps;
    if shell
        sheps = Fparams.parsh.eps;
        cond = cond || mindstsh < 0.1*sheps;
    end
    
    while cond && bis<=maxbis
    
        % Recompute dt and Ct{i+1} to avoid collision
        dt = dt/2;
        bis=bis+1; 
        Ctp = LOCAL_advance_center(Ct,dt,VW,Fparams); 
    
        %Check for collision between spheres (or sphere envelopes)
        [colevent,collist,mindst,mindstsh] = LOCAL_check_collision_sph(Ctp,Fparams);

        % Finer collision detection for non-spheres
        if ~strcmp(Shape,'')
            [colevent,collist,mindst,Xip,Xjp]=LOCAL_check_collision(collist,eps,Ctp,X2,MRot,Mt,VW,dt,np2); 
        end
        
        fprintf('\n bisection = %d: dt = %1.4e, mindst = %1.4e, mindstsh = %1.4e',bis,dt,mindst,mindstsh);
        
        % update condition
        cond = mindst < 0.1*eps;
        if shell
            sheps = Fparams.parsh.eps;
            cond = cond || mindstsh < 0.1*sheps;
        end
    end
    
    fprintf('\n Min pairwise relative distance after collision detection: %2.4f ',mindst);
    if shell
        fprintf('\n Min distance to geometry after collision detection: %2.4f ',mindstsh);
    end 
    
else
    collist=[]; 
    fprintf('\n Min pairwise relative distance: %2.4f',mindst);
    if isfield(Fparams,'parsh')
        fprintf('\n Min relative distance to geometry: %2.4f, absolute distance: %2.4f ',mindstsh,mindstsh*Fparams.parsh.rd);
    end
end

end

function [colevent,collist,mindst,mindstsh] = LOCAL_check_collision_sph(C,Fparams)

n3 = size(C,1); 

rd = Fparams.parbd.rd; 
mxrd = Fparams.parbd.mxrd; 
diam = Fparams.parbd.diam; 
eps = Fparams.parbd.eps;  

if n3>1
    % Compute center distances
    distC = LOCAL_CenterDistance(C);

    % Find pairs for which (C_i-C-j) <= (r_i+r_j)+1.1*eps*max(r_i,r_j)
    [ii,jj]=meshgrid(1:n3); 
    id = distC<=diam+1.1*eps*mxrd & ii<jj; 
    ip = ii(id); 
    jp = jj(id); 

    % Compute minimum relative distance between spheres
    mindst=min(reshape((distC-diam)./mxrd,[],1));
else
    ip=[]; jp=[]; 
    mindst=Inf; 
end

% If there is a spherical shell, compute signed distance to boundary
if isfield(Fparams,'parsh')
   rdsh = Fparams.parsh.rd; 
   sheps = Fparams.parsh.eps;
     
   NC = sqrt(sum(C.*C,2)); 
   distSh = rdsh - rd - NC; %(R-r) - ||C_i||
   mindstsh = min(distSh)/Fparams.parsh.rd; 
   iish = find(distSh <= 1.1*sheps*rdsh); 
   jjsh = (n3+1)*ones(size(iish)); %j=n3+1 -> collision with boundary
   
   colevent = mindst < 1.1*eps || mindstsh < 1.1*sheps;
   collist = [[ip ; iish] [jp ; jjsh]];
else
   colevent = mindst < 1.1*eps; 
   mindstsh = Inf;
   collist = [ip jp]; 
end

end

function [colevent,collist,mindst,Xip,Xjp]=LOCAL_check_collision(collist,eps,C,X,MRot,Mt,VW,dt,np)
% WARNING: This function needs to be updated to work w/ non-sphere objects 
idpw = collist(:,2)<size(C,1)+1;
ip = collist(idpw,1); jp = collist(idpw,2); 
distcol=zeros(size(ip)); 

Xip = zeros(length(ip),3); Xjp=Xip; 
        
for k=1:length(ip)
   xiind = (1:np)+np*(ip(k)-1); 
   Xi = repmat(C(ip(k),:),np,1) + X(xiind,:)*(MRot(VW(4:6,ip(k)),dt)*Mt{ip(k)}); 
   xjind = (1:np)+np*(jp(k)-1); 
   Xj = repmat(C(jp(k),:),np,1) + X(xjind,:)*(MRot(VW(4:6,jp(k)),dt)*Mt{jp(k)}); 
   dij = LOCAL_Distance(Xi,Xj);  
   [distcol(k),iik] = min(dij(:)); 
   jcol = ceil(iik/np); icol=iik-np*(jcol-1); 
   Xip(k,:) = Xi(icol,:); 
   Xjp(k,:) = Xj(jcol,:); 
end
        
id = distcol<=eps; 
colevent = ~isempty(id);

ip = ip(id); 
Xip = Xip(id,:); 
jp = jp(id); 
Xjp = Xjp(id,:); 
fprintf('\n Check collision \n')
display([ip jp distcol(id)])
mindst=min(distcol); 

% Have to add collision with wall to make this proper
collist = [[ip jp];collist(~idpw,:)]; 

end
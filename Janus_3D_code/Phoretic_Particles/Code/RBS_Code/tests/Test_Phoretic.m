function [Err,Utrg,Xtrg,U,Us,nrmX]=Test_Phoretic(fname,exp,p,C,M0,Nt,dt,mdist,tdisc,beta,flplot,A,M_theta,shell_data) 
%{
This function can run one of several tests, depending on the string 'exp': 
(0) 'shell' - Test of interior Dirichlet problem for Stokes for the outer
shell (currently not functional) 
(1) 'onesquirmer' - Test for one squirmer body in free flow. Velocity
should be (3/2)*B1 in the z direction (or whichever direction prescribed by
rotation matrix) (Not Implemented for Phoretic particles)
***************************************************************************
***************************************************************************
(2) 'free_phoretic' - Test for one or multiple phoretic particles bodies in free Stokes
flow
(3) 'fullsystem' - Test for one or multiple squirmer bodies in confined
flow inside a shell of radius shrd. The shell center is assumed to be
centered at the origin. For this case, parameters psh, shrd and epsh are
currently hard-coded. 

Inputs: 
fname    (string)        file header 
exp      (string)        experiment name ('shell','onesquirmer','squirmers','fullsystem')
p        (int)           spherical harmonics order for particles
C        (double ns x 3) particle center array
M0       (cell ns x 1)   initial rotation matrix cell. If empty, assumed to be M0{i}=eye(3) for i=1...ns
Nt       (int)           number of timesteps 
dt       (double)        initial timestep size
mdist    (double)        if centers are < this distance, particles are neighbors
tdisc    (string)        time discretization (euler,rk2,rk4)
beta     (double)        squirmer parameter. beta>0, beta=0, beta<0 are pullers, treadmill and pusher swimmers. 
flplot   (bool)          plot spheres every time-step
A        @(theta)        function handle with argument theta
M        @(theta)        function handle for computing Us
 
Outputs: 
Depending on the case, outputs are all empty or they return arrays needed
for debugging. 

Relevant outputs and timings profile are saved automatically using the
header provided in fname (appending a suffix with simulation-specific
parameters). A video or snapshots of the simulation can then be produced
using plotexample.m 
%}

addpath ../; 
addpath ../LCPsolvers/Num4LCP_MatLab/; 
addpath ../../../../FMMLIB/fmmlib3d-1.2/matlab;
addpath ../../../../FMMLIB/stfmmlib3d-1.2/matlab/;

Err = []; Utrg = []; Xtrg = []; nrmX = []; U = []; Us = []; 

switch exp
    case 'shell'
    params = LOCAL_setparams(p,[0 0 0],1,'SL_Stk_3D','',[],1,0,mdist,0);  

    %Test #1: Outer spherical shell
    [Err,Utrg,Xtrg,nrmX]=LOCAL_test_outershell(p,params,'stokeslets',flplot);  

    case 'onesquirmer'
    %Test #2.0: One squirmers at the origin
    Csq = [0 0 0]; %squirmer center at origin
    rd = 1; %squirmer radii
    params.parsq = LOCAL_setparams(p,Csq,rd,'DL_Stk_3D','',M0,1,1,mdist,1);  
    params.parsq.a=0.5; 
    [Err,Utrg,U,Us,Xtrg] = LOCAL_test_onesquirmer(p,params,beta,flplot);    
    
    case 'free_phoretic'
    %Test #2: Squirmers in free space
    Csq = C; %squirmer centers
    rd = 1*ones(size(C,1),1); %squirmer radii
    eps=0.3; %collision epsilon
    params.parsq = LOCAL_setparams(p,Csq,rd,'DL_Stk_3D','',M0,1,0,mdist,1);
    params.parsq.a=0.5; 
    params.A_label = A; params.M_label = M_theta;
    params.parLap = LOCAL_setparams(p,Csq,rd,'dSL_L_3D','',M0,1,1,mdist,1,1);
    params.parLap.a = -0.5; params.parsq.Mt = M0;
    
    LOCAL_test_free_phoretic(fname,p,params,M0,Nt,dt,eps,tdisc,beta,flplot);

    case 'fullsystem'
    %Test #3: Squirmers inside spherical shell of radius shrd
    shrd = shell_data.r;
    psh = shell_data.p;
    Ldense = shell_data.Ldense;
    Sdense = shell_data.Sdense;
    
    %generate shell spectra for S':
    Spectra = zeros((psh+1)^2,1);
    index = 1;
    for n=0:8
        for m = 1:2*n+1
            Spectra(index) = n/(2*n+1);
            index = index + 1;
        end
    end
             
    Csq = C; %squirmer centers
    rd = 1*ones(size(C,1),1); %squirmer radii
    eps=0.3; epsh=0.03; 
    params.parsq = LOCAL_setparams(p,Csq,rd,'DL_Stk_3D','',M0,Sdense,1,mdist,1);
    params.parsh = LOCAL_setparams(psh,[0 0 0],shrd,'SDL_Stk_3D','',[],1,1,0.75,0);
    params.parsh.a=-0.5; 
    params.parLap = LOCAL_setparams(p,Csq,rd,'dSL_L_3D','',M0,Ldense,1,mdist,1,1);
    params.parLap.a = -0.5; params.parsq.Mt = M0;
    params.parLapsh = LOCAL_setparams(psh,[0 0 0],shrd,'dSL_L_3D','',[],0,1,0.75,0,1);
    params.parLapsh.a= 0.5; params.parLapsh.Spectra = Spectra;
  %change this to 0 if you want to use fmm
    params.A_label = A; params.M_label = M_theta;
    
   
    ns = size(params.parsq.C,1);
    np = params.parsq.np; %size(params.parsq.Xrp,1)/ns;
    
    
    %integrates A(theta) over the particles to compute flux from particles
    init_dir = repmat([0 0 1], ns,1);
    surfacelabel = zeros(size(params.parsq.Xrp,1),1); 
    for i = 1:ns
       Xsphere=params.parsq.Xrp(np*(i-1)+1:np*i,:);
       Xbackrot=Xsphere*params.parsq.Mt{i}; 
       surfacelabel(np*(i-1)+1:np*i,:)=params.A_label(acos(Xbackrot*init_dir(i,:)'));
    end
    particle_flux = -params.parLap.W2*surfacelabel;
    params.parLapsh.shell_flux = particle_flux/(4*pi*shrd^2);
    
    LOCAL_test_fullsystem(fname,p,psh,shrd,params,M0,Nt,dt,eps,epsh,tdisc,beta,flplot);

end
end

function [ES,Utrg,Xtrg,nrmX] = LOCAL_test_outershell(p,params,flow,flplot)

% Outer shell setup 
shf = @(q) VshAna([q(1:3:end);q(2:3:end);q(3:3:end)],'VW'); 
sp=(p+1)^2; 
ii = (1:sp)'; 
nn = floor(sqrt(ii-1));

% Generate input flow u_inf Vsph density 
[Uh,params] = LOCAL_Uinf_Vsph(p,flow,params); 
%figure; plot(log10(abs(Uh)),'ob')

%Invert Stokes (first kind) 
out = 0; %interior flow
[SV,SW,SX] = LOCAL_get_eigen('SMat',out); 
eigS = [SV(nn);SW(nn);SX(nn)]; 
eigIS = 1./(eigS+(eigS==0)) - (eigS==0); 

sigmah = eigIS.*Uh; 
sigma  = VshSyn(sigmah,'VW');
sigma = reshape(reshape(sigma,[],3).',[],1);

if flplot
figure; plot(log10(abs(sigmah)),'or');
end

if isfield(params,'Q') && flplot
figure; plot(real(params.Vh),'ob'); hold on; plot(real(sigmah),'or'); hold on; plot(real(Uh),'ok');     
figure; plot(real(params.Q)); hold on; plot(real(sigma),'r'); 
end

%figure; plot(params.U,'k');

%Invert -1/2 I + S + D (second kind) 
[SDV,SDW,SDX] = LOCAL_get_eigen('SDMat',out);
eigSD = [SDV(nn);SDW(nn);SDX(nn)]; 
eigISD = 1./(eigSD+(eigSD==0)) - (eigSD==0); %1./eigSD; 

muh = eigISD.*Uh; 
mu  = VshSyn(muh,'VW'); 
mu = reshape(reshape(mu,[],3).',[],1);

%Evaluate at target points inside the sphere and compare
ptrg=4; nstrg = 9; 
np8 = 2*ptrg*(ptrg+1);
Sc = SurfaceSph(shape_gallery(ptrg,'')); 
% Model surface pts
Xp = reshape(Sc.cart.to_array,[],3);
Xtrg = repmat(Xp,nstrg,1); 
%Xtrg = [0.99999*Xp; 0.9999*Xp ; 0.999*Xp; 0.99*Xp ; 0.9*Xp; 0.75*Xp; 0.5*Xp; 0.25*Xp ; 0.1*Xp];
%Xtrg = [0.9999999*Xp ; 0.875*Xp; 0.75*Xp; 0.625*Xp; 0.5*Xp; 0.25*Xp ; 0.1*Xp];
Xtrg = (1-rand(size(Xtrg)).^6).*Xtrg;
nrmX = sqrt(sum(Xtrg'.*Xtrg')'); 
display(min(nrmX))
display(max(nrmX))
Nr = reshape(Sc.geoProp.nor.to_array,[],3); Nrtrg = repmat(Nr,nstrg,1); 

parS = params; parS.flag_pot = 'SL_Stk_3D'; parS.Vh = sigmah; 
[US,~,~]  = VSh_MatVec_RB_trg(sigma,Xtrg,Nrtrg,parS);

parSD = params; parSD.flag_pot = 'SDL_Stk_3D'; parSD.Vh = muh; %parS.Vh = muh; 
[USD,~,~] = VSh_MatVec_RB_trg(mu,Xtrg,Nrtrg,parSD);

% True velocity from example
Utrg = LOCAL_Uinf_trg(Xtrg,Nrtrg,flow,params); 

%Error
ES = abs(US-Utrg);
display(norm(ES)/norm(Utrg)); 
if flplot
figure; plot(log10(ES));
end
%figure; plot(abs(shf(Utrg(end-3*np8+1:end,:))),'ok');
%figure; plot((abs(shf(US(end-3*np8+1:end,:)))),'or');

if flplot
colrg = [min(abs(Utrg(:))) max(abs(Utrg(:)))]; 
for k=nstrg:-1:1
ind = np8*(k-1)+1:np8*k;
indV = 3*np8*(k-1)+1:3*np8*k;
LOCAL_plotsph(Xtrg(ind,:),abs(Utrg(indV)),colrg); 
hold on; 
pause; 
end
end

%{
%params.flag_pot = 'SL_Stk_3D'; 
%Xtrgv = reshape(repmat(0.1*Xp,1,3)',3,[])';  
%Xv = reshape(repmat(Xp,1,3)',3,[])';
%Uctr = Kernel_Eval(Xtrgv,Xv,params)*sigma; 
%norm(Uctr - Utrg(end-3*np8+1:end,:))
%figure; plot((abs(shf(Uctr(end-3*np8+1:end,:)))),'oc');
%}
 
ESD = abs(USD-Utrg); 
display(norm(ESD)/norm(Utrg));
if flplot
figure; plot(log10(ESD));  
end
%figure; plot(abs(shf(Utrg(end-3*np8+1:end,:))),'ok'); 
%figure; plot((abs(shf(USD(end-3*np8+1:end,:)))),'og'); 

%{
%params.flag_pot = 'SDL_Stk_3D'; 
%Uctr = Kernel_Eval(Xtrgv,Xv,params)*mu; 
%norm(Uctr - Utrg(end-3*np8+1:end,:))
%figure; plot((abs(shf(Uctr(end-3*np8+1:end,:)))),'oc');
%}

end

function [ES,Utrg,U,Us,Xtrg] = LOCAL_test_onesquirmer(p,params,beta,flplot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial params / setup
shf = @(q) VshAna([q(1:3:end);q(2:3:end);q(3:3:end)],'VW');

% Slip velocity parameters
B1 = 3/2; B2 = beta*B1; 
params.parsq.B1 = B1; params.parsq.B2 = B2; params.parsq.beta = beta; 

np = 2*p*(p+1); 
ns = 1; 
Xrp  = params.parsq.Xrp;
Xp   = params.parsq.Xp; 
Wg   = params.parsq.Wg; 
X0 = Xrp;   
params.parsq.a = 0.5; 

%create auxiliary matrices
[Cm,~,Dm] = LOCAL_BuildC(Wg,Xrp,[],np,ns); 
Kernels.C = full(Cm); Kernels.D = full(Dm); 

%Assume for now params.dense=1
A = real(LOCAL_MatVec([],params.parsq)); 
% Build completion flow matrix (rank 6*ns matrix)
parStk = params.parsq; parStk.flag_pot='SL_Stk_3D'; parStk.a=0; parStk.cj=(1:3)';
parStk.W2 = ones(1,3*ns);
C0 = repmat(params.parsq.C,3,1); 
Gij = Kernel_Eval(parStk.X,C0,parStk); %Stokeslet
parRot = parStk; parRot.flag_pot='Rotlet_3D'; 
Rij = Kernel_Eval(parRot.X,C0,parRot); %Rotlet

Km = zeros(size(Gij,1),2*size(Gij,2)); 
for k=1:ns
    Km(:,(1:6)+6*(k-1)) = [Gij(:,(1:3)+3*(k-1)) Rij(:,(1:3)+3*(k-1))]; 
end
Kernels.Km=Km;
Kernels.SDL = A;

CF = Km*Cm; %rank 6 PM completion flow
ss = svd(A+CF); 
figure; plot(log10(ss),'-o'); figure; 
display(sum(ss>1e-10))

M = [A+CF -Dm.' ; Cm zeros(6*ns,6*ns)]; 
Kernels.M = M; 

Schr = Cm*((A+CF)\Dm.'); 
%Prc = zeros(size(M)); 
%Prc(1:3*np,1:3*np)=inv(A+CF); Prc(3*np+1:end,3*np+1:end)=inv(Schr); 

fprintf('\n Condition number of 0.5*I+D+N = %1.2e',cond(A+CF));
fprintf('\n Condition number of full system = %1.2e',cond(M)); 
fprintf('\n Condition number Schur Complement = %1.2e',cond(Schr)); 
%fprintf('\n Cond number prec M = %1.2e',cond(Prc)); 

% Distance between centers
distC = LOCAL_CenterDistance(params.parsq.C);
rsum = repmat(params.parsq.rd,1,ns); rsum = rsum+rsum.'; 
distC = distC - rsum; 
mindst = min(distC(:)); 
fprintf('\n Minimum distance (spheres): %2.4f ',mindst);

colevent=false; 
colinfo = struct('col',colevent,'ip',[],'jp',[],'mindst',mindst,'eps',eps,'type','LCP');

Nt=1; 
Xt = cell(Nt+1,1); Ct=Xt;  

% Evolution  
Xt{1}=Xrp; Ct{1}=[0 0 0];
Mt = cell(Nt+1,ns); 
for k=1:ns
    Mt{1,k}=eye(3);   
end

%[Xtp,Mtp,Ctp,U,Us,mu,VW,Kernels,params,colinfo,dt]
fprintf('\n Explicit euler step \n')
[~,~,~,U,Us,mu,~,~,~,~,~] = ...
                LOCAL_euler_step(Xt{1},X0,Mt(1,:),Ct{1},Kernels,params,colinfo,0.1,1,flplot);

Xtrg = [Xp ; 2*Xp ; 1.5*Xp; 1.25*Xp; 1.1*Xp; 1.01*Xp; 1.0001*Xp];
Nrtrg = repmat(params.parsq.Nrp,7,1);             
            
% Target evaluation 
parSD = params.parsq; parSD.Vh = shf(mu); 
[Utrg,~,~] = VSh_MatVec_RB_trg(mu,Xtrg,Nrtrg,parSD);
Utrg = real(Utrg); 

Xtrgv = reshape(repmat(Xtrg,1,3)',3,[])';parStk.ci=repmat((1:3)',size(Xtrg,1),1);
parStk.W2 = ones(1,3*ns); 
Gij = Kernel_Eval(Xtrgv,C0,parStk); %Stokeslet
parRot = parStk; parRot.flag_pot='Rotlet_3D'; 
Rij = Kernel_Eval(Xtrgv,C0,parRot); %Rotlet

Kmtrg = zeros(size(Gij,1),2*size(Gij,2)); 
for k=1:ns
    Kmtrg(:,(1:6)+6*(k-1)) = [Gij(:,(1:3)+3*(k-1)) Rij(:,(1:3)+3*(k-1))]; 
end
Utrg = Utrg + Kmtrg*(Kernels.C*mu); 

fprintf('\n Self eval \n')
display(norm(U - (VSh_MatVec_RB2(mu,[],parSD) + Kernels.Km*(Kernels.C*mu)))); 
display(norm(U - (Utrg(1:3*np))));

% True velocity for squirmer
Utrue = LOCAL_Usq1_trg(Xtrg,params.parsq); 
VW = [0 0 2*B1/3 0 0 0].'; 

% Error 
Utrue = reshape(Utrue.',[],1);
Utrue = Utrue + repmat(Dm.'*VW,7,1);
ES = Utrue-Utrg;
display(norm(ES)/norm(Utrue));
if flplot
figure; plot(log10(abs(ES)));  
end 

end

function Utrg = LOCAL_Usq1_trg(Xtrg,params)
B1 = params.B1; B2 = params.B2; 

% spherical coordinates
[th,phi,r] = cart2sph(Xtrg(:,1),Xtrg(:,2),Xtrg(:,3));  
th(th<0)=th(th<0)+2*pi; 
phi=pi/2-phi;

x = 1./r; 
Ur = (2*B1/3)*(x.^3-1).*cos(phi) + (B2/2)*(x.^4-x.^2).*(3*cos(phi).^2-1); 
Uph = (2*B1/3)*(0.5*x.^3+1).*sin(phi) + (B2/2)*(x.^4).*sin(2*phi);

%Unit vectors
v = th; u = phi; 
er = [sin(u).*cos(v) sin(u).*sin(v) cos(u)];
eu = [cos(u).*cos(v) cos(u).*sin(v) -sin(u)];

%Add radial and azimuthal velocities
Utrg = repmat(Ur,1,3).*er + repmat(Uph,1,3).*eu; 

end

function LOCAL_test_free_phoretic(name,p,params,M0,Nt,dt,eps,tdisc,beta,flplot)
%rng(41);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial params / setup
np = 2*p*(p+1); 
ns = size(params.parsq.C,1); 
Nmu = 3*np*ns; 
Xrp  = params.parsq.Xrp; 
Wg   = params.parsq.Wg; 
X0 = Xrp; 
params.parsq.a = 0.5;
denseS = params.parsq.dense; 
denseL = params.parLap.dense; 

% Slip velocity parameters
B1 = 3/4; B2 = beta*B1; 
params.parsq.B1 = B1; params.parsq.B2 = B2; params.parsq.beta = beta;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create auxiliary matrices
[Cm,Bm,Dm] = LOCAL_Build_AuxMats(Wg,Xrp,[],np,ns);
if denseS
    Cm = full(Cm); Dm=full(Dm); Bm = full(Bm);  
end
Kernels.C = Cm; Kernels.D = Dm; Kernels.B=Bm; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Laplace matrix build (diffusion) 
Kernels.dSL = LOCAL_compute_BIE_matrix(params,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stokes matrix build 

% Power-Miranda completion flow 
parStk = params.parsq; parStk.flag_pot='SL_Stk_3D'; parStk.a=0; 
parStk.cj=repmat((1:3)',ns,1); parStk.W2 = ones(1,3*ns); 

Cv = reshape(repmat(params.parsq.C,1,3)',3,[])';  
Gij = Kernel_Eval(parStk.X,Cv,parStk); %Stokeslet
parRot = parStk; parRot.flag_pot='Rotlet_3D'; 
Rij = Kernel_Eval(parRot.X,Cv,parRot); %Rotlet

Km = zeros(size(Gij,1),2*size(Gij,2)); 
for k=1:ns
    Km(:,(1:6)+6*(k-1)) = [Gij(:,(1:3)+3*(k-1)) Rij(:,(1:3)+3*(k-1))]; 
end

%rank 6 PM completion flow
Kernels.Km = Km;  

% Stokes A = a*I+S+D matrix and augmented system for PM formulation
if denseS
   A = real(LOCAL_MatVec([],params.parsq)); 
   ACF = A + Km*Cm;
   M = [ACF -Dm.' ; Cm zeros(6*ns,6*ns)]; 
else
   A = @(V) real(LOCAL_MatVec(V,params.parsq)); 
   ACF = @(V) A(V) + Km*(Cm*V);
   M = @(V) [ACF(V(1:Nmu,:))-Dm.'*V(Nmu+1:end,:);Cm*V(1:Nmu,:)];
end

Kernels.SDL = A;
Kernels.ACF = ACF; 
Kernels.M = M; 

if ns==1 && denseS
ss = svd(M); 
figure; plot(log10(ss),'-o'); 
fprintf('\n Condition number of 0.5*I+D+N = %1.2e',cond(ACF));
fprintf('\n Condition number of full system = %1.2e',cond(M));
pause; 
end

% Distance between centers
distC = LOCAL_CenterDistance(params.parsq.C);
rsum = repmat(params.parsq.rd,1,ns); rsum = rsum+rsum.'; 
distC = distC - rsum; 
mindst = min(distC(:)); 
fprintf('\n Minimum distance (spheres): %2.4f ',mindst);

if mindst<=1.1*eps
    fprintf('\n Start at collision\n')
    colevent=true; 
    [ii,jj]=meshgrid(1:ns); 
    id = distC <= 1.1*eps & ii<jj; 
    ip = ii(id); 
    jp = jj(id); 
    [~,midx] = min(distC(id));  
    
    colinfo = struct('col',colevent,'ip',ip,'jp',jp,'mindst',mindst,'eps',eps,'midx',midx,'type','LCP'); 
else
    colevent=false; 
    colinfo = struct('col',colevent,'ip',[],'jp',[],'mindst',mindst,'eps',eps,'type','LCP'); 
end

% Rotation matrix (Rodrigues rotation formula)
MRot = @(wh,t) RotationMat(wh,t); 

t = 0; tt = zeros(Nt+1,1); 
mu = cell(Nt+1,1); U=mu; VW=mu; Xt=mu; Ct=mu; Us = U; 

% Can change to include initial orientations
Mt = cell(Nt+1,ns);
if isempty(M0)
    for k=1:ns
        Mt{1,k}=eye(3);   
    end
else
   Mt(1,:) = M0;  
end

for k=1:ns
   ind = (np*(k-1)+1):np*k;
   Xrp(ind,:) = Xrp(ind,:)*Mt{1,k}; 
end

% Evolution
tt(1)=0; %dt0=dt;  
Xt{1}=Xrp; 
Ct{1}=params.parsq.C; 
dt0 = dt; 

if flplot
    figure; 
end

for i=1:Nt
    dt = dt0; 
    
    switch tdisc
        case 'euler'
            fprintf('\n ---------------------------------------------------------- \n')
            fprintf('\n Explicit euler step \n')
            [Xt{i+1},Mt(i+1,:),Ct{i+1},U{i},Us{i},mu{i},VW{i},Kernels,params,colinfo,dt,Conc{i}] = ...
                LOCAL_euler_step(Xt{i},X0,Mt(i,:),Ct{i},Kernels,params,colinfo,dt,i,flplot); 
            
        case 'rk2'
            fprintf('\n ---------------------------------------------------------- \n')
            fprintf('\n (1) Trapezoidal, predictor step \n')
            [Xt1,Mt1,Ct1,U1,Us1,mu1,VW1,Kernels1,params1,colinfo1,dt] = ...
                LOCAL_euler_step(Xt{i},X0,Mt(i,:),Ct{i},Kernels,params,colinfo,dt,i,0); 
            
            fprintf('\n ---------------------------------------------------------- \n')
            fprintf('\n (2) Trapezoidal, corrector step \n')
            %[U,Us,mu,VW,params] = LOCAL_compute_velocities(Kernels,Mt,Xt,Ct,params,colinfo,dt,it,flplot);
            [U2,Us2,mu2,VW2,params] = LOCAL_compute_velocities(Kernels1,Mt1,Xt1,Ct1,params1,colinfo1,dt,i+1,0);  
            
            % Correct VW{i} as average of VW1 and VW2 (and associated quantities)
            VW{i}    = 0.5*(VW1+VW2); 
            mu{i}    = 0.5*(mu1+mu2);
            U{i}     = 0.5*(U1+U2);
            Us{i}     = 0.5*(Us1+Us2);
            
            [Xt{i+1},Mt(i+1,:),Ct{i+1},Kernels,params,colinfo,~] = LOCAL_advance_step(VW{i},[],Xt{i},X0,Mt(i,:),Ct{i},Kernels,params,colinfo,dt); 
        case 'rk4'
            fprintf('\n ---------------------------------------------------------- \n')
            fprintf('\n (1) Runge-Kutta 4th order, t=t0, x=x(t0) \n')
            dt41 = 0.5*dt; 
            [Xt1,Mt1,Ct1,U1,Us1,mu1,VW1,Kernels1,params1,colinfo1,dt4n] = ...
                LOCAL_euler_step(Xt{i},X0,Mt(i,:),Ct{i},Kernels,params,colinfo,dt41,i,0); 
            
            dt = 2*dt4n; 
            
            fprintf('\n ---------------------------------------------------------- \n')
            fprintf('\n (2) Runge-Kutta 4th order, t_{1/2}=t0+0.5*dt, x=x(t_{1/2}) \n')
            dt42 = 0.5*dt; 
            [U2,Us2,mu2,VW2,params] = LOCAL_compute_velocities(Kernels1,Mt1,Xt1,Ct1,params1,colinfo1,dt42,i+1,0);
            
            [Xt2,Mt2,Ct2,Kernels2,params2,colinfo2,~] = LOCAL_advance_step(VW2,[],Xt{i},X0,Mt(i,:),Ct{i},Kernels,params,colinfo,dt42); 
            
            fprintf('\n ---------------------------------------------------------- \n')
            fprintf('\n (3) Runge-Kutta 4th order, t_{1/2}=t0+0.5*dt, x=x(t_{1/2}) \n') 
            dt43 = dt; 
            
            [U3,Us3,mu3,VW3,params] = LOCAL_compute_velocities(Kernels2,Mt2,Xt2,Ct2,params2,colinfo2,dt43,i+1,0);
            
            [Xt3,Mt3,Ct3,Kernels3,params3,colinfo3,~] = LOCAL_advance_step(VW3,[],Xt{i},X0,Mt(i,:),Ct{i},Kernels,params,colinfo,dt43); 
            
            fprintf('\n ---------------------------------------------------------- \n')
            fprintf('\n (4) Runge-Kutta 4th order, t_{1}=t0+dt, x=x(t_{1}) \n')  
            
            [U4,Us4,mu4,VW4,params] = LOCAL_compute_velocities(Kernels3,Mt3,Xt3,Ct3,params3,colinfo3,dt,i+1,0);
            
            % Average velocities and related quantities using RK4 weights
            VW{i}    = (1/6)*(VW1+2*VW2+2*VW3+VW4); 
            mu{i}    = (1/6)*(mu1+2*mu2+2*mu3+mu4); 
            U{i}     = (1/6)*(U1+2*U2+2*U3+U4); 
            Us{i}     = (1/6)*(Us1+2*Us2+2*Us3+Us4); 
            
            [Xt{i+1},Mt(i+1,:),Ct{i+1},Kernels,params,colinfo,~] = LOCAL_advance_step(VW{i},[],Xt{i},X0,Mt(i,:),Ct{i},Kernels,params,colinfo,dt); 
            
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = t+dt;
tt(i+1)=t;
fprintf('\n Time: %2.2f ',t)

if params.parsq.beta > 0 
    fname = [name '_squirmers_' num2str(ns) '_puller' num2str(params.parsq.beta) ...
        '_p' num2str(p) '_' tdisc '_dt' num2str(dt0) '_mdist' num2str(params.parsq.mdist) '.mat'];
elseif params.parsq.beta < 0 
    fname = [name '_squirmers_' num2str(ns) '_pusher' num2str(abs(params.parsq.beta)) ...
        '_p' num2str(p) '_' tdisc '_dt' num2str(dt0) '_mdist' num2str(params.parsq.mdist) '.mat'];
else
    fname = [name '_squirmers_' num2str(ns) '_treadmill' '_p' num2str(p) ...
        tdisc '_dt' num2str(dt0) '_mdist' num2str(params.parsq.mdist) '.mat'];
end 
                                                                                                                                             
save(fname,'-v7.3','tt','Xt','Mt','Ct','mu','U','Us','VW','Conc');                                                                                         
  
end

end

function LOCAL_test_fullsystem(name,p,psh,shrd,params,M0,Nt,dt,eps,epsh,tdisc,beta,flplot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial params / setup
np = 2*p*(p+1); 
ns = size(params.parsq.C,1); 
Nmu = 3*np*ns; 
Xrp  = params.parsq.Xrp; 
Wg   = params.parsq.Wg; 
params.parsq.a=0.5;
X0 = Xrp; 
denseS = params.parsq.dense; 
denseL = params.parLap.dense; 

%Eigenvalues for Inverse Operator (at outer shell) 
out = 0; %interior flow
potsh = params.parsh.flag_pot; 
if strcmp(potsh(2),'L')
    Mat = [potsh(1) 'Mat']; 
    if strcmp(potsh(1),'S')
        rexp=1;
    elseif strcmp(potsh(1),'D')
       rexp=0;  
    end
elseif strcmp(potsh(1:2),'SD')
    Mat = 'SDMat'; rexp=0; 
end
[SV,SW,SX] = LOCAL_get_eigen(Mat,out); %Layer pot eigenvalues
spsh=(psh+1)^2; 
ii = (1:spsh)'; 
nn = floor(sqrt(ii-1));
eigS = [SV(nn);SW(nn);SX(nn)];
isre = (1/shrd).^rexp; 
eigIS = (isre)./(eigS+(eigS==0)) - (isre).*(eigS==0); 

%Initial outer flow is zero (pass empty density, potential and eigenvalues)
params.parsh.shellden = []; 
params.parsh.Vh = []; 
params.parsh.shellpot = potsh; 
params.parsh.eigI = eigIS; 

% Slip velocity parameters
B1 = 3/4; B2 = beta*B1; 
params.parsq.B1 = B1; params.parsq.B2 = B2; params.parsq.beta = beta;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create auxiliary matrices
[Cm,Bm,Dm] = LOCAL_Build_AuxMats(Wg,Xrp,[],np,ns);
if denseS
    Cm = full(Cm); Dm=full(Dm); Bm = full(Bm); 
end
Kernels.C = Cm; Kernels.D = Dm; Kernels.B = Bm; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Laplace matrix build (diffusion) 
Kernels.dSL = LOCAL_compute_BIE_matrix(params,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stokes matrix build 

% Power and Miranda Completion flow 
parStk = params.parsq; parStk.flag_pot='SL_Stk_3D'; parStk.a=0; 
parStk.cj=repmat((1:3)',ns,1); parStk.W2 = ones(1,3*ns); 

Cv = reshape(repmat(params.parsq.C,1,3)',3,[])';  
Gij = Kernel_Eval(parStk.X,Cv,parStk); %Stokeslet
parRot = parStk; parRot.flag_pot='Rotlet_3D'; 
Rij = Kernel_Eval(parRot.X,Cv,parRot); %Rotlet

Km = zeros(size(Gij,1),2*size(Gij,2)); 
for k=1:ns
    Km(:,(1:6)+6*(k-1)) = [Gij(:,(1:3)+3*(k-1)) Rij(:,(1:3)+3*(k-1))]; 
end

%rank 6 PM completion flow
Kernels.Km = Km;  

% Stokes A = a*I+S+D matrix and augmented system for PM formulation
if denseS
   A = real(LOCAL_MatVec([],params.parsq)); 
   ACF = A + Km*Cm;
   M = [ACF -Dm.' ; Cm zeros(6*ns,6*ns)]; 
else
   A = @(V) real(LOCAL_MatVec(V,params.parsq)); 
   ACF = @(V) A(V)+Km*(Cm*V);
   M = @(V) [ACF(V(1:Nmu,:))-Dm.'*V(Nmu+1:end,:);Cm*V(1:Nmu,:)];
end

Kernels.SDL = A;
Kernels.ACF=ACF; 
Kernels.M = M;

% Distance between centers
Ctr = params.parsq.C; rd = params.parsq.rd; 
distC = LOCAL_CenterDistance(Ctr);
rsum = repmat(rd,1,ns); rsum = rsum+rsum.'; 
distC = distC - rsum; 
mindstrel = min(reshape(distC./repmat(rd(:),1,ns),[],1)); 
mindst = min(distC(:)); 
fprintf('\n (Particles) Min relative distance: %2.4f, Min absolute distance: %2.4f ',mindstrel,mindst);

nrmC = sqrt(Ctr(:,1).^2+Ctr(:,2).^2+Ctr(:,3).^2); 
distCS = shrd-(rd(:)+nrmC); 
mindstsh = min(distCS); 
fprintf('\n (Shell) Min relative distance: %2.4f, absolute distance %2.4f ',mindstsh/shrd,mindstsh);

ctb = 4/3; 

if mindst<=ctb*eps || mindstsh <= ctb*epsh*shrd
    fprintf('\n Start at collision\n')
    colevent=true; 
    [ii,jj]=meshgrid(1:ns); 
    id = distC <= ctb*eps & ii<jj; 
    ip = ii(id); 
    jp = jj(id); 
    [~,midx] = min(distC(id));  
    
    ids = distCS <= ctb*epsh*shrd; 
    ips = ii(ids); jps=(ns+1)*ones(size(ips)); 
    ip = [ip;ips]; jp = [jp;jps]; 
    
    colinfo = struct('col',colevent,'ip',ip,'jp',jp,'mindst',mindst,'mindstsh',mindstsh,'eps',eps,'epsh',epsh,'midx',midx,'type','LCP'); 
else
    colevent=false; 
    colinfo = struct('col',colevent,'ip',[],'jp',[],'mindst',mindst,'eps',eps,'epsh',epsh,'type','LCP'); 
end

% Rotation matrix (Rodrigues rotation formula)
MRot = @(wh,t) RotationMat(wh,t); 

t = 0; tt = zeros(Nt+1,1); 
mu = cell(Nt+1,1); U=mu; VW=mu; Xt=mu; Ct=mu; Us = U; 
Usqr = U; Usqrh = U; gamma=mu; gammah = mu; 

% Evolution
tt(1)=0; %dt0=dt;  
Ct{1}=params.parsq.C;
Mt = cell(Nt+1,ns); 
dt0 = dt; 
 
if isempty(M0)
for k=1:ns
    Mt{1,k}=eye(3);   
end
else
   Mt(1,:) = M0;  
end

for k=1:ns
   ind = (np*(k-1)+1):np*k;
   Xrp(ind,:) = Xrp(ind,:)*Mt{1,k}; 
end

Xt{1}=Xrp; 

if flplot
    figure; 
end

for i=1:Nt
    dt = dt0; 
    
    switch tdisc
        case 'euler'
            fprintf('\n ---------------------------------------------------------- \n')
            fprintf('\n Explicit euler step \n')
            [Xt{i+1},Mt(i+1,:),Ct{i+1},U{i},Us{i},mu{i},VW{i},Kernels,params,colinfo,dt,Conc{i}] = ...
                LOCAL_euler_step(Xt{i},X0,Mt(i,:),Ct{i},Kernels,params,colinfo,dt,i,flplot); 
            
            Usqr{i} = params.parsh.U; Usqrh{i} = params.parsh.Uh; 
            gamma{i} = params.parsh.shellden; gammah{i} = params.parsh.Vh; 
            
        case 'rk2'
            fprintf('\n ---------------------------------------------------------- \n')
            fprintf('\n (1) Trapezoidal, predictor step \n')
            [Xt1,Mt1,Ct1,U1,Us1,mu1,VW1,Kernels1,params1,colinfo1,dt] = ...
                LOCAL_euler_step(Xt{i},X0,Mt(i,:),Ct{i},Kernels,params,colinfo,dt,i,0); 
            
            fprintf('\n ---------------------------------------------------------- \n')
            fprintf('\n (2) Trapezoidal, corrector step \n')
            [U2,Us2,mu2,VW2,params] = LOCAL_compute_velocities(Kernels1,Mt1,Xt1,Ct1,params1,colinfo1,dt,i+1,0);  
            
            % Correct VW{i} as average of VW1 and VW2 (and associated quantities)
            VW{i}    = 0.5*(VW1+VW2); 
            mu{i}    = 0.5*(mu1+mu2);
            U{i}     = 0.5*(U1+U2);
            Us{i}     = 0.5*(Us1+Us2);
            
            [Xt{i+1},Mt(i+1,:),Ct{i+1},Kernels,params,colinfo,mindst] = LOCAL_advance_step(VW{i},mu{i},Xt{i},X0,Mt(i,:),Ct{i},Kernels,params,colinfo,dt); 
            
            Usqr{i} = params.parsh.U; Usqrh{i} = params.parsh.Uh; 
            gamma{i} = params.parsh.shellden; gammah{i} = params.parsh.Vh; 
            
        case 'rk4'
            fprintf('\n ---------------------------------------------------------- \n')
            fprintf('\n (1) Runge-Kutta 4th order, t=t0, x=x(t0) \n')
            dt41 = 0.5*dt; 
            [Xt1,Mt1,Ct1,U1,Us1,mu1,VW1,Kernels1,params1,colinfo1,dt41] = ...
                LOCAL_euler_step(Xt{i},X0,Mt(i,:),Ct{i},Kernels,params,colinfo,dt41,i,0); 
            
            dt = 2*dt41; 
            
            fprintf('\n ---------------------------------------------------------- \n')
            fprintf('\n (2) Runge-Kutta 4th order, t_{1/2}=t0+0.5*dt, x=x(t_{1/2}) \n')
            dt42 = 0.5*dt; 
            [U2,Us2,mu2,VW2,params] = LOCAL_compute_velocities(Kernels1,Mt1,Xt1,Ct1,params1,colinfo1,i+1,0);
            
            [Xt2,Mt2,Ct2,Kernels2,params2,colinfo2,~] = LOCAL_advance_step(VW2,mu2,Xt{i},X0,Mt(i,:),Ct{i},Kernels,params,colinfo,dt42); 
            
            fprintf('\n ---------------------------------------------------------- \n')
            fprintf('\n (3) Runge-Kutta 4th order, t_{1/2}=t0+0.5*dt, x=x(t_{1/2}) \n') 
            dt43 = dt; 
            
            [U3,Us3,mu3,VW3,params] = LOCAL_compute_velocities(Kernels2,Mt2,Xt2,Ct2,params2,colinfo2,i+1,0);
            
            [Xt3,Mt3,Ct3,Kernels3,params3,colinfo3,~] = LOCAL_advance_step(VW3,mu3,Xt{i},X0,Mt(i,:),Ct{i},Kernels,params,colinfo,dt43); 
            
            fprintf('\n ---------------------------------------------------------- \n')
            fprintf('\n (4) Runge-Kutta 4th order, t_{1}=t0+dt, x=x(t_{1}) \n')  
            
            [U4,Us4,mu4,VW4,params] = LOCAL_compute_velocities(Kernels3,Mt3,Xt3,Ct3,params3,colinfo3,i+1,0);
            
            % Average velocities and related quantities using RK4 weights
            VW{i}    = (1/6)*(VW1+2*VW2+2*VW3+VW4); 
            mu{i}    = (1/6)*(mu1+2*mu2+2*mu3+mu4); 
            U{i}     = (1/6)*(U1+2*U2+2*U3+U4); 
            Us{i}     = (1/6)*(Us1+2*Us2+2*Us3+Us4); 
            
            [Xt{i+1},Mt(i+1,:),Ct{i+1},Kernels,params,colinfo,mindst] = LOCAL_advance_step(VW{i},mu{i},Xt{i},X0,Mt(i,:),Ct{i},Kernels,params,colinfo,dt); 
            
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = t+dt;
tt(i+1)=t;
fprintf('\n Time: %2.2f ',t)

if params.parsq.beta > 0 
    fname = [name '_squirmshell_' num2str(shrd) 'ns_' num2str(ns) '_puller' num2str(params.parsq.beta) ...
        '_p' num2str(p) '_' tdisc '_dt' num2str(dt0) '_mdist' num2str(params.parsq.mdist) '.mat'];
elseif params.parsq.beta < 0 
    fname = [name '_squirmshell_' num2str(shrd) 'ns_' num2str(ns) '_pusher' num2str(abs(params.parsq.beta)) ...
        '_p' num2str(p) '_' tdisc '_dt' num2str(dt0) '_mdist' num2str(params.parsq.mdist) '.mat'];
else
    fname = [name '_squirmshell_' num2str(shrd) 'ns_' num2str(ns) '_treadmill' ... 
        '_p' num2str(p) '_' tdisc '_dt' num2str(dt0) '_mdist' num2str(params.parsq.mdist) '.mat'];
end 

%['test2_squirmers_' num2str(ns) '_puller' num2str(params.parsq.beta) ...
      %  '_p' num2str(p) '_' tdisc '_dt' num2str(dt0) '_mdist' num2str(params.parsq.mdist) '.mat'];
                                                                                                                                             
save(fname,'-v7.3','tt','Xt','Mt','Ct','mu','U','Us','VW','Usqr','Usqrh','gamma','gammah','Conc');                                                                                         
  
end

end

function params = LOCAL_setparams(p,C,r,flag_pot,Shape,Mt,dense,doAna,mdist,out,kerd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%kernel dimension (3 or 1) 
if(nargin < 11)
kerd = 3; 
end
% Points on each swimmer
np = 2*p*(p+1); 
% Number of rigid bodies
nc=size(C,1); 
% Model Surfaces
Sc = SurfaceSph(shape_gallery(p,Shape)); 
% Model surface pts
X = reshape(Sc.cart.to_array,[],3);
Cg = reshape(repmat(C.',np,1),3,[]).';
rg = repmat(reshape(repmat(r.',np,1),1,[]).',1,3);
Nr = reshape(Sc.geoProp.nor.to_array,[],3);
Xrp = rg.*repmat(X,nc,1); 

Xg = Xrp + Cg;
Nrg = repmat(Nr,nc,1); 

% Smooth Quadrature Weights (GL x Trapezoidal)
[~, gwt]=g_grid(p+1);
wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
wt = wt(:);
% Area element
W = Sc.geoProp.W; 
W = W.*wt; 
Wg = (rg(:,1).^2).*repmat(W,nc,1); 

% Repeat kerd times for vector kernels
Xv = reshape(repmat(Xg,1,kerd)',3,[])'; 
Nrv = reshape(repmat(Nrg,1,kerd)',3,[])';    
Wv = repmat(Wg,1,kerd)'; Wv = Wv(:);    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Form and fill params struct
params = struct('flag_pot',flag_pot,'kh',0,'transinv',1,'sym',1,'proxy',0,...
    'keval','points','dim',3,'mu',1,'Xrp',Xrp,'Xp',Xg,'Nrp',Nrg,'X',Xv,'nor',Nrv,'Wg',Wg,'W2',Wv.',...
    'p',p,'np',np,'kerd',kerd,'n3',nc,'rd',r(:),'dense',dense,'C',C,'out',out,...
    'mdist',mdist,'a',0,'doAna',doAna,'Sc',Sc); 

if kerd>1
params.ci = repmat((1:kerd)',nc*np,1); params.cj = params.ci; 
end

end

%{
function Us = LOCAL_Uslip(B1,B2,u,v,p,np,ns,M)

%Us = B1*sin(u) + B2*sin(u).*cos(u);
dir = [cos(u).*cos(v) cos(u).*sin(v) -sin(u)];
Us=zeros(size(dir)); 

[u0,~] = gl_grid(p); u0=u0(:); 

for k=1:ns
   ind = (np*(k-1)+1):np*k; 
   dir(ind,:) = dir(ind,:)*M{k}';
   %u0 = u(ind); 
   Us(ind,:) = (B1*sin(u0) + B2*sin(u0).*cos(u0)).*dir(ind,:); 
end 

end
%}
function [Us, Concentration] = LOCAL_Solve_Neumann(params,Kernels)

% solves Neumann Problem
ns = size(params.parsq.C,1);
np = params.parsq.np; %size(params.parsq.Xrp,1)/ns;
init_dir = repmat([0 0 1], ns,1);
surfacelabel = zeros(size(params.parsq.Xrp,1),1); 
mobilitylabel = surfacelabel;
%Get angles
for i = 1:ns
       Xsphere=params.parsq.Xrp(np*(i-1)+1:np*i,:);
       Xbackrot=Xsphere*params.parsq.Mt{i}; 
       surfacelabel(np*(i-1)+1:np*i,:)=params.A_label(acos(Xbackrot*init_dir(i,:)'));
       mobilitylabel(np*(i-1)+1:np*i,:) = params.M_label(acos(Xbackrot*init_dir(i,:)'));
end

%Solve System

if(isfield(params,'parLapsh'))
    npshell = size(params.parLapsh.X,1); 
 %   densityAug = Lslv(Kernels.dSL, [zeros(npshell,1); -surfacelabel]);
    if(isfield(Kernels.dSL,'B_schur'))
    densityp = Lslv(Kernels.dSL.B_schur,Kernels.dSL.RHS_schur([params.parLapsh.shell_flux.*ones(npshell,1); -surfacelabel]));
    densitys = Kernels.dSL.Shell_Solve([params.parLapsh.shell_flux.*ones(npshell,1); -surfacelabel],densityp); 
    densityAug = [densitys;densityp];
    else
    densityAug = Lslv(Kernels.dSL, [params.parLapsh.shell_flux.*ones(npshell,1); -surfacelabel]);
    densitys = densityAug(1:npshell);
    densityp = densityAug(npshell+1:end);
    
    end
    fprintf('Prescribed flux from Particle= %1.2e \n',params.parLap.W2*(-surfacelabel));
    fprintf('Prescribed flux from Shell =%1.2e \n',sum(params.parLapsh.W2*params.parLapsh.shell_flux));
    
    %computes concentration on particles using density from shell, particle
    params.parLap.a=0; 
    Concentration = RBS_MatVec(densityp,[],'Vsh',params.parLap,1,0,'SL_L_3D',[]);
    a_hold = params.parLapsh.a; params.parLapsh.a = 0; params.parLap.flag_pot='SL_L_3D';
    Shell_Contrib = VSh_MatVec_RB_trg(densitys,params.parLap.X,params.parLap.Nrp,params.parLapsh);
    Concentration = Concentration + Shell_Contrib;
    
    params.parLapsh.a = a_hold;
    

    
    %compute Shell Concentration using density from shell, particle

    Particle_Contrib = VSh_MatVec_RB_trg(densityp,params.parLapsh.X,params.parLapsh.Nrp,params.parLap);
    Shell_Conc_Mat = RBS_MatVec([],[],'Vsh',params.parLapsh,1,0,'SL_L_3D',[]);
    Shell_Concentration = Shell_Conc_Mat(densitys) + Particle_Contrib;
else
    density = Lslv(Kernels.dSL, -surfacelabel); params.parLap.a=0; 
    Concentration = RBS_MatVec(density,[],'Vsh',params.parLap,1, 0,'SL_L_3D',[]);
end
    params.parLap.Concentration = Concentration;
%evaluate concentration on surfaces
GradC = zeros(size(Concentration,1),3); Us=GradC; 
for i = 1:ns
    indx = np*(i-1)+1:np*i;
    S = params.parsq.Sc;
    
    %G = S.geoProp.Grad(Concentration(indx));
    %GradC(indx,1) = G.x; GradC(indx,2) = G.y; GradC(indx,3) = G.z;  
    %GradC = real(GradC);
    
    GradC(indx,:) = reshape(Sph_spdiff(Concentration(indx),1),[],3); 
    
    Gdotn = sum((params.parsq.Nrp(indx,:)).*GradC(indx,:),2);
    %ProjGradC = Gdotn.*(params.parsq.Nrp(indx,:));
    %fprintf('\n |GradC*n| = %1.2e',norm(Gdotn));  
    Us(indx,:) = mobilitylabel(indx,:).*GradC(indx,:); %-ProjGradC);

end

if(isfield(params,'parLapsh'))
   Concentration = [Concentration;Shell_Concentration]; 
end

end

function [Uh,params] = LOCAL_Uinf_Vsph(p,flow,params)

shf = @(q) VshAna([q(1:3:end);q(2:3:end);q(3:3:end)],'VW'); 

switch flow 
    case 'stokeslets'
        parS = params; 
        parS.flag_pot = 'SL_Stk_3D'; 
        C = [-1.5 0 0;2 0 0;0 -1.5 0;0 2 0;0 0 -1.5;0 0 2]; % + 0.25*(2*rand(6,3)-1); 
        display(C)
        Cg = reshape(repmat(C,1,3)',3,[])';
        Xv = parS.X;   
        n3 = size(C,1); 
        parS.W2 = ones(1,3*n3); 
        parS.cj = repmat((1:3)',n3,1);
        alpha = ones(3*n3,1); 
        params.C2 = C; 
        params.alpha = alpha; 
        
        U = Kernel_Eval(Xv,Cg,parS)*alpha; 
        params.U = U; 
    case 'synthetic'
        fun = @(n) 2.^(-n); 
        sp = (p+1)^2; 
        ii = (1:sp)'; 
        nn = floor(sqrt(ii-1)); 
        pmx = max(p-1,4); 
        qVh = [rand(pmx^2,1); zeros((p+1)^2-pmx^2,1)];
        qWh = [rand(pmx^2,1); zeros((p+1)^2-pmx^2,1)];
        qXh = [rand(pmx^2,1); zeros((p+1)^2-pmx^2,1)];
        % normalize
        qVh = fun(nn).*(qVh); 
        qWh = fun(nn).*(qWh);
        qXh = fun(nn).*(qXh);
        
        %qVh = zeros(size(qVh)); qWh=qVh; qXh = qVh; qWh(4)=1; 
        
        qh = [qVh;qWh;qXh]; 
        Q = VshSyn(qh,'VW'); 
        Q = reshape(reshape(Q,[],3).',[],1);  
        U = VSh_MatVec_RB2(Q,[],params);
        
        params.Vh = qh; 
        params.Q = Q; 
        params.U = U; 
end

Uh = shf(U);
%plot(abs(Uh),'or'); hold off; 

end

function Utrg = LOCAL_Uinf_trg(Xtrg,Nrtrg,flow,params)

switch flow 
    case 'stokeslets'
        params.flag_pot = 'SL_Stk_3D'; 
        %C = [-1.5 0 0;2 0 0;0 -1.5 0;0 2 0;0 0 -1.5;0 0 2];
        C = params.C2; 
        Cg = reshape(repmat(C,1,3)',3,[])';
        
        Xtrgv = reshape(repmat(Xtrg,1,3)',3,[])';  
        params.ci = repmat((1:3)',size(Xtrg,1),1);
        
        n3 = size(C,1); 
        params.W2 = ones(1,3*n3); 
        params.cj = repmat((1:3)',n3,1);
        params.nor = Nrtrg;
        alpha = params.alpha;
        
        Utrg = Kernel_Eval(Xtrgv,Cg,params)*alpha;  
    case 'synthetic'
        params.flag_pot = 'SL_Stk_3D';
        Q = params.Q; 
        [Utrg,~,~] = VSh_MatVec_RB_trg(Q,Xtrg,Nrtrg,params);
end

end

function [FV,FW,FX] = LOCAL_get_eigen(type,out)
% Eigenvalues
%SVeg  = @(n) n./((2*n+1).*(2*n+3)); 
%SWeg  = @(n) (n+1)./((2*n+1).*(2*n-1));
%SXeg  = @(n) 1./(2*n+1);
if out
switch type 
    case 'SMat'
    % f(r,n) exterior (SMat) 
    FV = @(n) n./((2*n+1).*(2*n+3)); %SVeg; 
    FW = @(n) (n+1)./((2*n+1).*(2*n-1)); %SWeg; 
    FX = @(n) 1./(2*n+1); %SXeg; 
    case 'SpMat'
    % f'(r,n) exterior (SpMat)
    SpVext  = @(n) -(n+2).*(n./((2*n+1).*(2*n+3))); %-(n+2)*SVeg(n); 
    SpWVext = @(n) (-2*n./(4*n+2));
    SpWWext = @(n) -n.*((n+1)./((2*n+1).*(2*n-1))); %-n*SWeg(n);
    SpXext  = @(n) -(n+1)./(2*n+1); %-(n+1)*SXeg(n);
    PWext = @(n) -n; 

    FV = SpVext; 
    FW = {SpWWext,SpWVext,PWext}; 
    FX = SpXext;

    case 'TMat'
    % f(r,n) exterior (Traction)
    %TV
    TVext = @(n) 2*(-n-2).*(n./((2*n+1).*(2*n+3))); %-2(n+2)*SVeg(n);
    %TW 
    TWext = @(n) ((1+2*n.^2)./(1-4*n.^2)); 
    %TX
    TXext = @(n) -(n+2)./(2*n+1); %-(n+2)*SXeg(n); 

    FV = TVext; 
    FW = TWext; 
    FX = TXext; 
    case 'DMat'
    % f(r,n) exterior (double layer)
    %DV
    DVext = @(n) ((2*n.^2+4*n+3)./((2*n+1).*(2*n+3)));
    %DW
    DWext = @(n) ((2*n.^2-2)./(4*n.^2-1)); 
    %DX
    DXext = @(n) ((n-1)./(2*n+1)); 

    FV = DVext; 
    FW = DWext; 
    FX = DXext;  
    case 'DpMat'
    % f(r,n) exterior (double layer)
    %DV
    DVext = @(n) (-n-2).*((2*n.^2+4*n+3)./((2*n+1).*(2*n+3)));
    %DW
    DWVext = @(n) ((-2*n.*(n-1))./(2*n+1)); 
    DWWext = @(n) -n.*((2*n.^2-2)./(4*n.^2-1)); 
      
    %DX
    DXext = @(n) (-n-1).*((n-1)./(2*n+1)); 
    
    %Pressure 
    PWext = @(n) -n.*(n+1); 

    FV = DVext; 
    FW = {DWWext,DWVext,PWext}; 
    FX = DXext;  
    case 'SDMat'
    % f(r,n) exterior (single + double layer)
    %SDV
    SDVext = @(n) (n+1)./(2*n+1);
    %SDW
    SDWext = @(n) (n+1)./(2*n+1); 
      
    %SDX
    SDXext = @(n) n./(2*n+1); 

    FV = SDVext; 
    FW = SDWext; 
    FX = SDXext;  
end
else

switch type 
    case 'SMat'
    % f(r,n) exterior (SMat) 
    FV = @(n) n./((2*n+1).*(2*n+3)); %SVeg; 
    FW = @(n) (n+1)./((2*n+1).*(2*n-1)); %SWeg; 
    FX = @(n) 1./(2*n+1); %SXeg; 
    case 'SpMat'
    % f'(r,n) interior (SpMat)
    SpVVint = @(n) (n+1).*(n./((2*n+1).*(2*n+3))); %)SVeg(n).*r.^n; 
    SpVWint = @(n) (2*(n+1)./(4*n+2));
    SpWint  = @(n) (n-1).*((n+1)./((2*n+1).*(2*n-1))); %(n-1)*SWeg(n).*r.^(n-2);
    SpXint  = @(n) n./(2*n+1); %n*SXeg(n);
    PVint = @(n) -(n+1);

    FV = {SpVVint,SpVWint,PVint}; 
    FW = SpWint; 
    FX = SpXint;

    case 'TSMat'
    % f(r,n) interior (Traction of SL)
    %TV
    %Simplified these two as f(n)r^n+g(n)r^(n-2)  
    TVint = @(n) ((3+4*n+2*n.^2)./(3+8*n+4*n.^2)); 
    %TW
    TWint = @(n) 2*(n-1).*((n+1)./((2*n+1).*(2*n-1))); %SWeg(n); 
    %TX
    TXint = @(n) (n-1)./(2*n+1); %*SXeg(n);

    FV = TVint; 
    FW = TWint; 
    FX = TXint;
    case 'DMat'
    % f(r,n) interior (Double Layer)
    %DV 
    DVint = @(n) ((-2*n.*(n+2))./((2*n+1).*(2*n+3))); 
    %TW
    DWint = @(n) -((2*n.^2+1)./((2*n+1).*(2*n-1))); 
    %TX
    DXint = @(n) -((n+2)./(2*n+1));

    FV = DVint; 
    FW = DWint; 
    FX = DXint;    
    case 'DpMat'
    % f(r,n) interior (normal derivative Double Layer)
    %DV 
    DVVint = @(n) (n+1).*((-2*n.*(n+2))./((2*n+1).*(2*n+3))); 
    DVWint = @(n) -((2*(n+1).*(n+2))./(2*n+1)); 
    %TW
    DWint = @(n) -(n-1).*((2*n.^2+1)./((2*n+1).*(2*n-1))); 
    %TX
    DXint = @(n) -n.*((n+2)./(2*n+1));
    
    PVint = @(n) n.*(n+1);

    FV = {DVVint,DVWint,PVint}; 
    FW = DWint; 
    FX = DXint; 
    case 'SDMat'
    % f(r,n) interior (SL + DL)
    %[S+D]V 
    SDVint = @(n) -n./(2*n+1);  
    %[S+D]W
    SDWint = @(n) -n./(2*n+1); 
    %[S+D]X
    SDXint = @(n) -(n+1)./(2*n+1);

    FV = SDVint; 
    FW = SDWint; 
    FX = SDXint;
end
end

end

function y = Lapp(A,x)
if isnumeric(A)
    y=A*x; 
else
    y=real(A(x)); 
end
end

function x = Lslv(A,b,x0,restart,tol,maxit,pr)
global prec;

if nargin<6
   if ~isempty(prec) 
       pr=prec;  
   else
       pr=[]; 
   end
end

if isnumeric(A)
    x=A\b;
    %{
    [U,S,V] = svd(A); 
    ss = diag(S); 
    ks = sum(ss>(1e-12)*ss(1)); 
    ssi = zeros(size(ss)); 
    ssi(1:ks)=1./ss(1:ks);
    x = V*(diag(ssi)*(U'*b));
    %}
else
    x = zeros(size(b)); 
    for i=1:size(b,2)
    if nargin==2
        x(:,i)=gmres(A,b(:,i),8,1e-5,100,pr);
    elseif nargin==3
        x(:,i)=gmres(A,b(:,i),8,1e-5,100,pr,[],x0);
        %x = x + x0; 
    else
        x(:,i)=gmres(A,b(:,i),restart,tol,maxit,pr);
    end
    end
end
end

function [U,Us,mu,VW,params,dt,Concentration] = LOCAL_compute_velocities(Kernels,Mt,Xt,Ct,params,colinfo,dt,it,flplot)

denseS = params.parsq.dense; 
denseL = params.parLap.dense; 
shflg = isfield(params,'parsh');
colevent=colinfo.col; 

if shflg
    parsh = params.parsh; 
    shrd=parsh.rd; 
end
parsq = params.parsq; 

p = parsq.p; np = parsq.np; ns = parsq.n3; 

[Us0,Concentration] = LOCAL_Solve_Neumann(params,Kernels);
Us0 = reshape(Us0.',[],1);
Us = Us0; 

Nmu = 3*np*ns; Nvw = 6*ns;
rhs = [Us0 ; zeros(Nvw,1)];

if shflg
    if ~isfield(params,'Y0')
        Y0 = zeros(size(rhs)); 
    else
        Y0 = params.Y0; 
    end
    
    [U,mu,VW,params] = LOCAL_solve_joint_system_dir(rhs,Kernels,params,Y0,colinfo); 
    
    params.Y0 = [mu;VW(:)]; 
    Ushell=params.parsh.Ushell; 
    Us0 = Us0 - Ushell; % + Urb;
    rhs = [Us0 ; zeros(Nvw,1)];
    Err = abs(Lapp(Kernels.M,[mu;VW(:)])-rhs);  
else
    Y = Lslv(Kernels.M,rhs); 
    VW = real(Y(Nmu+1:Nmu+Nvw)); 
    Err = abs(Lapp(Kernels.M,Y)-rhs);
    
    % separate density mu from rigid body velocities VW  
    mu = Y(1:Nmu); 
    U = Lapp(Kernels.ACF,mu);   
end 

VW = reshape(VW,6,[]); 

mxV = max(sqrt(sum(VW(1:3,:).^2)));
if dt*mxV>0.1
   fprintf('\n Max velocity exceeded, max(V) = %2.2f, new dt=%1.2e \n',mxV,0.1/mxV)
   dt = 0.1/mxV;  
end

if colevent
    display([colinfo.ip colinfo.jp])
    ncol = sum(colinfo.jp<ns+1); %number of particle collisions
    
    % Setup LCP for contact forces
    ip = colinfo.ip(1:ncol); jp = colinfo.jp(1:ncol); 
    
    Em = zeros(ncol,Nvw);
    % Compute vectors and normal vectors for pairs
    R = Ct(ip,:)-Ct(jp,:);  %Ci - Cj numF x 3
    NR = sqrt(sum(R.*R,2)); %|Ci-Cj| numF x 1
    RSIJ = parsq.rd(ip)+parsq.rd(jp);  
    Rhat = repmat(1./NR,1,3).*R; %eij = (Ci - Cj)/|Ci-Cj|
        
    for k=1:ncol
        %Rhat = Ct(ip(k),:)-Ct(jp(k),:); Rhat=Rhat./norm(Rhat); 
        indk = [(1:3)+6*(ip(k)-1) (1:3)+6*(jp(k)-1)];   
        Em(k,indk) = [Rhat(k,:) -Rhat(k,:)]; 
    end
       
    Fm = Em.';
    phib = (1/dt)*(NR-RSIJ-colinfo.eps); 
    bvec = phib + real(Em*VW(:)); 
    
    if shflg
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Add collisions with boundary shell (if they exist)
        ids=colinfo.jp>ns; 
        ncols=sum(ids);
        
        display(ncols)
        if ncols>0
           ncolT = ncol+ncols; 
           ips = colinfo.ip(ids); %jps = colinfo.jp(ids); 
           % Compute normal vectors for shell collisions
           R = Ct(ips,:);  %Ci - Cshell (assumed [0 0 0]) numF x 3
           NR = sqrt(sum(R.*R,2)); %|Ci| numF x 1
           RSI = parsq.rd(ips);  
           Rhat = -repmat(1./NR,1,3).*R; %eij = -Ci/|Ci|
           
           if size(Em,1)==0
               Em = zeros(ncolT,Nvw); 
           end
           
           for k=ncol+1:ncolT
               indk = (1:3)+6*(ips(k-ncol)-1);   
               Em(k,indk) = Rhat(k-ncol,:);
           end
           
           Fm = Em.';
           phib = [phib; (1/dt)*(shrd-NR-RSI-colinfo.epsh*shrd)]; 
           bvec = phib + real(Em*VW(:));
           ncol = ncolT; 
        end 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isfield(params,'Y0c')
           Y0c = zeros(size(rhs)); 
           Y0c(Nmu+1:end) = -VW(:); 
        else
           Y0c = params.Y0c;  
        end
        
        FTlam = @(lam) Fm*lam; 
        rhs_c = @(lam) [zeros(Nmu,1); Fm*lam];
        VWmob = @(lam) getVW_joint_system(rhs_c(lam),Kernels,params,Y0c);
        Amat = @(lam) Em*VWmob(lam); 
        
        %LCP params
        max_iter=100; tol_rel=1e-5; tol_abs=1e-6; profile=1;
        % solve LCP
        [lam ,err ,iter, ~, ~, ~] = ...
    minmap_newton_matfree(Amat, bvec, zeros(size(bvec)), max_iter, tol_rel, tol_abs, profile );
    %BBPGD_matfree(Amat, bvec, zeros(size(bvec)), max_iter, tol_rel, tol_abs, profile );               
        
        fprintf('\n Minmap Newton Matfree LCP solution error = %e, iters = %d',err,iter); 
      
        % Compute FT, update mu and VW 
        FT = FTlam(lam); 
        [U_c,mu_c,VW_c,params] = LOCAL_solve_joint_system_dir(rhs_c(lam),Kernels,params,Y0c,colinfo);
        
        params.Y0c = [mu_c; VW_c(:)]; 
        
        % Update
        U = U + U_c; 
        mu = mu + mu_c;
        VW_c = reshape(VW_c,6,[]); 
        VW = VW + VW_c; 
        
        FT = reshape(FT,6,[]); 
        fprintf('\n Contact Forces');
        display(FT) 
    else    
        %LCP params
        max_iter=100; tol_rel=1e-6; tol_abs=1e-12; profile=1;
        lam0 = zeros(size(bvec));
        
        if denseS 
            %Build Amat densely from Schur complement  
            AD = real(Lslv(Kernels.ACF,Kernels.D.')); 
            Schr = Kernels.C*AD; 
            Amat = Em*(Schr\Fm);
            % solve LCP
            [lam ,err ,iter, ~, ~, ~] = ...
                minmap_newton(Amat, bvec, lam0, max_iter, tol_rel, tol_abs, profile );
            
            % Compute FT, update mu and VW 
            FT = Fm*lam; 
            VW_c = Schr\FT;
        else
            %Build Amat as function handle, use matfree versions
            FTlam = @(lam) Fm*lam;
            rhs_c = @(lam) [zeros(Nmu,1); Fm*lam];
            getVW = @(Y) Y(Nmu+1:end,:); 
            Amat = @(lam) Em*getVW(Lslv(Kernels.M,rhs_c(lam)));
            
            % solve LCP (matfree)
            [lam ,err ,iter, ~, ~, ~] = ...
                minmap_newton_matfree(Amat, bvec, zeros(size(bvec)), max_iter, tol_rel, tol_abs, profile );
                %BBPGD_matfree(Amat, bvec, zeros(size(bvec)), max_iter, tol_rel, tol_abs, profile );
            
            % Compute FT, update mu and VW 
            FT = FTlam(lam); 
            VW_c = getVW(Lslv(Kernels.M,rhs_c(lam)));    
        end
            
        fprintf('\n Minmap Newton LCP solution error = %e, iters = %d',err,iter);  
        
        % Update
        mu = mu + AD*VW_c;
        VW_c = reshape(VW_c,6,[]); 
        VW = VW + VW_c; 
        
        FT = reshape(FT,6,[]); 
        fprintf('\n Contact Forces');
        display(FT)
        
        U = Lapp(Kernels.ACF,mu); 
    end
end

fprintf('\n Rigid body velocities');
VW = reshape(VW,6,[]); 
display(VW)

fprintf('\n Full residual = %0.5g',norm(Err));

% plotb 
stride=1; 
if flplot && mod(it-1,stride)==0 
phi = abs(mu);
Cg = reshape(repmat(parsq.C.',np,1),3,[]).';
LOCAL_plotspheres(Xt+Cg,p,parsq.C,parsq.rd,phi)
pause(0.00001); 
end
end

function VW = getVW_joint_system(rhs,Kernels,params,Y0c)

[~,~,VW]=LOCAL_solve_joint_system_dir(rhs,Kernels,params,Y0c,[]);

end

function params = LOCAL_compute_shell_velocity(mu,Kernels,params,verb)

if nargin<4
    verb=1; 
end

nq=3*params.parsh.np; 
J = [1:3:nq 2:3:nq 3:3:nq]; 
shf = @(q) VshAna(q(J,:),'VW'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verb
    fprintf('\n Evaluation of squirmer flow on shell and boundary correction: \n');
end

Xtrg = params.parsh.Xp; Nrtrg = params.parsh.Nrp; 
% Evaluate DL + Completion from particles at boundary
partrg = params.parsq; partrg.a=0; 
Usqr = real(VSh_MatVec_RB_trg(mu,Xtrg,Nrtrg,partrg));

%Add completion flow
Cm = Kernels.C; 
Xtrgv = reshape(repmat(Xtrg,1,3)',3,[])';

ns=params.parsq.n3; 
parStk = params.parsq; parStk.flag_pot='SL_Stk_3D'; parStk.a=0; 
parStk.ci=repmat((1:3)',size(Xtrg,1),1);
parStk.cj=repmat((1:3)',ns,1); 
parStk.W2 = ones(1,3*ns); 

Cv = reshape(repmat(params.parsq.C,1,3)',3,[])';  
Gij = Kernel_Eval(Xtrgv,Cv,parStk); %Stokeslet
parRot = parStk; parRot.flag_pot='Rotlet_3D'; 
Rij = Kernel_Eval(Xtrgv,Cv,parRot); %Rotlet

Km = zeros(size(Gij,1),2*size(Gij,2)); 
for k=1:ns
    Km(:,(1:6)+6*(k-1)) = [Gij(:,(1:3)+3*(k-1)) Rij(:,(1:3)+3*(k-1))]; 
end
Usqr = Usqr + Km*(Cm*mu); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Uh = shf(Usqr); 
sigmah = -params.parsh.eigI.*Uh; 
sigmah(1) = -Uh(1); 
sigma  = real(VshSyn(sigmah,'VW'));
sigma = reshape(reshape(sigma,[],3).',[],1);

% Set params for shell density apply
params.parsh.U = Usqr; 
params.parsh.Uh = Uh; 
params.parsh.shellden = sigma; 
params.parsh.Vh = sigmah; 

V00 = (-1/sqrt(4*pi))*reshape(Nrtrg.',[],1); 

Ucorr = VSh_MatVec_RB2(sigma,[],params.parsh) - Uh(1)*V00; 
Uinf = Usqr + Ucorr;
params.parsh.Uinf = Uinf; 

if verb
    fprintf('\n Check flow at boundary = 0: ||Uinf||_2 = %0.5g , ||Uinf||_inf = %0.5g',norm(Uinf),max(abs(Uinf))); 
end 

end

function [Y,params] = LOCAL_apply_joint_system_dir(mu,VW,Kernels,params,colinfo)

% Parameters 
np = params.parsq.np; ns = params.parsq.n3;
Nsq = 3*np*ns; %Nvw = 6*ns;

% Evaluation at shell
params = LOCAL_compute_shell_velocity(mu,Kernels,params,0); 

% Evaluation at squirmer surfaces
Xtrg = params.parsq.Xp; Nrtrg = params.parsq.Nrp; 
        
% Evaluate V00(u,v) at particle discretization pts Xtrg
[thg,phig,~] = cart2sph(Xtrg(:,1),Xtrg(:,2),Xtrg(:,3));  
thg(thg<0)=thg(thg<0)+2*pi; phig=pi/2-phig; 
thg=thg(:); phig=phig(:);  
V00 = -sqrt(1/(4*pi))*[sin(phig).*cos(thg) sin(phig).*sin(thg) cos(phig)];
V00 = reshape(V00.',[],1); 
        
Y = Lapp(Kernels.M,[mu;VW]);

% Correction Ushell = (D+N_inf)[mu] = D[mu] - Uh(1)*V_0^0
Ushell = real(VSh_MatVec_RB_trg(params.parsh.shellden,Xtrg,Nrtrg,params.parsh)) - params.parsh.Uh(1)*V00;
params.parsh.Ushell=Ushell; 

Y(1:Nsq) = Y(1:Nsq)+Ushell;

end

function [U,mu,VW,params] = LOCAL_solve_joint_system_dir(rhs,Kernels,params,y0,colinfo)

np = params.parsq.np; ns = params.parsq.n3; 
Nsq = 3*np*ns; 

Mapp = @(V) LOCAL_apply_joint_system_dir(V(1:Nsq),V(Nsq+1:end),Kernels,params,colinfo);  
Mprec = @(V) Lslv(Kernels.M,V,zeros(size(V)),20,1e-5,5,[]); 

%gmres solve of joint system 
restart=30; tol=1e-6; maxit=5;

tic; 
Sol = Lslv(Mapp,rhs,y0,restart,tol,maxit,[]);
fprintf('\n Joint system solve time = %1.2e',toc);  

% Separate density and rigid-body velocities 
mu = Sol(1:Nsq); 
VW = Sol(Nsq+1:end); 

[USol,params] = LOCAL_apply_joint_system_dir(mu,VW,Kernels,params,colinfo); 
U = USol(1:Nsq); 

end

function B = LOCAL_compute_BIE_matrix(params,shell)

    shflg = isfield(params,'parLapsh');
    ns = size(params.parLap.C,1);
    np = 2*(params.parLap.p+1)*params.parLap.p;
    NE = np*ns;
    typeMV = 'Vsh';
    params.parLap.out = 1;
    PM = make_projection(params.parLap.p);

    if(shflg)
        NI = 2*(params.parLapsh.p+1)*params.parLapsh.p;
          if(params.parLap.dense==1)
            Bp2p= RBS_MatVec([],[],typeMV,params.parLap,1,-0.5,'dSL_L_3D',[]);
            PME = kron(eye(ns),PM); %block-diag of projections for exterior cases
            Bp2p = Bp2p + params.parLap.a*(eye(NE)-PME);
            PMs = make_projection(params.parLapsh.p); 
            Bs2s = RBS_MatVec('Mat',[],typeMV,params.parLapsh,1,0.5,'dSL_L_3D',[]);
            Bs2s = Bs2s + params.parLapsh.a*(eye(NI)-PMs);
            
            %one dimensional kernel in shell-to-shell matrix removed view
            %completion flow
            W = params.parLapsh.W2; W = W(:);
            Y00 = (1/sqrt(4*pi))*ones(NI,1); 
            Q_orth = (Y00)*((W/params.parLapsh.rd^2).*Y00)'; 
            Bs2s = Bs2s + Q_orth;

           
           %make shell-to-particle interaction
           a_hold = params.parLapsh.a;
           Bs2p = VSh_MatVec_RB_trg('Mat',params.parLap.X,params.parLap.Nrp,params.parLapsh);
           params.parLapsh.a=a_hold;
           a_hold = params.parLap.a; params.parLap.a=0;
           %make particle-to-shell interaction
           Bp2s = VSh_MatVec_RB_trg('Mat',params.parLapsh.X,params.parLapsh.Nrp,params.parLap);
           
           %augmented system
           B = [Bs2s Bp2s; Bs2p Bp2p];

          elseif(params.parLap.dense == 0)
          
                Projhandle = @(V) params.parLap.a*(V - reshape(PM*(reshape(V,[],ns)),[],1));
                Bp2p= @(V) RBS_MatVec(V,[],typeMV,params.parLap,1,-0.5,'dSL_L_3D',[]) + Projhandle(V);
                Dinv = @(V) Local_invert_Shell(V,params.parLapsh.Spectra);
                a_hold = params.parLapsh.a; params.parLapsh.a=0;
           
           %make shell-to-particle interaction
                params.parLapsh.doAna = 1;
                Bs2p = @(V) VSh_MatVec_RB_trg(V,params.parLap.X,params.parLap.Nrp,params.parLapsh);
                params.parLapsh.a=a_hold;
                a_hold = params.parLap.a; params.parLap.a=0;
           %make particle-to-shell interaction
                Bp2s = @(V) VSh_MatVec_RB_trg(V,params.parLapsh.X,params.parLapsh.Nrp,params.parLap);    
               
             
           %Schur Complement
           B.B_schur = @(V) Bp2p(V) - Bs2p((Dinv(Bp2s(V))));
           B.RHS_schur = @(V) V(end-NE+1:end) - Bs2p(Dinv(V(1:NI)));
           B.Shell_Solve = @(V,Den) Dinv(V(1:NI) - Bp2s(Den)); 
          end
    else
           if(params.parLap.dense==1)
                Bp2p= RBS_MatVec([],[],typeMV,params.parLap,1,-0.5,'dSL_L_3D',[]);
                PME = kron(eye(ns),PM); %block-diag of projections for exterior cases
                Bp2p = Bp2p + params.parLap.a*(eye(NE)-PME);
           else
                Projhandle = @(V) params.parLap.a*(V - reshape(PM*(reshape(V,[],ns)),[],1));
                Bp2p= @(V) RBS_MatVec(V,[],typeMV,params.parLap,1,-0.5,'dSL_L_3D',[]) + Projhandle(V);
           end
           B = Bp2p; 
     
    end

end


function[V] = Local_invert_Shell(U,Spectra)
%Computes V = (A + QY00)_inv*U given vector U, Spectra of A.
%note: Assumes that QY00 is such that QY00*Y00 = Y00
           ShU = shAna(U);
           ShUEig = ShU./(Spectra + (Spectra==0));
           ShUEig(1) = ShU(1);
           V = shSyn(ShUEig);


end

function [Kernels,params] = LOCAL_update_operators(Ct,Mt,params,Kernels)

denseS = params.parsq.dense; denseL = params.parLap.dense; 
B1 = params.parsq.B1; B2 = params.parsq.B2; beta = params.parsq.beta;  
ns = params.parsq.n3; a = params.parsq.a;  

%Update params
params.parsq = LOCAL_setparams(params.parsq.p,Ct,params.parsq.rd,params.parsq.flag_pot,'',Mt,...
    params.parsq.dense,params.parsq.doAna,params.parsq.mdist,params.parsq.out);
params.parLap = LOCAL_setparams(params.parLap.p,Ct,params.parLap.rd,params.parLap.flag_pot,'',Mt,...
    params.parLap.dense,params.parLap.doAna,params.parLap.mdist,params.parLap.out,params.parLap.kerd);
params.parsq.B1 = B1; params.parsq.B2 = B2; params.parsq.beta = beta;  
params.parsq.a = a; 
params.parsq.Mt = Mt;
params.parLap.a = -0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create auxiliary matrices 
[Cm,Bm,Dm] = LOCAL_Build_AuxMats(params.parsq.Wg,params.parsq.Xrp,[],params.parsq.np,ns);
if denseS
    Kernels.C = full(Cm); Kernels.D = full(Dm); Kernels.B = full(Bm);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Laplace matrix update 
Kernels.dSL = LOCAL_compute_BIE_matrix(params,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Update operators A, ACF and M 

% Power-Miranda Completion flow 
parStk = params.parsq; parStk.flag_pot='SL_Stk_3D'; parStk.a=0; 
parStk.cj=repmat((1:3)',ns,1); parStk.W2 = ones(1,3*ns); 

Cv = reshape(repmat(params.parsq.C,1,3)',3,[])';  
Gij = Kernel_Eval(parStk.X,Cv,parStk); %Stokeslet
parRot = parStk; parRot.flag_pot='Rotlet_3D'; 
Rij = Kernel_Eval(parRot.X,Cv,parRot); %Rotlet

Km = zeros(size(Gij,1),2*size(Gij,2)); 
for k=1:ns
    Km(:,(1:6)+6*(k-1)) = [Gij(:,(1:3)+3*(k-1)) Rij(:,(1:3)+3*(k-1))]; 
end

%rank 6 PM completion flow
Kernels.Km = Km;  

Nmu = 3*params.parsq.np*ns;
% Stokes A = a*I+S+D matrix and augmented system for PM formulation
if denseS
   A = real(LOCAL_MatVec([],params.parsq)); 
   ACF = A + Km*Cm;
   M = [ACF -Dm.' ; Cm zeros(6*ns,6*ns)]; 
else
   A = @(V) real(LOCAL_MatVec(V,params.parsq)); 
   ACF = @(V) A(V) + Km*(Cm*V);
   M = @(V) [ACF(V(1:Nmu,:))-Dm.'*V(Nmu+1:end,:);Cm*V(1:Nmu,:)];
end

Kernels.SDL = A;
Kernels.ACF = ACF; 
Kernels.M = M;

end

function den = LOCAL_CenterDistance(C)

[Y_g1,  X_g1  ] = meshgrid(C(:,1), C(:,1));
[Y_g2,  X_g2  ] = meshgrid(C(:,2), C(:,2));
[Y_g3,  X_g3  ] = meshgrid(C(:,3), C(:,3));
d1 = (X_g1 - Y_g1); d2 = (X_g2 - Y_g2); d3 = (X_g3 - Y_g3); 
den = sqrt(d1.^2 + d2.^2 + d3.^2);
den = 10000*(den==0)+den; 

end

function [C,B,D,L] = LOCAL_BuildC(Wg,Xg,Xc,np,n3)

N = 3*np*n3; 
IC = []; JC = []; VC = []; 
VD = []; VB = [];  

if nargout==4
IL = []; JL = []; VL = []; 
end

for i=1:n3
    
idx = (1:np)+np*(i-1); 
X = Xg(idx,:);
W = Wg(idx); 

oW = ones(size(W));
sumW = sum(W); 

W2 = zeros(1,3*np); 
W2(1:3:3*np)=W; W2(2:3:3*np)=W; W2(3:3:3*np)=W; 
    
if ~isempty(Xc)
% Center X
X = X - repmat(Xc(i,:),np,1); 
end

indx =(1:3:3*np)+3*np*(i-1); 
indy =(2:3:3*np)+3*np*(i-1); 
indz =(3:3:3*np)+3*np*(i-1);

% C*sigma = [int{sigma} ; int{X \times sigma}]    
IC = [IC; reshape(repmat(6*(i-1)+(1:6),np,1),[],1) ; reshape(repmat(6*(i-1)+(4:6),np,1),[],1)];
JC = [JC; indx'; indy'; indz'; indy'; indz'; indx'; indz'; indx'; indy']; 
VC = [VC; W; W; W ; -W.*X(:,3); -W.*X(:,1); -W.*X(:,2); W.*X(:,2);  W.*X(:,3); W.*X(:,1)];
% First three vectors are just integrals of f for each coordinate
VD = [VD; oW; oW; oW ; -X(:,3); -X(:,1); -X(:,2); X(:,2);  X(:,3); X(:,1)];

tau1 = sum(W.*X(:,3).^2)+sum(W.*X(:,2).^2); 
tau2 = sum(W.*X(:,1).^2)+sum(W.*X(:,3).^2); 
tau3 = sum(W.*X(:,2).^2)+sum(W.*X(:,3).^2); 

VB = [VB; oW./sumW; oW./sum(W); oW./sum(W) ; ...
    -X(:,3)./tau1; -X(:,1)./tau2; -X(:,2)./tau3; ...
     X(:,2)./tau1;  X(:,3)./tau2; X(:,1)./tau3];

% We divide by the integral of ((X-X_c) (x) 1)^2 over Sc. 
idxM = (1:3*np)+3*np*(i-1); 
C = sparse(IC,JC,VC,6*n3,N); 
B = sparse(IC,JC,VB,6*n3,N);

if nargout==4  
   [JJ,II] = meshgrid(idxM,idxM); 
   IL = [IL ; II(:)]; 
   JL = [JL ; JJ(:)]; 
   Lb = B(6*(i-1)+(1:6),idxM)'*C(6*(i-1)+(1:6),idxM);   
   VL = [VL ; Lb(:)]; 
end

end

C = sparse(IC,JC,VC,6*n3,N); 
B = sparse(IC,JC,VB,6*n3,N);
D = sparse(IC,JC,VD,6*n3,N);

if nargout==4
   L = sparse(IL,JL,VL,N,N);  
end

end

function Y = LOCAL_MatVec(V,params,p,C,rd,Shape,dense,doAna,mdist,out)
   
if ~isstruct(params) && nargin>2
    params = LOCAL_setparams(p,C,rd,params,Shape,dense,doAna,mdist,out); 
else
    p=params.p; 
    C=params.C;
end

n3 = size(C,1); 

% Build dense matrix
if isempty(V) && params.dense==1
    V = 'Mat'; 
end

a = params.a; 
Y = VSh_MatVec_RB2(V,[],params);
    
np = 2*p*(p+1); 
VProj = @(Y) reshape(VshProj(reshape(Y,3*np,[]),'VW'),3*np*n3,[]); 
    
if isnumeric(V)
    Y = VProj(Y)+a*(V-VProj(V)); 
elseif strcmp(V,'Mat')
    P = make_vshprojection(p,1);
    Y = Y + a*kron(eye(n3),eye(3*np)-P); 
end
    
end

function Ctp = LOCAL_advance_center(Ct,dt,VW)

Ctp = Ct + dt*VW(1:3,:).';

end

function [colinfo,dt,Ctp] = LOCAL_collision_info(Ct,colinfo,params,dt,VW)

params = params.parsq; 
Ctp = LOCAL_advance_center(Ct,dt,VW);
distC0 = LOCAL_CenterDistance(Ct); 
distC = LOCAL_CenterDistance(Ctp); 
rsum = repmat(params.rd,1,params.n3); rsum = rsum+rsum.'; 
distC = distC - rsum; 
distC0 = distC0 - rsum; 
mindst0 = min(distC0(:)); 
[mindst,mid] = min(distC(:)); 

ceps = min(colinfo.eps,mindst0/2); 
display(mindst0)
display(ceps)

if isempty(colinfo.ip)
    [ip,jp]=meshgrid(1:params.n3);
    display(mid)
    display(mindst)
else
     mid = colinfo.midx; ip = colinfo.ip; jp = colinfo.jp; 
end

Cd = Ct(ip(mid),:) - Ct(jp(mid),:); 
VWt = VW(1:3,:).'; 
Vd = VWt(ip(mid),:) - VWt(jp(mid),:);

c = sum(Cd.^2) - ceps^2; 
b = 2*Cd*Vd.'; 
a = sum(Vd.^2); 

dt = (-b + sqrt(b^2 - 4*a*c))/(2*a);

Ctp = LOCAL_advance_center(Ct,dt,VW); 

end

function colinfo = LOCAL_collision_update(Ct,colinfo,params)

ctb = 4/3;
shflg = isfield(params,'parsh'); 

ns = params.parsq.n3;
distC = LOCAL_CenterDistance(Ct); 
rsum = repmat(params.parsq.rd,1,ns); rsum = rsum+rsum.'; 
distC = distC - rsum; 
mindst = min(distC(:)); 
fprintf('\n Minimum distance (spheres): %2.4f ',mindst);
colevent = mindst<=ctb*colinfo.eps; 

if shflg
    shrd = params.parsh.rd; 
    nrmC = sqrt(Ct(:,1).^2+Ct(:,2).^2+Ct(:,3).^2); 
    distCS = shrd-(params.parsq.rd(:)+nrmC);
    mindstsh = min(distCS); 
    fprintf('\n (Shell) Min relative distance: %2.4f, absolute distance %2.4f ',mindstsh/shrd,mindstsh);
    colevent = colevent || mindstsh<=ctb*colinfo.epsh*shrd; 
end
    
if colevent
    [ii,jj]=meshgrid(1:ns); 
    id = distC <= ctb*colinfo.eps & ii<jj; 
    ip = ii(id); ip=ip(:); 
    jp = jj(id); jp=jp(:); 
    [~,midx] = min(distC(id));  
    
    if shflg
        ii = (1:ns)'; 
        ids = distCS <= ctb*colinfo.epsh*shrd; 
        ips = ii(ids); jps=(ns+1)*ones(size(ips)); 
        ip = [ip(:);ips(:)]; jp = [jp(:);jps(:)];
        colinfo = struct('col',colevent,'ip',ip,'jp',jp,'mindst',mindst,'mindstsh',mindstsh,'eps',colinfo.eps,'epsh',colinfo.epsh,'midx',midx,'type',colinfo.type);
    else
        colinfo = struct('col',colevent,'ip',ip,'jp',jp,'mindst',mindst,'eps',colinfo.eps,'midx',midx,'type',colinfo.type);
    end 
else
    if isfield(params,'parsh')
        colinfo = struct('col',colevent,'ip',[],'jp',[],'mindst',mindst,'mindstsh',mindstsh,'eps',colinfo.eps,'epsh',colinfo.epsh,'type',colinfo.type); 
    else
         colinfo = struct('col',colevent,'ip',[],'jp',[],'mindst',mindst,'eps',colinfo.eps,'type',colinfo.type);
    end
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

function [Xtp,Mtp,Ctp,U,Us,mu,VW,Kernels,params,colinfo,dt,Conc] = LOCAL_euler_step(Xt,X0,Mt,Ct,Kernels,params,colinfo,dt,it,flplot)

%(1) Get velocities (solve BIEs + contact problem)
[U,Us,mu,VW,params,dt,Conc] = LOCAL_compute_velocities(Kernels,Mt,Xt,Ct,params,colinfo,dt,it,flplot);

%(2) Advance centers, rotation, update collision info and operators
[Xtp,Mtp,Ctp,Kernels,params,colinfo,~] = LOCAL_advance_step(VW,mu,Xt,X0,Mt,Ct,Kernels,params,colinfo,dt);

end

function [Xtp,Mtp,Ctp,Kernels,params,colinfo,mindst] = LOCAL_advance_step(VW,mu,Xt,X0,Mt,Ct,Kernels,params,colinfo,dt)
MRot = @(wh,t) RotationMat(wh,t); 
ns = params.parsq.n3; 
np = params.parsq.np;
ctb=(4/3); 

%(2) Advance center
Ctp = LOCAL_advance_center(Ct,dt,VW);
            
%(3) Check for collisions 
% Distance between centers
distC = LOCAL_CenterDistance(Ctp);
rsum = repmat(params.parsq.rd,1,ns); rsum = rsum+rsum.'; 
distC = distC - rsum; 
mindst = min(distC(:)); 

nocolevent = mindst > ctb*colinfo.eps; 

if isfield(params,'parsh')
   shelldst = min(params.parsh.rd-(sqrt(sum(Ctp.^2,2))+params.parsq.rd)); 
   fprintf('\n Minimum distance to wall: %2.4f ',shelldst);
   nocolevent = nocolevent && shelldst > ctb*colinfo.epsh*params.parsh.rd; 
end

if nocolevent
    fprintf('\n Minimum distance (spheres): %2.4f ',mindst);
    if isfield(params,'parsh')
        colinfo = struct('col',false,'ip',[],'jp',[],'mindst',mindst,'eps',colinfo.eps,'epsh',colinfo.epsh,'type',colinfo.type); 
    else
        colinfo = struct('col',false,'ip',[],'jp',[],'mindst',mindst,'eps',colinfo.eps,'type',colinfo.type); 
    end
else    
    if mindst < 0
        [colinfo,dt,Ctp] = LOCAL_collision_info(Ct,colinfo,params,dt,VW); 
    end
    colinfo = LOCAL_collision_update(Ctp,colinfo,params);
end
            
%(4) Advance rotation 
[Mtp,Xtp,~] = LOCAL_advance_rotation(MRot,VW,Mt,Xt,X0,dt,np,ns); 
            
%(5) Update operators
[Kernels,params] = LOCAL_update_operators(Ctp,Mt,params,Kernels);

end

function LOCAL_plotsph(X,phi,colrg,axvec)

warning off;
p = (sqrt(1+2*size(X,1))-1)/2; 
np = 2*p*(p+1); 
mxc=0; 

if ~isempty(phi)
    px = phi(1:3:3*np);
    py = phi(2:3:3*np);
    pz = phi(3:3:3*np);
    color = max(([px py pz]'))'; 
    mxc = max(color(:)); 
    mnc = min(color(:)); 
else
    color=[]; 
end


if mxc>0 && ~isempty(color)
   %color = (1/mxc)*color;  
end

if nargin<4
diam = max(max(X)-min(X)); 
axvec = repmat([-1.5*diam 1.5*diam],1,3); 
end

%figure; 
plotb(X(:),color);
axis(axvec); 
if isempty(colrg)
    caxis([mnc,mxc]);
else
    caxis(colrg); 
end
view([1 1 1])

end

function LOCAL_plotspheres(Xp,p,C,rd,phi)

np = 2*p*(p+1); ns = size(C,1); 

if ns>1
mxc = max(C); mnc = min(C); 
else
mxc = C; mnc = C;     
end

diam=2*max(rd); 
axvec = [mnc(1) mxc(1) mnc(2) mxc(2) mnc(3) mxc(3)] + repmat([-1.5*diam 1.5*diam],1,3); 

px = phi(1:3:end); py = phi(2:3:end); pz = phi(3:3:end);
color = max(([px py pz]'))'; 
mxc = max(color(:)); 
mnc = min(color(:)); 
colrg = [mnc mxc]; 

for k=1:ns
ind = np*(k-1)+1:np*k;
indV = 3*np*(k-1)+1:3*np*k;
LOCAL_plotsph(Xp(ind,:),phi(indV),colrg,axvec); 
hold on;  
end
hold off; 

end

function [C,CA,D]=LOCAL_Build_AuxMats(Wg,Xg,Xc,np,n3)

N = 3*np*n3; 
IC = []; JC = []; VC = []; 
VD = []; VCA = [];  

for i=1:n3
    
idx = (1:np)+np*(i-1); 
X = Xg(idx,:);
W = Wg(idx); 

oW = ones(size(W));
sumW = sum(W); 

%W2 = zeros(1,3*np); 
%W2(1:3:3*np)=W; W2(2:3:3*np)=W; W2(3:3:3*np)=W; 
    
if ~isempty(Xc)
% Center X
X = X - repmat(Xc(i,:),np,1); 
end

indx =(1:3:3*np)+3*np*(i-1); 
indy =(2:3:3*np)+3*np*(i-1); 
indz =(3:3:3*np)+3*np*(i-1);

% C*sigma = [int{sigma} ; int{X \times sigma}]    
IC = [IC; reshape(repmat(6*(i-1)+(1:6),np,1),[],1) ; reshape(repmat(6*(i-1)+(4:6),np,1),[],1)];
JC = [JC; indx'; indy'; indz'; indy'; indz'; indx'; indz'; indx'; indy']; 
VC = [VC; W; W; W ; -W.*X(:,3); -W.*X(:,1); -W.*X(:,2); W.*X(:,2);  W.*X(:,3); W.*X(:,1)];
% First three vectors are just integrals of f for each coordinate
VD = [VD; oW; oW; oW ; -X(:,3); -X(:,1); -X(:,2); X(:,2);  X(:,3); X(:,1)];

tau1 = sum(W.*X(:,3).^2)+sum(W.*X(:,2).^2); 
tau2 = sum(W.*X(:,1).^2)+sum(W.*X(:,3).^2); 
tau3 = sum(W.*X(:,2).^2)+sum(W.*X(:,3).^2); 

VCA = [VCA; W./sumW; W./sumW; W./sumW ; ...
    -W.*X(:,3)./tau1; -W.*X(:,1)./tau2; -W.*X(:,2)./tau3; ...
     W.*X(:,2)./tau1;  W.*X(:,3)./tau2; W.*X(:,1)./tau3];

end

C = sparse(IC,JC,VC,6*n3,N); 
CA = sparse(IC,JC,VCA,6*n3,N);
D = sparse(IC,JC,VD,6*n3,N);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
function Y = LOCAL_apply_joint_system(sigma,mu,Kernels,Ct,params,colinfo)

shf = @(q) real(VshAna([q(1:3:end);q(2:3:end);q(3:3:end)],'VW'));
Ishf = @(qh) reshape(reshape(real(VshSyn(qh,'VW')),[],3).',[],1);

np = params.parsq.np; ns = params.parsq.n3;
nps = params.parsh.np; 
Nsh = 3*nps; Nsq = 3*np*ns; Nvw = 6*ns;
if colinfo.col && ~isempty(colinfo.ip)
   ncol = length(colinfo.ip); N = Nsh+Nsq+Nvw+ncol;  
else
   N = Nsh+Nsq+Nvw;  
end

Y = zeros(N,1); 

% Evaluation at shell
Xtrs = params.parsh.Xp; Nrtrs = params.parsh.Nrp;
Y(1:Nsh) = sigma + Ishf(params.parsh.eigI.*shf(VSh_MatVec_RB_trg(mu(1:Nsq),Xtrs,Nrtrs,params.parsq)));
% Evaluation at squirmer surfaces

% Modify operator if collisions are present
if colinfo.col && ~isempty(colinfo.ip)
    fprintf('\n Computing contact force at collision sites:')
    ip = colinfo.ip; jp = colinfo.jp; 
    
    if params.parsq.dense 
       Em = zeros(ncol,6*ns);  
       
       for k=1:ncol
          % R = x_i - x_j 
          Rhat = Ct(ip(k),:)-Ct(jp(k),:); Rhat=Rhat./norm(Rhat); 
          indk = [(1:3)+6*(ip(k)-1) (1:3)+6*(jp(k)-1)];   
          Em(k,indk) = [Rhat -Rhat]; 
       end
       
       Fm = Em.';
       
       Kernels.M = [Kernels.M [zeros(Nmu,ncol) ; -Fm] ...
           ; zeros(ncol,Nmu) Em zeros(ncol,ncol)]; 
    else 
       Em = zeros(ncol,Nvw);  
        
       for k=1:ncol
          % R = x_i - x_j 
          Rhat = Ct(ip(k),:)-Ct(jp(k),:); Rhat=Rhat./norm(Rhat); 
          indk = [(1:3)+6*(ip(k)-1) (1:3)+6*(jp(k)-1)];   
          Em(k,indk) = [Rhat -Rhat]; 
       end
       
       Fm = Em.'; 
        
       nm = size(Kernels.M,2);  
       Kernels.M = @(V) Kernels.M(V(1:nm)) + ...
           [zeros(nA,1) ; -Fm*V(Nmu+1:Nmu+Nvw) ; Em*V(Nmu+Nvw+1:end)];  
    end
end

%Y(Nsh+1:end) = Lapp(Kernels.M,mu);
Xtrsq = params.parsq.Xp; Nrtrsq = params.parsq.Nrp; params.parsh.Vh = shf(sigma); 
Nmu = size(mu,1); 
Y(Nsh+1:end) = mu + Lslv(Kernels.M,[real(VSh_MatVec_RB_trg(sigma,Xtrsq,Nrtrsq,params.parsh)); zeros(Nmu-Nsq,1)]);
%Y(Nsh+1:Nsh+Nsq) = Y(Nsh+1:Nsh+Nsq) + real(VSh_MatVec_RB_trg(sigma,Xtrsq,Nrtrsq,params.parsh)); 

end

function [U,Us,VW,mu,params] = LOCAL_solve_joint_system(Kernels,Mt,Xt,Ct,params,colinfo,dt,it,flplot)

shf = @(q) real(VshAna([q(1:3:end);q(2:3:end);q(3:3:end)],'VW'));
p = params.parsq.p; 
np = params.parsq.np; ns = params.parsq.n3;
B1 = params.parsq.B1; B2 = params.parsq.B2;
nps = params.parsh.np; 
Nsh = 3*nps; Nsq = 3*np*ns; Nvw = 6*ns;
if colinfo.col && ~isempty(colinfo.ip)
   ncol = length(colinfo.ip); 
   N = Nsh+Nsq+Nvw+ncol;  
else
   N = Nsh+Nsq+Nvw;  
end

rhs = zeros(N,1);  

%compute rhs and solve augmented system 
Xrp = params.parsq.Xrp; 
for k=1:ns
   ind = np*(k-1)+1:np*k; 
   Xrp(ind,:) = Xrp(ind,:)*Mt{k}; 
end
    
% spherical coordinates
[th,phi,~] = cart2sph(Xrp(:,1),Xrp(:,2),Xrp(:,3));  
th(th<0)=th(th<0)+2*pi; phi=pi/2-phi;
    
Us0 = LOCAL_Uslip(B1,B2,phi,th,p,ns); 
Us0 = reshape(Us0.',[],1);
Us = Us0; 

rhs(Nsh+1:Nsh+Nsq) = Us0;
rhs(Nsh+1:end) = Lslv(Kernels.M,rhs(Nsh+1:end)); 

[~,~,mui,VWi] = LOCAL_compute_velocities(Kernels,Mt,[],Ct,params,colinfo,dt,1,0);
pari = LOCAL_compute_shell_velocity(mui,Kernels,params);
x0 = zeros(N,1); 
x0(1:(Nsh+Nsq+Nvw)) = [pari.parsh.shellden ; mui ; VWi(:)]; 

% Solve joint system
JS = @(V) LOCAL_apply_joint_system(V(1:Nsh),V(Nsh+1:end),Kernels,Ct,params,colinfo);
Vr = rand(N,1); 
JS(Vr); 
den = Lslv(JS,rhs,x0); 

% Get densities and other solution components 
sigma = den(1:Nsh); 
params.parsh.shellden = sigma; 
params.parsh.Vh = shf(sigma); 

Y = den(Nsh+1:end); 
mu = Y(1:Nsq); 
VW = reshape(real(Y(Nsq+1:Nsq+Nvw)),6,[]);

Xtrs = params.parsh.Xp; Nrtrs = params.parsh.Nrp;
Usqr = real(VSh_MatVec_RB_trg(mu,Xtrs,Nrtrs,params.parsq));
Uh = shf(Usqr); 
params.parsh.U = Usqr; params.parsh.Uh = Uh; 

if colinfo.col && ~isempty(colinfo.ip)
   %Compute FT and display 
   gamma = Y(Nsq+Nvw+1:end);
   FT = real(reshape(Fm*gamma,6,[])); 
   fprintf('\n Contact Forces \n');
   display(FT)
end

fprintf('\n Rigid body velocities \n');
display(VW)

% surface velocity U 
U = Lapp(Kernels.SDL,mu); 

end
%}
function Test_Spharm_Stk(p,n,m,out)

if nargin==3
    out=2; 
end

% Setup 
shf = @(q) VshAna([q(1:3:end);q(2:3:end);q(3:3:end)],'VW'); 
Sc = SurfaceSph(shape_gallery(p,''));
Xp = reshape(Sc.cart.to_array,[],3);  
np=2*p*(p+1); sp=(p+1)^2;  
ii = (1:sp)'; 
nn = floor(sqrt(ii-1)); 
mm=ii-nn.^2-nn-1;
nnv = repmat(nn,3,1); 
mmv = repmat(mm,3,1); 
nmid = reshape(repmat(1:3,sp,1),[],1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stokes

% Stokes SL, S'
SMat = kernelS([],Sc); 
SpMat = kerneldS_RI(Sc,[]); 
DMat = kernelD([],Sc);
TMat = SpMat; 

[u,v]=gl_grid(p); 

%First example: Vnm,Wnm
V  = Vnm('Vnm',Sc,n,m,u,v); V = reshape(V,3,[]).'; 
W  = Vnm('Wnm',Sc,n,m,u,v); W = reshape(W,3,[]).'; 
X  = Vnm('Xnm',Sc,n,m,u,v); X = reshape(X,3,[]).'; 

%{
MD = -0.5*eye(3*np)+DMat;
[UM,SM,VM] = svd(MD); 
Uc = UM(:,end); Vc=VM(:,end); 
J = [1:3:3*np 2:3:3*np 3:3:3*np]; 
Uch = VshAna(Uc(J,:),'VW'); 
Vch = VshAna(Vc(J,:),'VW'); 
figure; plot(log10(abs(Uch)),'ok'); 
figure; plot(log10(abs(Vch)),'or'); 
pause; 

% Smooth Quadrature Weights (GL x Trapezoidal)
[~, gwt]=g_grid(p+1);
wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
wt = wt(:);
% Area element
W = Sc.geoProp.W; 
W = W.*wt; 
Wv = repmat(W,1,3)'; Wv = Wv(:);
%SW = sum(Wv); 
Vv = reshape(V.',[],1);
MD = MD + Vv*(Vv'.*Wv'); 

display(size(MD))
display(rank(MD))
pause; 
sig = MD\reshape(V.',[],1);
plot(real(sig))
J = [1:3:3*np 2:3:3*np 3:3:3*np]; 
sigh = VshAna(sig(J,:),'VW'); 
figure; plot(log10(abs(sigh)),'o'); 
pause; 
%}

% Check eigenvalues
SVeg  = n/((2*n+1)*(2*n+3)); 
SWeg  = (n+1)/((2*n+1)*(2*n-1));
SpVeg = (3/2)/((2*n+1)*(2*n+3));
SpWeg = -(3/2)/((2*n+1)*(2*n-1));
SXeg  = 1/(2*n+1); 
SpXeg = -(3/2)/(2*n+1); 

fprintf('\n Check eigenvalues for S and Sp on sphere \n')
SV=reshape(SMat*reshape(V.',[],1),3,[]).'; 
SW=reshape(SMat*reshape(W.',[],1),3,[]).';
SX=reshape(SMat*reshape(X.',[],1),3,[]).';
fprintf('\n S[V],S[W] and S[X] \n')
fprintf('l_nm error for S[V]: %e\n',norm(SV-SVeg*V))
fprintf('l_nm error for S[W]: %e\n',norm(SW-SWeg*W))
fprintf('l_nm error for S[X]: %e\n',norm(SX-SXeg*X))
 
SpV=reshape(SpMat*reshape(V.',[],1),3,[]).'; 
SpW=reshape(SpMat*reshape(W.',[],1),3,[]).';  
SpX=reshape(SpMat*reshape(X.',[],1),3,[]).';  
fprintf('\n T[V],T[W] and T[X] \n')
fprintf('l_nm error for T[V]: %e\n',norm(SpV-SpVeg*V))
fprintf('l_nm error for T[W]: %e\n',norm(SpW-SpWeg*W))
fprintf('l_nm error for T[X]: %e\n',norm(SpX-SpXeg*X))

M = [reshape(V.',[],1) reshape(W.',[],1) reshape(X.',[],1)]; 


pause; 
% Far-evaluation check

%Kernel parameters
par = Kernel_Eval_parameters('SL_Stk_3D',0,1,1,1,1e-8,2,400,1);
par.dim = 3; par.mu=1; 
[~, gwt]=g_grid(p+1);
wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
wt = wt(:);
WT = Sc.geoProp.W; WT= WT.*wt; 
Nr = reshape(Sc.geoProp.nor.to_array,[],3);
Xv = reshape(repmat(Xp,1,3)',3,[])'; 
Nrv = reshape(repmat(Nr,1,3)',3,[])';    
Wv = repmat(WT,1,3)'; Wv = Wv(:); 
par.X = Xv; par.nor = Nrv; par.W2 = Wv.'; 
par.ci = repmat((1:3)',np,1); 
par.cj = repmat((1:3)',np,1); 

if out>0
% Target points outside sphere
rt=4;
Xtrg = rt*Xv;
Rtrg = sqrt(sum(Xtrg.'.*Xtrg.')).'; 

%Compare kernel eval with spherical harmonic expansion for V 
KV = Kernel_Eval(Xtrg,Xv,par)*reshape(V.',[],1); 

SVext = @(R,n) SVeg*R.^(-n-2); 
KV2 = reshape(V.',[],1).*SVext(Rtrg,n); 

fprintf('Far eval S[V] outside sphere r=%1.1f ->',rt)
fprintf('error for S[V]: %e\n',norm(KV-KV2))

%Compare kernel eval with spherical harmonic expansion for W 
KW = Kernel_Eval(Xtrg,Xv,par)*reshape(W.',[],1); 

SWVext = @(R,n) (n/(4*n+2))*(R.^(-n-2)-R.^(-n));
SWWext = @(R,n) SWeg*R.^(-n);
KW2 = reshape(V.',[],1).*SWVext(Rtrg,n)+reshape(W.',[],1).*SWWext(Rtrg,n); 

fprintf('Far eval S[W] outside sphere r=%1.1f ->',rt)
fprintf('error for S[W]: %e\n',norm(KW-KW2))

% Compare kernel eval with spherical harmonic expansion for X
KX = Kernel_Eval(Xtrg,Xv,par)*reshape(X.',[],1);

fprintf('Far eval S[X] outside sphere r=%1.1f ->',rt)
SXext = @(R,n) SXeg*R.^(-n-1);
KX2 = reshape(X.',[],1).*SXext(Rtrg,n); 
fprintf('error for S[X]: %e\n',norm(KX-KX2))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Far evaluation for Sp
par.flag_pot='dSL_Stk_3D'; 
dKV = Kernel_Eval(Xtrg,Xv,par)*reshape(V.',[],1); 

SpVext = @(R,n) (-n-2)*SVeg*R.^(-n-3);
dKV2 = reshape(V.',[],1).*SpVext(Rtrg,n); 

fprintf('Far eval dSL[V] outside sphere r=%1.1f ->',rt)
fprintf('error for dSL[V]: %e\n',norm(dKV-dKV2))

dKW = Kernel_Eval(Xtrg,Xv,par)*reshape(W.',[],1);  

SpWVext = @(R,n) (n/(4*n+2))*((-n-2)*R.^(-n-3)-(-n)*R.^(-n-1));
SpWWext = @(R,n) (-n)*SWeg*R.^(-n-1);
dKW2 = reshape(V.',[],1).*SpWVext(Rtrg,n)+reshape(W.',[],1).*SpWWext(Rtrg,n); 

fprintf('Far eval dSL[W] outside sphere r=%1.1f ->',rt)
fprintf('error for dSL[W]: %e\n',norm(dKW-dKW2))

% Check what dSL[W] is like far away from sphere
dKX = Kernel_Eval(Xtrg,Xv,par)*reshape(X.',[],1);

SpXext = @(R,n) -(n+1)*SXeg*R.^(-n-2);
dKX2 = reshape(X.',[],1).*SpXext(Rtrg,n); 

fprintf('Far eval dSL[X] outside sphere r=%1.1f ->',rt)
fprintf('error for dSL[X]: %e\n',norm(dKX-dKX2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Far evaluation for Traction kernel
par.flag_pot='TSL_Stk_3D'; 
TKV = Kernel_Eval(Xtrg,Xv,par)*reshape(V.',[],1); 

TVext = @(R,n) 2*(-n-2)*SVeg*R.^(-n-3);
%cfV = ( ((3*n+2)/(2*n+1))*SpVext(Rtrg,n)+( (n*(-n-2)./((2*n+1)*Rtrg))).*SVext(Rtrg,n)); 
cfV = TVext(Rtrg,n);  
%cfW = ((n+1)/(2*n+1))*(-SpVext(Rtrg,n) + ((-n-2)./Rtrg).*SVext(Rtrg,n)) = 0 

TKV2 = cfV.*reshape(V.',[],1); 

fprintf('Far eval T[V] outside sphere r=%1.1f ->',rt)
fprintf('error for T[V]: %e\n',norm(TKV-TKV2))

TKW = Kernel_Eval(Xtrg,Xv,par)*reshape(W.',[],1); 

SpWVext = @(R,n) (n/(4*n+2))*((-n-2)*R.^(-n-3)-(-n)*R.^(-n-1));
SpWWext = @(R,n) (-n)*SWeg*R.^(-n-1);    
PWext = @(R,n) -n*R.^(-n-1); 

cfV = ( (1+(n+1)/(2*n+1))*SpWVext(Rtrg,n) - (1/(2*n+1))*(n*SpWWext(Rtrg,n)+PWext(Rtrg,n)) + ...
      (n/((2*n+1)*Rtrg))*((-n-2)*SWVext(Rtrg,n) + (n-1)*SWWext(Rtrg,n))); 
cfW = ( (1+(n)/(2*n+1))*SpWWext(Rtrg,n) + (1/(2*n+1))*((-n-1)*SpWVext(Rtrg,n)+PWext(Rtrg,n)) + ...
      ((n+1)/((2*n+1)*Rtrg))*((-n-2)*SWVext(Rtrg,n) + (n-1)*SWWext(Rtrg,n))); 
TKW2 = cfV.*reshape(V.',[],1) + cfW.*reshape(W.',[],1); 

fprintf('Far eval T[W] outside sphere r=%1.1f ->',rt)
fprintf('error for T[W]: %e\n',norm(TKW-TKW2))

% Check what T[W] is like far away from sphere
TKX = Kernel_Eval(Xtrg,Xv,par)*reshape(X.',[],1);

TXext = @(R,n) -(n+2)*SXeg*R.^(-n-2); 
TKX2 = TXext(Rtrg,n).*reshape(X.',[],1); 

fprintf('Far eval T[X] outside sphere r=%1.1f ->',rt)
fprintf('error for T[X]: %e\n',norm(TKX-TKX2))


%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kernel Eval test for random density
%3 random densities of order p-1
pmx=p-1; 
qVh = [rand(pmx^2,1); zeros((p+1)^2-pmx^2,1)];
qWh = [rand(pmx^2,1); zeros((p+1)^2-pmx^2,1)];
qXh = [rand(pmx^2,1); zeros((p+1)^2-pmx^2,1)];
qh = [qVh;qWh;qXh]; 
Q = VshSyn(qh,'VW'); Q = reshape(Q,[],3);   
q = reshape(Q.',[],1);
q2 = [q(1:3:end);q(2:3:end);q(3:3:end)]; 
nrmq=norm(q); 

Ntrg = reshape(Sc.cart.to_array,[],3);  

TqVsh = Vsh_Kernel_Eval_off(qh,'SpMat',0,1,ones(size(u)),u,v,Ntrg);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kernel Eval (self-interaction, need to upsample)
pnew=p+1; 
% Setup 
Scn = SurfaceSph(shape_gallery(pnew,''));
[un,vn] = gl_grid(pnew); 
SpMatn = kerneldS_RI(Scn,[]); 
spn = (pnew+1)^2; 
zr = zeros(spn-sp,1); 
qhn = [qVh;zr;qWh;zr;qXh;zr];
Qn = VshSyn(qhn,'VW'); Qn = reshape(Qn,[],3);  
qn = reshape(Qn.',[],1);
qn2 = [qn(1:3:end);qn(2:3:end);qn(3:3:end)]; 

Ntrgn = reshape(Scn.cart.to_array,[],3);  

a = -0.5; 
tmp=a*qn+SpMatn*qn; 
Tqn = zeros(3*np,1); TqnVsh = Sqn; 
tmpVsh = Vsh_Kernel_Eval_off(qhn,'SpMat',0,1,ones(size(un)),un,vn,Ntrgn);

for k=1:3
Tqn(k:3:end)=interpsh(tmp(k:3:end),p); 
TqnVsh(k:3:end)=interpsh(tmpVsh(k:3:end),p);
end

fprintf('\n T[Q] on sphere (upsampled) from outside r=%1.1f ->',1)
if nrmq>1e-10
    fprintf('\n error for upsampled T[Q] (Vsh_KE_off)): %e',norm(Tqn-TqnVsh)/nrmq)
    fprintf('\n error for T[Q] (Vsh_KE_off)): %e',norm(Tqn-TqVsh)/nrmq)
else
    fprintf('\n error for upsampled T[Q] (Vsh_KE_off)): %e',norm(Tqn-TqnVsh))
    fprintf('\n error for T[Q] (Vsh_KE_off)): %e',norm(Tqn-TqVsh))
end

pause; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Points on another sphere
rnd = 0.1*repmat(rand(np,1),1,3); 
Xtrg2=(rt+rnd).*Xp; Xvtrg2=reshape(repmat(Xtrg2,1,3)',3,[])'; 

[th,phi,rho] = cart2sph(Xtrg2(:,1),Xtrg2(:,2),Xtrg2(:,3)); 
th(th<0)=th(th<0)+2*pi; phi=pi/2-phi; 

Nrtrg = Xtrg2./repmat(rho,1,3); 
par.nor = reshape(repmat(Nr,1,3)',3,[])'; 

% Kernel Eval
TKq = Kernel_Eval(Xvtrg2,Xv,par)*q; 

% Spherical harmonic eval
TKqVsh = Vsh_Kernel_Eval_off(q2,'TSMat',1,1,rho,phi,th,[]);

fprintf('Far eval T[Q] outside sphere r=%1.1f ->',rt)
fprintf('\n error for T[Q] (Vsh_KE_off) normal perturbation): %e',norm(TKq-TKqVsh)/nrmq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q = reshape(Q.',[],1);
q2 = [q(1:3:end);q(2:3:end);q(3:3:end)]; 
nrmq=norm(q); 

% Points on another sphere
cB = ((rt+1)/norm(cBr))*cBr; 
Xtrg2 = Xp+repmat(cB,np,1); 
Xvtrg2 = Xv+repmat(cB,3*np,1);

[th,phi,rho] = cart2sph(Xtrg2(:,1),Xtrg2(:,2),Xtrg2(:,3)); 
th(th<0)=th(th<0)+2*pi; phi=pi/2-phi; 

Nrtrg = Xtrg2./repmat(rho,1,3); 
par.nor = reshape(repmat(Nrtrg,1,3)',3,[])'; 

% Kernel Eval
TKq = Kernel_Eval(Xvtrg2,Xv,par)*q; 

% Spherical harmonic eval

TKqVsh = Vsh_Kernel_Eval_off(q2,'TSMat',1,1,rho,phi,th,Nrtrg);

%TKqVsh2 = reshape(TVext(repmat(rho,1,3),n).',[],1).*Vnm('Vnm',Sc,n,m,phi,th,Nrtrg); 

fprintf('Far eval T[Q] outside sphere cB=[%1.1f,%1.1f,%1.1f] ->',cB(1),cB(2),cB(3))
fprintf('\n error for T[Q] (Vsh_KE_off) another sph: -> %e',norm(TKq-TKqVsh)/nrmq)

GY = Vnm('GY',[],1,-1,phi,th,ones(size(phi)));
GY = reshape(GY,[],3); 
NrG = sqrt(sum(GY*GY',2)); 
GY = GY./repmat(NrG,1,3); 
par.nor = reshape(repmat(GY,1,3)',3,[])'; 

% Kernel Eval
TKq = Kernel_Eval(Xvtrg2,Xv,par)*q; 
% Spherical harmonic eval
TKqVsh = Vsh_Kernel_Eval_off(q2,'TSMat',1,1,rho,phi,th,GY);

fprintf('\n error for T[Q] diff normal: -> %e',norm(TKq-TKqVsh)/nrmq)
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %{
%Experiments with traction kernel step by step
fprintf('\n Experiments with T[V] (ext) \n')
% u = S[V]
%SV = reshape(SV',[],1); SpV = reshape(SpV',[],1); 
%SW = reshape(SW',[],1); SpW = reshape(SpW',[],1); 

[DuVn,DuVtn] = LOCAL_Find_stress(KV2,dKV2,Sc,Nr,Rtrg); 
[DuWn,DuWtn] = LOCAL_Find_stress(KW2,dKW2,Sc,Nr,Rtrg); 

fprintf('\n DuVt*n in terms of Yn and GY ->')
Y = Ynm(n,m,u,v);
GY = Sc.geoProp.Grad(Y);
GY = reshape(GY.to_array,[],3);

YNr = repmat(Y,1,3).*Nr; 
cfp = reshape((-(n+1)*YNr)',[],1); 
cf = ((-n-2)./Rtrg).*reshape(GY',[],1); 
DuVtn2 = cf.*SVext(Rtrg,n) + cfp.*SpVext(Rtrg,n); 

fprintf('error for DuVtn: %e\n',norm(DuVtn2-DuVtn))

M = [reshape(YNr',[],1) reshape(GY',[],1)]; 
c=M\DuVtn;
c2=M\DuVtn2;
%display(norm(c-c2))

cu = M\DuVn; 
ct = M\TKV; 
%display(ct-c2-cu)

fprintf('\n Check (-pn+Du*n+Dut*n) - T[V] error: %e \n',norm(M*(c2+cu)-TKV))
%DmDt = DuVn+DuVtn; 

fprintf('\n DuWt*n in terms of Yn and GY ->')
cfp = reshape((-(n+1)*YNr)',[],1); 
cf = ((-n-2)./Rtrg).*reshape(GY',[],1); 
cgp = reshape((n*YNr)',[],1); 
cg = ((n-1)./Rtrg).*reshape(GY',[],1); 
DuWtn2 = cf.*SWVext(Rtrg,n) + cfp.*SpWVext(Rtrg,n)+...
         cg.*SWWext(Rtrg,n) + cgp.*SpWWext(Rtrg,n); 

fprintf('error for DuWtn: %e\n',norm(DuWtn2-DuWtn))     

c=M\DuWtn;
c2=M\DuWtn2;
%display(norm(c-c2))
cu = M\DuWn; 
ct = M\TKW; 
%display(ct-c2-cu)

%Pressure
pNr2 = reshape(YNr',[],1).*(-n*Rtrg.^(-n-1));
cp2 = M\pNr2; 

fprintf('\n Check (-pn+Du*n+Dut*n) - T[W] error: %e \n',norm(M*(c2+cu+cp2)-TKW))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n Experiments with pressure \n')
par.flag_pot='PSL_Stk_3D'; 
PKV = Kernel_Eval(Xtrg,Xv,par)*reshape(V',[],1);
PVN = PKV.*reshape(Nr',[],1); 
cPV = M\PVN; 
fprintf('\n Pressure coefs for S[V] (zero) [%e,%e]\n',real(cPV(1)),real(cPV(2)))

PKW = Kernel_Eval(Xtrg,Xv,par)*reshape(W',[],1);
PWN = PKW.*reshape(Nr',[],1); 
cPW = M\PWN; 
fprintf('\n Pressure coefs for S[W] [%e,%e]\n',real(cPW(1)),real(cPW(2)))

PKX = Kernel_Eval(Xtrg,Xv,par)*reshape(X',[],1);
PXN = PKX.*reshape(Nr',[],1); 
cPX = M\PXN; 
fprintf('\n Pressure coefs for S[X] (zero) [%e,%e]\n',real(cPX(1)),real(cPX(2)))
display(norm(M*cPW-PWN))
% %}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{ 
% Far evaluation for Double layer kernel

par2=par; 
par2.flag_pot='DL_Stk_3D'; 
DV = Kernel_Eval(Xtrg,Xv,par2)*reshape(V.',[],1); 

DVext = @(R,n) ((2*n.^2+4*n+3)./((2*n+1).*(2*n+3))).*R.^(-n-2); 
cDV = M\DV; 

fprintf('Far eval D[V] outside sphere r=%1.1f ->',rt); 
%fprintf('error for D[V]: %e\n',norm((cDV(1) - DVext(rt,n))*reshape(V.',[],1)))
DV2 = DVext(rt,n)*reshape(V.',[],1); 
fprintf('error for D[V]: %e\n',norm(DV-DV2));

par2=par; 
par2.flag_pot='DL_Stk_3D'; 
DW = Kernel_Eval(Xtrg,Xv,par2)*reshape(W.',[],1); 
cDW = M\DW; 
%display(cDW)
%display(cDW(2) + (2*(n+1).*(n-1)./((2*n+1).*(2*n-1))) * rt.^(-n) )
%display(cDW(1) + (2*n*(n-1)./(4*n+2)) * (rt.^(-n-2)-rt.^(-n)) )

DWVext = @(R,n) (2*n*(n-1)./(4*n+2)) * (R.^(-n-2)-R.^(-n)); 
DWWext = @(R,n) (2*(n+1).*(n-1)./((2*n+1).*(2*n-1))) * R.^(-n); 

DW2 = DWVext(rt,n).*reshape(V.',[],1) + DWWext(rt,n).*reshape(W.',[],1); 
  
fprintf('Far eval D[W] outside sphere r=%1.1f ->',rt); 
fprintf('error for D[W]: %e\n',norm(DW-DW2));

par2=par; 
par2.flag_pot='DL_Stk_3D'; 
DX = Kernel_Eval(Xtrg,Xv,par2)*reshape(X.',[],1); 
cDX = M\DX; 
%display(cDX)
%display(cDX(3) + ((n-1)./(2*n+1)) * rt.^(-n-1) )

DXext = @(R,n) ((n-1)./(2*n+1)) * R.^(-n-1);  
DX2 = DXext(rt,n).*reshape(X.',[],1); 

fprintf('Far eval D[X] outside sphere r=%1.1f ->',rt); 
fprintf('error for D[X]: %e\n',norm(DX-DX2));
pause;

par2.flag_pot='PDL_Stk_3D'; 
PKV = Kernel_Eval(Xtrg,Xv,par2)*reshape(V.',[],1);
PVN = PKV.*reshape(Nr.',[],1); 
cPV = M\PVN; 
fprintf('\n Pressure coefs for S[V] (zero) [%e,%e,%e]\n',real(cPV(1)),real(cPV(2)),real(cPV(3)))

PKW = Kernel_Eval(Xtrg,Xv,par2)*reshape(W.',[],1);
PWN = PKW.*reshape(Nr.',[],1); 
cPW = M\PWN; 
fprintf('\n Pressure coefs for S[W] [%e,%e,%e]\n',real(cPW(1)),real(cPW(2)),real(cPW(3)))
%display(norm(M*cPW-PWN))

Y = Vnm('YNr',Sc,n,m,u,v); Y = reshape(Y,3,[]).'; Y=reshape(Y.',[],1); 
cPY = Y\PWN; 
%display(cPY)
%display(norm(Y*cPY-PWN))

PKX = Kernel_Eval(Xtrg,Xv,par2)*reshape(X.',[],1);
PXN = PKX.*reshape(Nr.',[],1); 
cPX = M\PXN; 
fprintf('\n Pressure coefs for S[X] (zero) [%e,%e,%e]\n',real(cPX(1)),real(cPX(2)),real(cPX(3)))
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kernel Eval test for random density
%3 random densities of order p-1
pmx = max(n-3,min(3,n)); 
display(pmx)
qVh = [rand(pmx^2,1); zeros((p+1)^2-pmx^2,1)];
qWh = [rand(pmx^2,1); zeros((p+1)^2-pmx^2,1)];
qXh = [rand(pmx^2,1); zeros((p+1)^2-pmx^2,1)]; 
qh = [qVh;qWh;qXh]; 
Q = VshSyn(qh,'VW'); 
Q = reshape(Q,[],3);  
q = reshape(Q.',[],1);
q2 = [q(1:3:end);q(2:3:end);q(3:3:end)]; 
nrmq=norm(q); 

% Self-interaction on sphere (limit of Vsh_KE_off from outside)
a=0; 
Sq=a*q+SMat*q; 
Sqsh = VshAna([Sq(1:3:end);Sq(2:3:end);Sq(3:3:end)],'VW'); 
SqVsh = Vsh_Kernel_Eval_off(q2,'SMat',1,1,ones(size(u)),u,v,[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pnew=p+1; 
% Setup 
Scn = SurfaceSph(shape_gallery(pnew,''));
[un,vn] = gl_grid(pnew); 
SMatn = kernelS([],Scn); 
spn = (pnew+1)^2; 
zr = zeros(spn-sp,1); 
qhn = [qVh;zr;qWh;zr;qXh;zr];
Qn = VshSyn(qhn,'VW'); Qn = reshape(Qn,[],3);  
qn = reshape(Qn.',[],1);
qn2 = [qn(1:3:end);qn(2:3:end);qn(3:3:end)]; 

tmp=a*qn+SMatn*qn; 
Sqn = zeros(3*np,1); SqnVsh = Sqn; 
tmpVsh = Vsh_Kernel_Eval_off(qn2,'SMat',1,1,ones(size(un)),un,vn,[]);

for k=1:3
Sqn(k:3:end)=interpsh(tmp(k:3:end),p); 
SqnVsh(k:3:end)=interpsh(tmpVsh(k:3:end),p);
end

fprintf('\n S[Q] on sphere (upsampled) from outside r=%1.1f ->',1)
if nrmq>1e-10
    fprintf('\n error for upsampled S[Q] (Vsh_KE_off)): %e',norm(Sqn-SqnVsh)/nrmq)
    fprintf('\n error for S[Q] (Vsh_KE_off)): %e',norm(Sqn-SqVsh)/nrmq)
else
    fprintf('\n error for upsampled S[Q] (Vsh_KE_off)): %e',norm(Sqn-SqnVsh))
    fprintf('\n error for S[Q] (Vsh_KE_off)): %e',norm(Sq-SqnVsh))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Random density eval 
pause; 
par.flag_pot = 'SL_Stk_3D'; 
% Points on another sphere
%rnd = 0.1*repmat(rand(np,1),1,3); 
Xtrg2=(rt+0).*Xp; 
Xvtrg2=reshape(repmat(Xtrg2,1,3)',3,[])'; 

% SL Kernel Eval
Kq = Kernel_Eval(Xvtrg2,Xv,par)*q; 

% Spherical harmonic eval
[th,phi,rho] = cart2sph(Xtrg2(:,1),Xtrg2(:,2),Xtrg2(:,3)); 
th(th<0)=th(th<0)+2*pi; phi=pi/2-phi; 
KqVsh = Vsh_Kernel_Eval_off(q2,'SMat',1,1,rho,phi,th,[]);

fprintf('\n Far eval S[Q] outside sphere r=%1.1f ->',rt)
fprintf('\n error for S[Q] (Vsh_KE_off) normal perturbation): %e',norm(Kq-KqVsh)/nrmq)

par.flag_pot = 'DL_Stk_3D';
% DL Kernel Eval
Dq = Kernel_Eval(Xvtrg2,Xv,par)*q; 
DqVsh = Vsh_Kernel_Eval_off(q2,'DMat',1,1,rho,phi,th,par.nor);
fprintf('\n error for D[Q] (Vsh_KE_off) normal perturbation): %e',norm(Dq-DqVsh)/nrmq)

SDq = Kq+Dq; 
SDqVsh = Vsh_Kernel_Eval_off(q2,'SDMat',1,1,rho,phi,th,par.nor);
fprintf('\n error for (S+D)[Q] (Vsh_KE_off) normal perturbation): %e',norm(SDq-SDqVsh)/nrmq)

par.flag_pot = 'TSL_Stk_3D';
% TSL Kernel Eval
Tq = Kernel_Eval(Xvtrg2,Xv,par)*q; 
TqVsh = Vsh_Kernel_Eval_off(q2,'TSMat',1,1,rho,phi,th,par.nor);
fprintf('\n error for T[Q] (Vsh_KE_off) normal perturbation): %e',norm(Tq-TqVsh)/nrmq)

par.flag_pot = 'dSL_Stk_3D';
% dSL Kernel Eval
dSq = Kernel_Eval(Xvtrg2,Xv,par)*q; 
dSqVsh = Vsh_Kernel_Eval_off(q2,'SpMat',1,1,rho,phi,th,par.nor);
fprintf('\n error for Sp[Q] (Vsh_KE_off) normal perturbation): %e',norm(dSq-dSqVsh)/nrmq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.flag_pot = 'SL_Stk_3D';

% Points on another sphere
cBr = 2*rand(1,3)-1; cB = ((rt+1)/norm(cBr))*cBr; 
Xtrg2 = Xp+repmat(cB,np,1); 
Xvtrg2 = Xv+repmat(cB,3*np,1); 

% Kernel Eval
Kq = Kernel_Eval(Xvtrg2,Xv,par)*q; 

% Spherical harmonic eval
[th,phi,rho] = cart2sph(Xtrg2(:,1),Xtrg2(:,2),Xtrg2(:,3)); 
th(th<0)=th(th<0)+2*pi; phi=pi/2-phi; 
Nrtrg = Xtrg2./repmat(rho,1,3); 

KqVsh = Vsh_Kernel_Eval_off(q2,'SMat',1,1,rho,phi,th,Nrtrg);
%KqVsh2 = reshape(SVext(repmat(rho,1,3),n).',[],1).*Vnm('Vnm',Sc,n,m,phi,th,Nrtrg); 

fprintf('\n Far eval S[Q] outside sphere cB=[%1.1f,%1.1f,%1.1f] ->',cB(1),cB(2),cB(3))
fprintf('\n error for S[Q] (Vsh_KE_off) another sph: -> %e',norm(Kq-KqVsh)/nrmq)

par.flag_pot = 'DL_Stk_3D';
% DL Kernel Eval
Dq = Kernel_Eval(Xvtrg2,Xv,par)*q; 
DqVsh = Vsh_Kernel_Eval_off(q2,'DMat',1,1,rho,phi,th,par.nor);
fprintf('\n error for D[Q] (Vsh_KE_off) another sph): %e',norm(Dq-DqVsh)/nrmq)

SDq = Kq+Dq; 
SDqVsh = Vsh_Kernel_Eval_off(q2,'SDMat',1,1,rho,phi,th,par.nor);
SDqVsh2 = Vsh_Kernel_Eval_off(qh,'SDMat',0,1,rho,phi,th,par.nor);
norm(SDqVsh - SDqVsh2)

parSD = par; parSD.flag_pot = 'SDL_Stk_3D'; 
parSD.p = p; parSD.n3 = 1; parSD.kerd=3; parSD.dense=1; parSD.C = [0 0 0]; parSD.out=1; parSD.mdist=4; parSD.a=0; 
parSD.Xp = Xp; 
Utrg = VSh_MatVec_RB_trg(q,Xtrg2,Nr,parSD);
%Uslf = Vsh_MatVec_RB2(q,[],parSD);
fprintf('\n Check for trg eval (S+D)[Q]: %e',norm(Utrg-SDqVsh)/nrmq);

SDM = VSh_MatVec_RB2('Mat',[],parSD);
%v = rand(size(SDM,1)); 
fprintf('\n S+D matrix ')
norm(SMat*q+DMat*q-SDM*q)

%Nrtrgv = reshape(repmat(par.nor,1,kerd)',3,[])';
Uslf = VSh_MatVec_RB2(q,[],parSD)+0.5*q; 
Uslf2 = VSh_MatVec_RB_trg(q,Xp,Nr,parSD);
fprintf('\n Check for self eval (S+D)[Q]: %e',norm(Uslf-Uslf2)/nrmq);
%plot(abs(Uslf)); hold on; plot(abs(Uslf2),'r')

parSD.flag_pot='SL_Stk_3D'; 
SM = VSh_MatVec_RB2('Mat',[],parSD);
%v = rand(size(SDM,1)); 
fprintf('\n S matrix ')
norm(SMat*q-SM*q)

Uslf = VSh_MatVec_RB2(q,[],parSD); 
Uslf2 = VSh_MatVec_RB_trg(q,Xp,Nr,parSD);
fprintf('\n Check for self eval S[Q]: %e',norm(Uslf-Uslf2)/nrmq);

parSD.flag_pot = 'DL_Stk_3D';
DM = VSh_MatVec_RB2('Mat',[],parSD);
ss = svd(-0.5*eye(size(DM,1))+DM); %ss=diag(S); 
ss2 = svd(-0.5*eye(size(DM,1))+DMat); 
plot(log10(ss),'-o'); hold on; plot(log10(ss2),'-or')
v = rand(size(DM,1),1); v = v./norm(v);  
fprintf('\n D matrix ')
err = DMat*v-DM*v;
err2 = DMat*q-DM*q;

sherr = shf(err);
sherr2 = shf(err2);
find(abs(sherr)>1e-10)

figure; 
plot(abs(sherr),'o'); 
figure; 
plot(abs(sherr2),'or'); 

Uslf = VSh_MatVec_RB2(q,[],parSD)+0.5*q; 
Uslf2 = VSh_MatVec_RB_trg(q,Xp,Nr,parSD);
fprintf('\n Check for self eval D[Q]: %e',norm(Uslf-Uslf2)/nrmq);


fprintf('\n error for (S+D)[Q] (Vsh_KE_off) another sph): %e',norm(SDq-SDqVsh)/nrmq)

parSD.flag_pot = 'TSL_Stk_3D'; 
TM = VSh_MatVec_RB2('Mat',[],parSD);
%parSD2 = parSD; parSD2.out=0; 
%TM2 = VSh_MatVec_RB2('Mat',[],parSD2);
%v = rand(size(SDM,1)); 
fprintf('\n T matrix ')
TMat = SpMat; 
norm(TMat*q-TM*q)
%{
sherr = shf(TMat*q-TM2*q);
figure; plot(log10(abs(sherr)),'o'); hold on; 
plot(log10(abs(shf(q))),'or'); hold off; 
norm(TMat*q-TM2*q)
%}

Uslf = VSh_MatVec_RB2(q,[],parSD)-0.5*q; 
Uslf2 = VSh_MatVec_RB_trg(q,Xp,Nr,parSD);
Uslf3 = TM*q - 0.5*q; 
fprintf('\n Check for self eval T[Q]: %e',norm(Uslf-Uslf2)/nrmq);
fprintf('\n Check for self eval T[Q]: %e',norm(Uslf-Uslf3)/nrmq);
pause;  

par.flag_pot = 'TSL_Stk_3D';
% TSL Kernel Eval
Tq = Kernel_Eval(Xvtrg2,Xv,par)*q; 
TqVsh = Vsh_Kernel_Eval_off(q2,'TSMat',1,1,rho,phi,th,par.nor);
fprintf('\n error for T[Q] (Vsh_KE_off) another sph): %e',norm(Tq-TqVsh)/nrmq)

par.flag_pot = 'dSL_Stk_3D';
% dSL Kernel Eval 
dSq = Kernel_Eval(Xvtrg2,Xv,par)*q; 
dSqVsh = Vsh_Kernel_Eval_off(q2,'SpMat',1,1,rho,phi,th,par.nor);
fprintf('\n error for Sp[Q] (Vsh_KE_off) another sph): %e',norm(dSq-dSqVsh)/nrmq)
pause;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if out~=1
% Target points inside sphere
rt=0.1; 
Xtrg = rt*Xv;   
Rtrg = sqrt(sum(Xtrg'.*Xtrg'))'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compare kernel eval with spherical harmonic expansion for V
par.flag_pot='SL_Stk_3D';
KV = Kernel_Eval(Xtrg,Xv,par)*reshape(V',[],1); 

SVVint = @(R,n) SVeg*R.^(n+1); 
SVWint = @(R,n) -((n+1)/(4*n+2))*(R.^(n-1)-R.^(n+1));
KV2 = reshape(V',[],1).*SVVint(Rtrg,n)+reshape(W',[],1).*SVWint(Rtrg,n);  

fprintf('\n Far eval S[V] inside sphere r=%1.2f ->',rt)
fprintf('error for S[V]: %e\n',norm(KV-KV2))

%Compare kernel eval with spherical harmonic expansion for W
KW = Kernel_Eval(Xtrg,Xv,par)*reshape(W',[],1); 

SWint = @(R,n) SWeg*R.^(n-1);
KW2 = reshape(W',[],1).*SWint(Rtrg,n); 

fprintf('Far eval S[W] inside sphere r=%1.2f ->',rt)
fprintf('error for S[W]: %e\n',norm(KW-KW2))

% Compare kernel eval with spherical harmonic expansion for X
KX = Kernel_Eval(Xtrg,Xv,par)*reshape(X',[],1);

fprintf('Far eval S[X] inside sphere r=%1.2f ->',rt)

SXint = @(R,n) SXeg*R.^n;
KX2 = reshape(X',[],1).*SXint(Rtrg,n); 
fprintf('error for S[X]: %e\n',norm(KX-KX2))
%}
%{
M = [reshape(V.',[],1) reshape(W.',[],1) reshape(X.',[],1)];
cX = M\KX;
%display(cX)
%display(cX./SXeg)
%display(norm(KX-M*cX))
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
% Far evaluation for Sp

par.flag_pot='dSL_Stk_3D'; par.nor = Nrv; 
dKV = Kernel_Eval(Xtrg,Xv,par)*reshape(V.',[],1); 

SpVVint = @(R,n) (n+1)*SVeg*R.^(n); 
SpVWint = @(R,n) -((n+1)/(4*n+2))*((n-1)*R.^(n-2)-(n+1)*R.^(n));
dKV2 = reshape(V.',[],1).*SpVVint(Rtrg,n)+reshape(W.',[],1).*SpVWint(Rtrg,n);

fprintf('Far eval dSL[V] inside sphere ->')
fprintf('error for dSL[V]: %e\n',norm(dKV-dKV2))

dKW = Kernel_Eval(Xtrg,Xv,par)*reshape(W.',[],1); 

SpWint = @(R,n) (n-1)*SWeg*R.^(n-2);
dKW2 = reshape(W.',[],1).*SpWint(Rtrg,n); 

fprintf('Far eval dSL[W] inside sphere ->')
fprintf('error for dSL[W]: %e\n',norm(dKW-dKW2))

fprintf('Far eval dSL[X] inside sphere ->')

dKX = Kernel_Eval(Xtrg,Xv,par)*reshape(X.',[],1); 
SpXint = @(R,n) n*SXeg*R.^(n-1);
dKX2 = reshape(X.',[],1).*SpXint(Rtrg,n); 
fprintf('error for dSL[X]: %e\n',norm(dKX-dKX2))
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
% Far evaluation for T and D
par.flag_pot='TSL_Stk_3D'; 
TKV = Kernel_Eval(Xtrg,Xv,par)*reshape(V.',[],1); 

par2=par; 
par2.flag_pot='DL_Stk_3D'; 
DV = Kernel_Eval(Xtrg,Xv,par2)*reshape(V.',[],1); 
%display(M\TKV)
cDV = M\DV; 
%display(cDV)
%display(rt)
%display((cDV(1))./( ((2*n.*(n+2))./((2*n+1).*(2*n+3))) * rt.^(n+1) ))
%display((cDV(2))./( (((n+1).*(n+2))./((2*n+1))) * (-rt.^(n-1)+rt.^(n+1)) ))
DVVint = @(R,n) -(2*n.*(n+2))./((2*n+1).*(2*n+3)) * R.^(n+1);
DVWint = @(R,n) -(((n+1).*(n+2))./((2*n+1))) * (-R.^(n-1)+R.^(n+1));
DV2 = DVVint(rt,n).*reshape(V.',[],1) + DVWint(rt,n).*reshape(W.',[],1);

fprintf('Far eval D[V] inside sphere ->')
fprintf('error for D[V]: %e\n',norm(DV-DV2))

DV3 = Vsh_Kernel_Eval_off(reshape(V,[],1),'DMat',1,0,rt*ones(size(u)),u,v,Nrv);
cD3 = M\DV3; 
fprintf('error for D[V]: %e\n',norm(DV-DV3))

SpVVint = @(R,n) (n+1)*SVeg*R.^(n); 
SpVWint = @(R,n) -((n+1)/(4*n+2))*((n-1)*R.^(n-2)-(n+1)*R.^(n));
PVint = @(R,n) -(n+1)*R.^n; 

cfV = ( (1+(n+1)/(2*n+1))*SpVVint(Rtrg,n) - (1/(2*n+1))*(n*SpVWint(Rtrg,n)+PVint(Rtrg,n)) + ...
      (n/((2*n+1)*Rtrg))*((-n-2)*SVVint(Rtrg,n) + (n-1)*SVWint(Rtrg,n))); 
cfW = ( (1+(n)/(2*n+1))*SpVWint(Rtrg,n) + (1/(2*n+1))*((-n-1)*SpVVint(Rtrg,n)+PVint(Rtrg,n)) + ...
      ((n+1)/((2*n+1)*Rtrg))*((-n-2)*SVVint(Rtrg,n) + (n-1)*SVWint(Rtrg,n))); 
TKV2 = cfV.*reshape(V.',[],1) + cfW.*reshape(W.',[],1); 

%M = [reshape(V',[],1) reshape(W',[],1)];

fprintf('Far eval T[V] inside sphere ->')
fprintf('error for T[V]: %e\n',norm(TKV-TKV2))
  
TKV3 = Vsh_Kernel_Eval_off(reshape(V,[],1),'TSMat',1,0,rt*ones(size(u)),u,v,Nrv);
fprintf('error for T[V]: %e\n',norm(TKV-TKV3))

cTV = M\TKV; 
%display(cTV); 
%display(M\TKV2)
cTV3 = M\TKV3; 
%display(cTV3);
%display((cTV(1)-cTV3(1))/(rt.^n)); 
%display((cTV(2)-cTV3(2))/(rt.^n)); 

pause; 

TKW = Kernel_Eval(Xtrg,Xv,par)*reshape(W.',[],1); 

par2=par; 
par2.flag_pot='DL_Stk_3D'; 
DW = Kernel_Eval(Xtrg,Xv,par2)*reshape(W.',[],1); 
%display(M\TKW)
cDW = M\DW; 
%display(cDW)
%display(rt)
%display((cDW(2))./( ((2*n.^2+1)./((2*n+1).*(2*n-1))) * rt.^(n-1) ))
DWint = @(R,n) -((2*n.^2+1)./((2*n+1).*(2*n-1))) * R.^(n-1); 
DW2 = DWint(rt,n).*reshape(W.',[],1);

SpWint = @(R,n) (n-1)*SWeg*R.^(n-2);
TKW2 = 2*reshape(W.',[],1).*SpWint(Rtrg,n); 

fprintf('Far eval D[W] inside sphere ->')
fprintf('error for D[W]: %e\n',norm(DW-DW2))

DW3 = Vsh_Kernel_Eval_off(reshape(W,[],1),'DMat',1,0,rt*ones(size(u)),u,v,Nrv);
fprintf('error for D[W]: %e\n',norm(DW-DW3))
%norm(DW-DW3)

fprintf('Far eval T[W] inside sphere ->')
fprintf('error for T[W]: %e\n',norm(TKW-TKW2))
TKW3 = Vsh_Kernel_Eval_off(reshape(W,[],1),'TSMat',1,0,rt*ones(size(u)),u,v,Nrv);
fprintf('error for T[W]: %e\n',norm(TKW-TKW3))
%norm(TKW-TKW3)

TKX = Kernel_Eval(Xtrg,Xv,par)*reshape(X.',[],1); 

par2=par; 
par2.flag_pot='DL_Stk_3D'; 
DX = Kernel_Eval(Xtrg,Xv,par2)*reshape(X.',[],1); 
%display(M\TKX)
cDX = M\DX; 
%display(cDX)
%display(rt)
%display((cDX(3))./( ((n+2)./(2*n+1)) * rt.^(n) ))
DXint = @(R,n) -((n+2)./(2*n+1)) * R.^(n);
DX2 = DXint(rt,n)*reshape(X.',[],1); 

TXint = @(R,n) (n-1)*SXeg*R.^(n-1);
TKX2 = TXint(Rtrg,n).*reshape(X.',[],1); 

fprintf('Far eval D[X] inside sphere ->')
fprintf('error for D[X]: %e\n',norm(DX-DX2))

DX3 = Vsh_Kernel_Eval_off(reshape(X,[],1),'DMat',1,0,rt*ones(size(u)),u,v,Nrv);
fprintf('error for D[X]: %e\n',norm(DX-DX3))
%norm(DX-DX3)

fprintf('Far eval T[X] inside sphere ->')
fprintf('error for T[X]: %e\n',norm(TKX-TKX2))

TKX3 = Vsh_Kernel_Eval_off(reshape(X,[],1),'TSMat',1,0,rt*ones(size(u)),u,v,Nrv);
fprintf('error for T[X]: %e\n',norm(TKX-TKX3))
%norm(TKX-TKX3)

par2.flag_pot='PDL_Stk_3D'; 
PKV = Kernel_Eval(Xtrg,Xv,par2)*reshape(V.',[],1);
PVN = PKV.*reshape(Nr.',[],1); 
cPV = M\PVN; 
fprintf('\n Pressure coefs for S[V] [%e,%e,%e]\n',real(cPV(1)),real(cPV(2)),real(cPV(3)))

Y = Vnm('YNr',Sc,n,m,u,v); Y = reshape(Y,3,[]).'; Y=reshape(Y.',[],1); 
cPY = Y\PVN; 
%display(cPY)
%display(norm(Y*cPY-PVN))

PKW = Kernel_Eval(Xtrg,Xv,par2)*reshape(W.',[],1);
PWN = PKW.*reshape(Nr.',[],1); 
cPW = M\PWN; 
fprintf('\n Pressure coefs for S[W] (zero) [%e,%e,%e]\n',real(cPW(1)),real(cPW(2)),real(cPW(3)))
%display(norm(M*cPW-PWN))

PKX = Kernel_Eval(Xtrg,Xv,par2)*reshape(X.',[],1);
PXN = PKX.*reshape(Nr.',[],1); 
cPX = M\PXN; 
fprintf('\n Pressure coefs for S[X] (zero) [%e,%e,%e]\n',real(cPX(1)),real(cPX(2)),real(cPX(3)))
pause; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kernel Eval test for random density
%3 random densities of order p-1
pmx = p-3; 
qVh = [rand(pmx^2,1); zeros((p+1)^2-pmx^2,1)]; 
qWh = [rand(pmx^2,1); zeros((p+1)^2-pmx^2,1)];
qXh = [rand(pmx^2,1); zeros((p+1)^2-pmx^2,1)];

%qVh = zeros(size(qVh)); 
%qWh = zeros(size(qWh));
%qWh(3)=1;
%qXh = zeros(size(qXh));
%qXh(3)=1;

qh = [qVh;qWh;qXh]; 
Q = VshSyn(qh,'VW');
Q = reshape(Q,[],3);
q = reshape(Q.',[],1);
q2 = [q(1:3:end);q(2:3:end);q(3:3:end)]; 
nrmq=norm(q); 

% Points on another sphere
rnd = 0; %0.1*repmat(rand(np,1),1,3); 
Xtrg2=(rt+rnd).*Xp; 
Xvtrg2=reshape(repmat(Xtrg2,1,3)',3,[])'; 

% Kernel Eval
par.flag_pot = 'SL_Stk_3D'; 
Kq = Kernel_Eval(Xvtrg2,Xv,par)*q; 

% Spherical harmonic eval
[th,phi,rho] = cart2sph(Xtrg2(:,1),Xtrg2(:,2),Xtrg2(:,3)); 
th(th<0)=th(th<0)+2*pi; phi=pi/2-phi; 
KqVsh = Vsh_Kernel_Eval_off(q2,'SMat',1,0,rho,phi,th,[]);

fprintf('Far eval S[Q] inside sphere r=%1.1f ->',rt)
fprintf('\n error for S[Q] (Vsh_KE_off) normal perturbation): %e',norm(Kq-KqVsh)/nrmq)

par.flag_pot = 'DL_Stk_3D';
% DL Kernel Eval
Dq = Kernel_Eval(Xvtrg2,Xv,par)*q; 
DqVsh = Vsh_Kernel_Eval_off(q2,'DMat',1,0,rho,phi,th,par.nor);
fprintf('\n error for D[Q] (Vsh_KE_off) normal perturbation): %e',norm(Dq-DqVsh)/nrmq)

SDq = Kq+Dq; 
SDqVsh = Vsh_Kernel_Eval_off(q2,'SDMat',1,0,rho,phi,th,par.nor);
fprintf('\n error for (S+D)[Q] (Vsh_KE_off) normal perturbation): %e',norm(SDq-SDqVsh)/nrmq)

par.flag_pot = 'TSL_Stk_3D';
% TSL Kernel Eval
Tq = Kernel_Eval(Xvtrg2,Xv,par)*q; 
TqVsh = Vsh_Kernel_Eval_off(q2,'TSMat',1,0,rho,phi,th,par.nor);
fprintf('\n error for T[Q] (Vsh_KE_off) normal perturbation): %e',norm(Tq-TqVsh)/nrmq)

par.flag_pot = 'dSL_Stk_3D';
% dSL Kernel Eval
dSq = Kernel_Eval(Xvtrg2,Xv,par)*q; 
dSqVsh = Vsh_Kernel_Eval_off(q2,'SpMat',1,0,rho,phi,th,par.nor);
fprintf('\n error for Sp[Q] (Vsh_KE_off) normal perturbation): %e',norm(dSq-dSqVsh)/nrmq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%q = reshape(Q.',[],1);
%q2 = [q(1:3:end);q(2:3:end);q(3:3:end)]; 
%nrmq=norm(q); 

% Points on another sphere
cBr = 2*rand(1,3)-1; cB = (0.1/norm(cBr))*cBr; 
Xtrg2 = rt*Xp+repmat(cB,np,1); 
Xvtrg2 = rt*Xv+repmat(cB,3*np,1); 

% Kernel Eval
par.flag_pot = 'SL_Stk_3D'; 
Kq = Kernel_Eval(Xvtrg2,Xv,par)*q; 

% Spherical harmonic eval
[th,phi,rho] = cart2sph(Xtrg2(:,1),Xtrg2(:,2),Xtrg2(:,3)); 
th(th<0)=th(th<0)+2*pi; phi=pi/2-phi; 
Nrtrg = Xtrg2./repmat(rho,1,3); 

KqVsh = Vsh_Kernel_Eval_off(q2,'SMat',1,0,rho,phi,th,Nrtrg);

%KqVsh2 = reshape(SVext(repmat(rho,1,3),n).',[],1).*Vnm('Vnm',Sc,n,m,phi,th,Nrtrg); 

fprintf('\n Far eval S[Q] inside sphere r=%1.1f ,cB=[%1.2f,%1.2f,%1.2f] ->',rt,cB(1),cB(2),cB(3))
fprintf('\n error for S[Q] (Vsh_KE_off) another sph: -> %e',norm(Kq-KqVsh)/nrmq)

par.flag_pot = 'DL_Stk_3D';
% DL Kernel Eval
Dq = Kernel_Eval(Xvtrg2,Xv,par)*q; 
DqVsh = Vsh_Kernel_Eval_off(q2,'DMat',1,0,rho,phi,th,par.nor);
fprintf('\n error for D[Q] (Vsh_KE_off) another sph): %e',norm(Dq-DqVsh)/nrmq)

SDq = Kq+Dq; 
SDqVsh = Vsh_Kernel_Eval_off(q2,'SDMat',1,0,rho,phi,th,par.nor);
fprintf('\n error for (S+D)[Q] (Vsh_KE_off) another sph): %e',norm(SDq-SDqVsh)/nrmq)

par.flag_pot = 'TSL_Stk_3D';
% TSL Kernel Eval
Tq = Kernel_Eval(Xvtrg2,Xv,par)*q; 
TqVsh = Vsh_Kernel_Eval_off(q2,'TSMat',1,0,rho,phi,th,par.nor);
fprintf('\n error for T[Q] (Vsh_KE_off) another sph): %e',norm(Tq-TqVsh)/nrmq)

par.flag_pot = 'dSL_Stk_3D';
% dSL Kernel Eval 
dSq = Kernel_Eval(Xvtrg2,Xv,par)*q; 
dSqVsh = Vsh_Kernel_Eval_off(q2,'SpMat',1,0,rho,phi,th,par.nor);
fprintf('\n error for Sp[Q] (Vsh_KE_off) another sph): %e \n',norm(dSq-dSqVsh)/nrmq)
pause; 
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
pause; 
parSD = par; parSD.flag_pot = 'SDL_Stk_3D'; 
parSD.p = p; parSD.n3 = 1; parSD.kerd=3; parSD.dense=1; parSD.C = [0 0 0]; parSD.out=0; parSD.mdist=0.75; parSD.a=0; 
parSD.Xp = Xp; 
Utrg = VSh_MatVec_RB_trg(q,Xtrg2,Nr,parSD);
%Uslf = Vsh_MatVec_RB2(q,[],parSD);
fprintf('\n Check for trg eval (S+D)[Q]: %e',norm(Utrg-SDqVsh)/nrmq);

SDM = VSh_MatVec_RB2('Mat',[],parSD);
%v = rand(size(SDM,1)); 
fprintf('\n S+D matrix ')
norm(SMat*q+DMat*q-SDM*q)

%Nrtrgv = reshape(repmat(par.nor,1,kerd)',3,[])';
Uslf = VSh_MatVec_RB2(q,[],parSD)-0.5*q; 
Uslf2 = VSh_MatVec_RB_trg(q,Xp,Nr,parSD);
fprintf('\n Check for self eval (S+D)[Q]: %e',norm(Uslf-Uslf2)/nrmq);
%plot(abs(Uslf)); hold on; plot(abs(Uslf2),'r')

parSD.flag_pot='SL_Stk_3D'; 
SM = VSh_MatVec_RB2('Mat',[],parSD);
%v = rand(size(SDM,1)); 
fprintf('\n S matrix ')
norm(SMat*q-SM*q)
sherr = shf(SMat*q-SM*q);
figure; 
%plot(abs(sherr),'o');
ier = find(abs(sherr)>1e-10); 
display([nnv(ier) mmv(ier) nmid(ier)])
pause; 

Uslf = VSh_MatVec_RB2(q,[],parSD); 
Uslf2 = VSh_MatVec_RB_trg(q,Xp,Nr,parSD);
fprintf('\n Check for self eval S[Q]: %e',norm(Uslf-Uslf2)/nrmq);

parSD.flag_pot = 'DL_Stk_3D';
DM = VSh_MatVec_RB2('Mat',[],parSD);
%ss = svd(0.5*eye(size(DM,1))+DM); %ss=diag(S); 
%ss2 = svd(0.5*eye(size(DM,1))+DMat); 
%plot(log10(ss),'-o'); hold on; plot(log10(ss2),'-or')
v = rand(size(DM,1),1); v = v./norm(v);  
fprintf('\n D matrix ')
err = DMat*v-DM*v;
err2 = DMat*q-DM*q;
norm(err2)

sherr = shf(err);
ier=find(abs(sherr)>1e-10);
display([nnv(ier) mmv(ier) nmid(ier)])
sherr2 = shf(err2);
ier2=find(abs(sherr2)>1e-10);
display([nnv(ier2) mmv(ier2) nmid(ier2)])
%{
figure; 
plot(abs(sherr),'o'); 
figure; 
plot(abs(sherr2),'or'); 
%}
Uslf = VSh_MatVec_RB2(q,[],parSD)-0.5*q; 
Uslf2 = VSh_MatVec_RB_trg(q,Xp,Nr,parSD);
fprintf('\n Check for self eval D[Q]: %e',norm(Uslf-Uslf2)/nrmq);

fprintf('\n error for (S+D)[Q] (Vsh_KE_off) another sph): %e',norm(SDq-SDqVsh)/nrmq)

parSD.flag_pot = 'TSL_Stk_3D'; 
TM = VSh_MatVec_RB2('Mat',[],parSD);
TMq = VSh_MatVec_RB2(q,[],parSD); 

fprintf('\n T[W] ')
w = reshape(W.',[],1);
TKWtrue = TMat*w; 
TKWq = (2*SpWint(1,n) - 0.5).*w;
TKWq2 = TM*w; 
TKWq3 = VSh_MatVec_RB2(w,[],parSD); 
norm(TKWq-TKWtrue)
norm(TKWq2-TKWtrue)
norm(TKWq3-TKWtrue)
norm(TKWq2-TKWq3)

M = [reshape(V.',[],1) reshape(W.',[],1) reshape(X.',[],1)];
ct = (2*n+1)*(2*n-1); 
display(ct*(M\TKWtrue))
display(ct*(M\TKWq))
display(ct*(M\TKWq2))
display(ct*(M\TKWq3))
pause; 

%v = rand(size(SDM,1)); 
fprintf('\n T matrix ')
TMat = SpMat; 
norm(TMat*q-TM*q)
norm(TMat*q-TMq)
norm(TM*q-TMq)
sherr = shf(TMat*q-TM*q);
%ier=find(abs(sherr)>1e-10);
%display([nnv(ier) mmv(ier) nmid(ier)])
figure; 
plot(log10(abs(shf(q))),'ob'); 
figure; 
plot(log10(abs(sherr)),'ok')
figure; 
plot(log10(abs(shf(TMat*q))),'diamondg'); hold on; 
plot(log10(abs(shf(TM*q))),'oc')
plot(log10(abs(shf(TMq))),'squarer'); 

Uslf = VSh_MatVec_RB2(q,[],parSD)+0.5*q; 
Uslf2 = VSh_MatVec_RB_trg(q,Xp,Nr,parSD);
Uslf3 = TM*q + 0.5*q; 
fprintf('\n Check for self eval T[Q]: %e',norm(Uslf-Uslf2)/nrmq);
fprintf('\n Check for self eval T[Q]: %e',norm(Uslf-Uslf3)/nrmq);
pause; 

par.flag_pot = 'TSL_Stk_3D';
% TSL Kernel Eval
Tq = Kernel_Eval(Xvtrg2,Xv,par)*q; 
TqVsh = Vsh_Kernel_Eval_off(q2,'TSMat',1,0,rho,phi,th,par.nor);
fprintf('\n error for T[Q] (Vsh_KE_off) another sph): %e',norm(Tq-TqVsh)/nrmq)
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
% Kernel Eval test for random density
%3 random densities of order p-1
qVh = [rand(p^2,1); zeros((p+1)^2-p^2,1)];
qWh = [rand(p^2,1); zeros((p+1)^2-p^2,1)];
qXh = [rand(p^2,1); zeros((p+1)^2-p^2,1)];
qh = [qVh;qWh;qXh]; 
Q = VshSyn(qh,'VW');   
q = reshape(Q.',[],1);
q2 = [q(1:3:end);q(2:3:end);q(3:3:end)]; 
nrmq=norm(q); 

% Points on another sphere
rnd = 0.1*repmat(rand(np,1),1,3); 
Xtrg2=(rt+rnd).*Xp; Xvtrg2=reshape(repmat(Xtrg2,1,3)',3,[])'; 

[th,phi,rho] = cart2sph(Xtrg2(:,1),Xtrg2(:,2),Xtrg2(:,3)); 
th(th<0)=th(th<0)+2*pi; phi=pi/2-phi; 

Nrtrg = Xtrg2./repmat(rho,1,3); 
par.nor = reshape(repmat(Nr,1,3)',3,[])'; 

% Kernel Eval
TKq = Kernel_Eval(Xvtrg2,Xv,par)*q; 

% Spherical harmonic eval
TKqVsh = Vsh_Kernel_Eval_off(q2,'TSMat',1,0,rho,phi,th,par.nor);

fprintf('Far eval T[Q] inside sphere r=%1.1f ->',rt)
fprintf('\n error for T[Q] (Vsh_KE_off) normal perturbation): %e',norm(TKq-TKqVsh)/nrmq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q = reshape(Q.',[],1);
q2 = [q(1:3:end);q(2:3:end);q(3:3:end)]; 
nrmq=norm(q); 

% Points on another sphere
cB = (0.1/norm(cBr))*cBr; 
Xtrg2 = rt*Xp+repmat(cB,np,1); 
Xvtrg2 = rt*Xv+repmat(cB,3*np,1); 

[th,phi,rho] = cart2sph(Xtrg2(:,1),Xtrg2(:,2),Xtrg2(:,3)); 
th(th<0)=th(th<0)+2*pi; phi=pi/2-phi; 

Nrtrg = Xtrg2./repmat(rho,1,3); 
par.nor = reshape(repmat(Nrtrg,1,3)',3,[])'; 

% Kernel Eval
TKq = Kernel_Eval(Xvtrg2,Xv,par)*q; 

% Spherical harmonic eval

TKqVsh = Vsh_Kernel_Eval_off(q2,'TSMat',1,0,rho,phi,th,par.nor);

%TKqVsh2 = reshape(TVext(repmat(rho,1,3),n).',[],1).*Vnm('Vnm',Sc,n,m,phi,th,Nrtrg); 

fprintf('Far eval T[Q] inside sphere rt=%1.1f, cB=[%1.2f,%1.2f,%1.2f] ->',rt,cB(1),cB(2),cB(3))
fprintf('\n error for T[Q] (Vsh_KE_off) another sph: -> %e',norm(TKq-TKqVsh)/nrmq)
%fprintf('\n rel error for S another sph (Sh_Kernel_Eval_off2): -> %e',norm(TKq-TKqVsh2)/nrmq)
%fprintf('\n rel error for S another sph (difference): -> %e',norm(TKqVsh-TKqVsh2)/nrmq)
%}

%{
M = [reshape(V',[],1) reshape(W',[],1) reshape(X',[],1)];
cX = M\TKX;
display(cX)
display(cX./SXeg)
display(norm(TKX-M*cX))
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%Experiments with traction kernel step by step
fprintf('\n Experiments with T[V] (ext) \n')
[DuVn,DuVtn] = LOCAL_Find_stress(KV2,dKV2,Sc,Nr,Rtrg); 
[DuWn,DuWtn] = LOCAL_Find_stress(KW2,dKW2,Sc,Nr,Rtrg); 

fprintf('\n DuVt*n in terms of Yn and GY ->')
Y = Ynm(n,m,u,v);
GY = Sc.geoProp.Grad(Y);
GY = reshape(GY.to_array,[],3);

YNr = repmat(Y,1,3).*Nr; 
cfp = reshape((-(n+1)*YNr)',[],1); 
cf = ((-n-2)./Rtrg).*reshape(GY',[],1); 
cgp = reshape((n*YNr)',[],1); 
cg = ((n-1)./Rtrg).*reshape(GY',[],1); 
DuVtn2 = cf.*SVVint(Rtrg,n) + cfp.*SpVVint(Rtrg,n)+...
         cg.*SVWint(Rtrg,n) + cgp.*SpVWint(Rtrg,n);  

fprintf('error for DuVtn: %e\n',norm(DuVtn2-DuVtn))

M = [reshape(YNr',[],1) reshape(GY',[],1)]; 
c=M\DuVtn;
c2=M\DuVtn2;
%display(norm(c-c2))

cu = M\DuVn; 
ct = M\TKV; 

%display(ct-c2-cu)

%fprintf('\n Pressure \n')
pNr2 = reshape(YNr',[],1).*(-(n+1)*Rtrg.^n);
cp2 = M\pNr2; 
%display(cp2)

fprintf('\n Check (-pn+Du*n+Dut*n) - T[V] error: %e \n',norm(M*(c2+cu+cp2)-TKV));   

fprintf('\n DuWt*n in terms of Yn and GY ->')
cgp = reshape((n*YNr)',[],1); 
cg = ((n-1)./Rtrg).*reshape(GY',[],1); 
DuWtn2 = cg.*SWint(Rtrg,n) + cgp.*SpWint(Rtrg,n); 

fprintf('error for DuWtn: %e\n',norm(DuWtn2-DuWtn))

c=M\DuWtn;
c2=M\DuWtn2;
%display(norm(c-c2))

cu = M\DuWn; 
ct = M\TKW; 

%display(ct-c2-cu)

fprintf('\n Check (-pn+Du*n+Dut*n) - T[W] error: %e \n',norm(M*(c2+cu)-TKW)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n Experiments with pressure \n')
par.flag_pot='PSL_Stk_3D'; 
PKV = Kernel_Eval(Xtrg,Xv,par)*reshape(V',[],1);
PVN = PKV.*reshape(Nr',[],1); 
cPV = M\PVN; 
%display(norm(M*cPV-PVN))
fprintf('\n Pressure coefs for S[V] [%e,%e]\n',real(cPV(1)),real(cPV(2)))
%display(real(cPV'))
%display(norm(cPV+cp2))

PKW = Kernel_Eval(Xtrg,Xv,par)*reshape(W',[],1);
PWN = PKW.*reshape(Nr',[],1); 
cPW = M\PWN; 
%display(norm(M*cPW-PWN))
fprintf('\n Pressure coefs for S[W] (zero) [%e,%e]\n',real(cPW(1)),real(cPW(2)))

PKX = Kernel_Eval(Xtrg,Xv,par)*reshape(X',[],1);
PXN = PKX.*reshape(Nr',[],1); 
cPX = M\PXN; 
fprintf('\n Pressure coefs for S[X] (zero) [%e,%e]\n',real(cPX(1)),real(cPX(2)))
%}
%}

end

end

function [Dun,Dutn] = LOCAL_Find_stress(u,dK,Sc,Nr,Rtrg)

u1 = u(1:3:end); u2=u(2:3:end); u3 = u(3:3:end); 
Rtrg = reshape(Rtrg,3,[])'; 

Dgu1 = Sc.geoProp.Grad(u1); Dgu1 = (1./(Rtrg+(Rtrg==0))).*reshape(Dgu1.to_array,[],3);
Dgu2 = Sc.geoProp.Grad(u2); Dgu2 = (1./(Rtrg+(Rtrg==0))).*reshape(Dgu2.to_array,[],3);
Dgu3 = Sc.geoProp.Grad(u3); Dgu3 = (1./(Rtrg+(Rtrg==0))).*reshape(Dgu3.to_array,[],3);

Gu1 = repmat(dK(1:3:end),1,3).*Nr + Dgu1; 
Gu2 = repmat(dK(2:3:end),1,3).*Nr + Dgu2; 
Gu3 = repmat(dK(3:3:end),1,3).*Nr + Dgu3; 

Dun = zeros(size(u)); 
Dun(1:3:end) = sum(Gu1'.*Nr')'; Dun(2:3:end) = sum(Gu2'.*Nr')'; Dun(3:3:end) = sum(Gu3'.*Nr')'; 

Gutx = [Gu1(:,1) Gu2(:,1) Gu3(:,1)]; 
Guty = [Gu1(:,2) Gu2(:,2) Gu3(:,2)];
Gutz = [Gu1(:,3) Gu2(:,3) Gu3(:,3)]; 

Dutn = zeros(size(u)); 
Dutn(1:3:end) = sum(Gutx'.*Nr')'; Dutn(2:3:end) = sum(Guty'.*Nr')'; Dutn(3:3:end) = sum(Gutz'.*Nr')'; 

end
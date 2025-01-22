function [dKq,dKq4]=Test_Spharm_Lap(p,n,m)

% Setup 
Sc = SurfaceSph(shape_gallery(p,''));
X = reshape(Sc.cart.to_array,[],3); 
np = 2*p*(p+1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Laplace

% Laplace SL, S' and DL
[SMat, SpMat, ~] = kernelDLap(Sc); 

[u,v]=gl_grid(p); 

%First example: Ynm
%n=4;m=3;
Y = Ynm(n,m,u,v); 

% Check eigenvalues
Seg  = @(n) (1./(2*n+1)); 
Speg = @(n) -(1./(4*n+2));  
fprintf('\n Check eigenvalues for S and Sp \n')
fprintf('l_nm error for S: %e\n',norm(SMat*Y-Seg(n)*Y))
fprintf('l_nm error for Sp: %e\n',norm(SpMat*Y-Speg(n)*Y))

fprintf('\n Check evaluation for S and Sp on sphere \n')
%random density of order p-1
qh = [rand(p^2,1); zeros((p+1)^2-p^2,1)];
q = shSyn(qh); 
nn=floor(sqrt((0:((p+1)^2-1))))';
Sqh = Seg(nn).*qh;
Spqh = Speg(nn).*qh; 
Tq = shSyn([Sqh Spqh]); 
Sq = Tq(:,1); 
Spq = Tq(:,2);

nrmq = norm(q); 
fprintf('rel error for S: %e\n',norm(Sq-SMat*q)/nrmq)
fprintf('rel error for Sp: %e\n',norm(Spq-SpMat*q)/nrmq)

% Far-evaluation check

%Kernel parameters
par = Kernel_Eval_parameters('SL_L_3D',0,1,1,1,1e-8,2,400,1);
par.dim = 3;
[~, gwt]=g_grid(p+1);
wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
wt = wt(:);
W = Sc.geoProp.W; W= W.*wt; 
Nr = reshape(Sc.geoProp.nor.to_array,[],3); 
par.X = X; par.nor = Nr; par.W2 = W.'; 

% Target points outside sphere
rt=4; 
Xtrg = rt*X;
Rtrg = sqrt(sum(Xtrg'.*Xtrg'))'; 

%Compare kernel eval with spherical harmonic expansion
Y = real(Y); 
KY = Kernel_Eval(Xtrg,X,par)*Y; 
K2 = Seg(n)*Y.*Rtrg.^(-n-1); 

fprintf('\n Check far eval for S[Y] outside sphere r=%1.1f \n',rt)
fprintf('error for S: %e\n',norm(KY-K2))

fprintf('\n Check far eval for S[q] outside sphere r=%1.1f \n',rt)
Kq = Kernel_Eval(Xtrg,X,par)*q; 
Kq2 = shSyn((rt.^(-nn-1)).*Sqh); %(all pts have same radius)

fprintf('rel error for S: %e\n',norm(Kq-Kq2)/nrmq)

fun = @(r,n) r.^(-n-1); 
Xtrg2 = Xtrg+repmat(rand(np,1)-0.5,1,3).*X; 
Rtrg2 = sqrt(sum(Xtrg2'.*Xtrg2'))';
Kq = Kernel_Eval(Xtrg2,X,par)*q; 
Kq3 = LOCAL_off_shSyn(p,fun,Sqh,Rtrg2,u,v); 

fprintf('rel error for S (local shSyn / normal perturbation): %e\n',norm(Kq-Kq3)/nrmq)

Xtrg3 = X+repmat([2 2.5 2],np,1); 
%Rtrg3 = sqrt(sum(Xtrg3'.*Xtrg3'))';
Kq = Kernel_Eval(Xtrg3,X,par)*q;
[th,phi,rho] = cart2sph(Xtrg3(:,1),Xtrg3(:,2),Xtrg3(:,3)); 
th(th<0)=th(th<0)+2*pi; phi=pi/2-phi; 
Kq4 = LOCAL_off_shSyn(p,fun,Sqh,rho,phi,th); 

fprintf('rel error for S (from another sph): %e\n',norm(Kq-Kq4)/nrmq)

Kq5 = Sh_Kernel_Eval_off(q,'SMat',1,1,rho,phi,th,Nr);

fprintf('rel error for S (Sh_Kernel_Eval_off): %e\n',norm(Kq-Kq5)/nrmq)

% Far evaluation for Sp
par.flag_pot='dSL_L_3D'; 
dKY = Kernel_Eval(Xtrg,X,par)*Y; 
dK2 = Seg(n)*(-n-1)*Y.*Rtrg.^(-n-2); 

fprintf('\n Check far eval for dSL[Y] outside sphere r=%1.1f \n',rt)
fprintf('error for Sp: %e\n',norm(dKY-dK2))

fprintf('\n Check far eval for dSL[q] outside sphere r=%1.1f \n',rt)
dKq = Kernel_Eval(Xtrg,X,par)*q; 
dKq2 = shSyn(((-nn-1).*rt.^(-nn-2)).*Sqh); %(all pts have same radius)

fprintf('rel error for Sp: %e\n',norm(dKq-dKq2)/nrmq)

dKq = Kernel_Eval(Xtrg2,X,par)*q; 
fun = @(r,n) (-n-1)*r.^(-n-2);
dKq3 = LOCAL_off_shSyn(p,fun,Sqh,Rtrg2,u,v); 

fprintf('rel error for Sp (local shSyn / normal perturbation): %e\n',norm(dKq-dKq3)/nrmq)
%{
figure; 
plot(real(dKq-dKq3))
figure; 
plot(imag(dKq-dKq3))
%}

%q=Y;
%qh=shAna(q); 
%Sqh = Seg(nn).*qh;

cB = [2 3 2.5]; cB = (3/norm(cB))*cB;  
Xtrg3 = X+repmat(cB,np,1);  
[th,phi,rho] = cart2sph(Xtrg3(:,1),Xtrg3(:,2),Xtrg3(:,3)); 
par.nor = Nr; 
dKq = Kernel_Eval(Xtrg3,X,par)*q;
th(th<0)=th(th<0)+2*pi; 
phi=pi/2-phi; 
dKq4 = LOCAL_off_shSyn(p,fun,Sqh,rho,phi,th); 

fprintf('rel error for Sp (from another sph): %e\n',norm(dKq-dKq4)/nrmq)
 
dKq5 = Sh_Kernel_Eval_off(q,'SpMat',1,1,rho,phi,th,Nr);

figure; 
plot(real(dKq))
hold on; 
plot(real(dKq5),'r')

fprintf('rel error for Sp (Sh_Kernel_Eval_off): %e\n',norm(dKq-dKq5)/nrmq)
fprintf('Dif between 2 approaches for (Sh_Kernel_Eval_off): %e\n',norm(dKq4-dKq5)/nrmq)


% Target points inside sphere
rt=0.25; 
Xtrg = rt*X;
Rtrg = sqrt(sum(Xtrg'.*Xtrg'))'; 

%Compare kernel eval with spherical harmonic expansion
par.flag_pot='SL_L_3D'; par.nor = X; 
Y = real(Y); 
KY = Kernel_Eval(Xtrg,X,par)*Y; 
K2 = Seg(n)*Y.*Rtrg.^n; 

fprintf('\n Check far eval for S[Y] inside sphere r=%1.2f \n',rt)
fprintf('error for S: %e\n',norm(KY-K2))

fprintf('\n Check far eval for S[q] inside sphere r=%1.2f \n',rt)
Kq = Kernel_Eval(Xtrg,X,par)*q; 
Kq2 = shSyn((rt.^nn).*Sqh); %(all pts have same radius)

fprintf('rel error for S: %e\n',norm(Kq-Kq2)/nrmq)

fun = @(r,n) r.^n; 
Xtrg2 = Xtrg+repmat(min(rt/2,1-rt/2)*(2*rand(np,1)-1),1,3).*X; 
Rtrg2 = sqrt(sum(Xtrg2'.*Xtrg2'))';
Kq = Kernel_Eval(Xtrg2,X,par)*q; 
Kq3 = LOCAL_off_shSyn(p,fun,Sqh,Rtrg2,u,v); 

fprintf('rel error for S (local shSyn / normal perturbation): %e\n',norm(Kq-Kq3)/nrmq)

Kq5 = Sh_Kernel_Eval_off(q,'SMat',1,0,Rtrg2,u,v);

fprintf('rel error for S (Sh_Kernel_Eval_off): %e\n',norm(Kq-Kq5)/nrmq)

% Far evaluation for Sp
par.flag_pot='dSL_L_3D'; 
dKY = Kernel_Eval(Xtrg,X,par)*Y; 
dK2 = Seg(n)*n*Y.*Rtrg.^(n-1); 

fprintf('\n Check far eval for dSL[Y] inside sphere r=%1.2f \n',rt)
fprintf('error for Sp: %e\n',norm(dKY-dK2))

fprintf('\n Check far eval for dSL[q] inside sphere r=%1.2f \n',rt)
dKq = Kernel_Eval(Xtrg,X,par)*q; 
dKq2 = shSyn((nn.*rt.^(nn-1)).*Sqh); %(all pts have same radius)

fprintf('rel error for Sp: %e\n',norm(dKq-dKq2)/nrmq)

dKq = Kernel_Eval(Xtrg2,X,par)*q; 
fun = @(r,n) n*r.^(n-1);
dKq3 = LOCAL_off_shSyn(p,fun,Sqh,Rtrg2,u,v); 

fprintf('rel error for Sp (local shSyn / normal perturbation): %e\n',norm(dKq-dKq3)/nrmq)

dKq5 = Sh_Kernel_Eval_off(q,'SpMat',1,0,Rtrg2,u,v);

fprintf('rel error for Sp (Sh_Kernel_Eval_off): %e\n',norm(dKq-dKq5)/nrmq)
%}
end

function Mq = LOCAL_off_shSyn(p,fun,qh,r,u,v)

Mq = zeros(size(r)); 

for n=0:p
   ind = n^2+1:(n+1)^2;  
   Yq = Ynm(n,[],u,v)*qh(ind,:);
   Mq = Mq+fun(r,n).*Yq; 
end

end
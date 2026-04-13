function Y = Sh_Mod_Kernel_Eval_off(phi,type,doAna,out,rho,u,v,Nr,lambda)
%{
This code evaluates the integral kernel given by type and density given by
phi on spherical coordinates (rho,u,v). 

Inputs
phi - (double (np / sp) x nv) integral density or corresponding vector
spherical harmonics coefficients 
type - kernel to be evaluated (SMat, SpMat, DMat, etc). 
doAna - (bool) indicates whether phi is density (np x nv) or spherical harmonic coefs (sp x nv)
out - (bool) exterior vs interior problem
(rho,u,v) - (double ntrg x 1) spherical coordinates arrays. Code evaluates
for rho=1, [u v] = gl_grid(p) if empty using fast transforms. 
Nr - (double ntrg x 3) Normal vector (if needed). 
lambda - (double) parameter for modified Laplace
Output 
Y - (double np x nv) integral kernel K[phi] evaluated at target points
%}
if(nargin==0), testShKEvaloff(); return;end

if(nargin<8)
    Nr=[]; 
end

LSeg  = @(n) (1./(2*n+1)); 

if size(phi,1)>1
% do Analysis to obtain shc or get from input
if doAna
    [np,d2] = size(phi); 
    p = (sqrt(2*np+1)-1)/2;
    sp= (p+1)^2; 
    
    % phi = Sum(sh_nm*Ynm(u,v))
    shc = shAna(phi);
else
    [sp,d2] = size(phi); 
    p = sqrt(sp)-1;    
    shc = phi; 
end

%nn = floor(sqrt((1:sp)'-1));
%Lam = repmat(LSeg(nn),1,d2); 

% density for S[phi]
Sqh = shc;

% Synthesize
if d2>1 || ~isempty(u)
    Y = LOCAL_off_shSyn(p,Sqh,out,type,rho,u,v,Nr,lambda); 
else
    Y = LOCAL_self_shSyn(p,Sqh,out,type,rho,lambda); 
end

else
   p = phi(1); d2 = phi(2);  
   if size(phi,2)>2
      Rdiag = phi(3:end).';  
   else
      Rdiag=[]; 
   end
   % Build dense matrix
   Y = LOCAL_shSyn_matrix(p,d2,out,type,rho,u,v,Nr,Rdiag,lambda); 
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Yq = LOCAL_self_shSyn(p,qh,out,type,r,lambda)

% get function 
%fun = LOCAL_get_Frn(type,out); 
Fun=Get_Spectra_Vec(p,lambda,r,type,out);
sp = (p+1)^2; 
ii = (1:sp)'; 
nn = floor(sqrt(ii-1)); 
Fr = Fun([nn+1])'; 

Mqh = repmat(Fr,1,size(qh,2)).*qh; 
Yq = shSyn(Mqh);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Yq = LOCAL_off_shSyn(p,qh,out,type,r,u,v,Nr,lambda)

d2 = size(qh,2);  
ntrg = length(u)/d2; 
Yq = zeros(ntrg,p+1); 

tic;

Fr=Get_Spectra_Vec(p,lambda,r',type,out);
   

if strcmp(type,'SMat') || strcmp(type,'DMat')

if d2==1
Yq2 = zeros(ntrg,p+1);     
for k=0:p
   ind = k^2+1:(k+1)^2;  
   Yq2(:,k+1) =  Fr(:,k+1).*(reshape(Ynm(k,[],u,v)*qh(ind,:),[],1));  
end
else
for k=0:p
   ind = k^2+1:(k+1)^2;
   FY = repmat(Fr(:,k+1),1,2*k+1).*Ynm(k,[],u,v); 
   FY = reshape(permute(reshape(FY,[ntrg d2 2*k+1]),[1 3 2]),ntrg,[]);
   Yq(:,k+1) =  FY*reshape(qh(ind,:),[],1);  
end
Yq2=Yq; 
end

else
    
% Normal Unit vector
er = [sin(u).*cos(v) sin(u).*sin(v) cos(u)];
dotp = @(x,y,nv) x(1:nv,:).*conj(y(1:nv,:)) +...
    x(nv+1:2*nv,:).*conj(y(nv+1:2*nv,:)) + x(2*nv+1:3*nv,:).*conj(y(2*nv+1:3*nv,:));


% er dot Nr(x) 
rdotN = dotp(er(:),Nr(:),d2*ntrg); 

SFr = Get_Spectra_Vec(p,lambda,r',[type(1) 'Mat'],out)./repmat(r+(r==0),1,p+1); 


if d2==1
Yq2 = zeros(ntrg,p+1);     
for k=0:p
   ind = k^2+1:(k+1)^2;
   GYk = Vnm('GY',[],k,[],u,v,er); 
   GYdotN = zeros(ntrg,size(GYk,2)); 
   for mk=1:size(GYk,2)
      GYdotN(:,mk) = dotp(GYk(:,mk),Nr(:),ntrg); 
   end

   Yq2(:,k+1) =  rdotN.*(Fr(:,k+1).*(reshape(Ynm(k,[],u,v)*qh(ind,:),[],1))) ...
      + SFr(:,k+1).*(reshape(GYdotN*qh(ind,:),[],1));  
end
else
for k=0:p
   ind = k^2+1:(k+1)^2;
   GYk = Vnm('GY',[],k,[],u,v,er); 
   GYdotN = zeros(d2*ntrg,size(GYk,2)); 
   for mk=1:size(GYk,2)
      GYdotN(:,mk) = dotp(GYk(:,mk),Nr(:),d2*ntrg); 
   end
   
   FY = repmat(Fr(:,k+1).*rdotN,1,2*k+1).*Ynm(k,[],u,v); 
   GY = repmat(SFr(:,k+1),1,2*k+1).*GYdotN; 
   
   for nd = 1:d2
       idY = (1:ntrg)+(nd-1)*ntrg; 
   Yq(:,k+1) = Yq(:,k+1) + (FY(idY,:)+GY(idY,:))*qh(ind,nd);
   end
end
Yq2=Yq; 
end    
end

Yq = reshape(sum(Yq2,2),ntrg,1);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K = LOCAL_shSyn_matrix(p,d2,out,type,r,u,v,Nr,Rdiag,lambda)

LSeg  = @(p) 1;

sp = (p+1)^2; 
np = 2*p*(p+1); 
ntrg = length(u)/d2; 
K = zeros(ntrg,sp*d2); 

Fr=Get_Spectra_Vec(p,lambda,r',type,out);

T = shAna(eye(np)); 

if strcmp(type,'SMat') || strcmp(type,'DMat')

if d2==1
  
for k=0:p
   ind = k^2+1:(k+1)^2;  
   K(:,ind) =  repmat(LSeg(k)*Fr(:,k+1),1,2*k+1).*Ynm(k,[],u,v);  
end

if ~isempty(Rdiag)
   T = diag(Rdiag)*T;  
end

K = K*T; 

else

for k=0:p
   ind = k^2+1:(k+1)^2;
   FY = repmat(LSeg(k)*Fr(:,k+1),1,2*k+1).*Ynm(k,[],u,v); 
   FY = reshape(permute(reshape(FY,[ntrg d2 2*k+1]),[1 3 2]),ntrg,[]);
   indrep = repmat(ind.',1,d2)+repmat((0:d2-1)*sp,2*k+1,1);
   indrep = indrep(:); 
   K(:,indrep) =  FY;  
end

tmp = K; 
K = zeros(ntrg,d2*np); 

for i=1:d2
   indt = (1:np)+(i-1)*np; 
   indf = (1:sp)+(i-1)*sp; 
   if ~isempty(Rdiag)
       K(:,indt) = tmp(:,indf)*diag(Rdiag(indf))*T;  
   else
       K(:,indt) = tmp(:,indf)*T; 
   end
end

end

else
    
% Normal Unit vector
er = [sin(u).*cos(v) sin(u).*sin(v) cos(u)];
dotp = @(x,y,nv) x(1:nv,:).*conj(y(1:nv,:)) +...
    x(nv+1:2*nv,:).*conj(y(nv+1:2*nv,:)) + x(2*nv+1:3*nv,:).*conj(y(2*nv+1:3*nv,:));


% er dot Nr(x) 
rdotN = dotp(er(:),Nr(:),d2*ntrg); 


SFr = Get_Spectra_Vec(p,lambda,r',[type(1) 'Mat'],out)./repmat(r+(r==0),1,p+1); 
    
if d2==1
K = zeros(ntrg,sp);     
for k=0:p
   ind = k^2+1:(k+1)^2;
   GYk = Vnm('GY',[],k,[],u,v,er); 
   GYdotN = zeros(ntrg,size(GYk,2)); 
   for mk=1:size(GYk,2)
      GYdotN(:,mk) = dotp(GYk(:,mk),Nr(:),ntrg); 
   end

   K(:,ind) =  repmat(LSeg(k)*rdotN.*Fr(:,k+1),1,2*k+1).*Ynm(k,[],u,v) ...
      + repmat(LSeg(k)*SFr(:,k+1),1,2*k+1).*GYdotN;  
end

K = K*T;

else
K = zeros(ntrg*d2,sp);      
for k=0:p
   ind = k^2+1:(k+1)^2;
   GYk = Vnm('GY',[],k,[],u,v,er); 
   GYdotN = zeros(d2*ntrg,size(GYk,2)); 
   for mk=1:size(GYk,2)
      GYdotN(:,mk) = dotp(GYk(:,mk),Nr(:),d2*ntrg); 
   end
   
   FY = repmat(LSeg(k)*Fr(:,k+1).*rdotN,1,2*k+1).*Ynm(k,[],u,v); 
   GY = repmat(LSeg(k)*SFr(:,k+1),1,2*k+1).*GYdotN; 
   K(:,ind) =  FY + GY; 
end
 
tmp = K; K = zeros(ntrg,np*d2); 
for i=1:d2
    ind = (1:np) + (i-1)*np; 
    K(:,ind) = tmp(ind,:)*T; 
end

end    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fun = LOCAL_get_Frn(type,out)
SYext = @(r,n) r.^(-n-1); 
SYint = @(r,n) r.^n; 
SpYext = @(r,n) -(n+1).*r.^(-n-2);
SpYint = @(r,n) n.*r.^(n-1); 
DYext = @(r,n) n.*r.^(-n-1);
DYint = @(r,n) -(n+1).*r.^n; 

if out
%Outside sphere
    switch type 
        case 'SMat'
            fun = SYext; 
        case 'SpMat'
            fun = SpYext;
        case 'DMat'
            fun =  DYext; 
    end
else
%Inside sphere
    switch type 
        case 'SMat'
            fun = SYint;
        case 'SpMat'
            fun = SpYint;
        case 'DMat'
            fun = DYint; 
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testShKEvaloff()
p=10; 
Sc = SurfaceSph(shape_gallery(p,''));
X = reshape(Sc.cart.to_array,[],3); 
np = 2*p*(p+1);  
par = Kernel_Eval_parameters('SL_L_3D',0,1,1,1,1e-8,2,400,1);
par.dim = 3;
[~, gwt]=g_grid(p+1);
wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
wt = wt(:);
W = Sc.geoProp.W; W= W.*wt; 
Nr = reshape(Sc.geoProp.nor.to_array,[],3); 
par.X = X; par.nor = Nr; par.W2 = W.'; 

dpar = par; dpar.flag_pot='dSL_L_3D'; 

% Random density of freq p-1
%qh = [rand(p^2,1); zeros((p+1)^2-p^2,1)]; 
qh = zeros((p+1)^2,1); 
qh(2)=1; 
q = shSyn(qh); 
nrmq=norm(q); 

% exterior examples for S and Sp
fprintf('\n Outside the sphere \n')
out=true; 

cB = [2 2.5 2]; dst = 2;  
cB = (dst+2)/(norm(cB))*cB; 
Xtrg = X+repmat(cB,np,1); 
% Find rho, u,v
[v,u,rho] = cart2sph(Xtrg(:,1),Xtrg(:,2),Xtrg(:,3)); 
v(v<0)=v(v<0)+2*pi; u=pi/2-u; 
dpar.nor = Xtrg./repmat(rho,1,3); 

% Check with regular Kernel Eval
Kq = Kernel_Eval(Xtrg,X,par)*q;
dKq = Kernel_Eval(Xtrg,X,dpar)*q;

ShKq = Sh_Mod_Kernel_Eval_off(qh,'SMat',0,out,rho,u,v);
ShdKq = Sh_Mod_Kernel_Eval_off(qh,'SpMat',0,out,rho,u,v);
ShdKq2 = Sh_Mod_Kernel_Eval_off(qh,'SpMat',0,out,rho,u,v,dpar.nor);

fprintf('Check w/ Nr = er: %e\n',norm(ShdKq2-ShdKq)/nrmq)

fprintf('rel error for S outside: %e\n',norm(Kq-ShKq)/nrmq)
fprintf('rel error for Sp outside: %e\n',norm(dKq-ShdKq)/nrmq)

dotp = @(x,y) conj(x(1:np,:)).*y(1:np,:)+conj(x(np+1:2*np,:)).*y(np+1:2*np,:)+conj(x(2*np+1:3*np,:)).*y(2*np+1:3*np,:);
[u0,v0] = gl_grid(p); 
GY = Vnm('GY',[],1,-1,u,v,ones(size(u)));
GY = reshape(GY,[],3); 
NrG = sqrt(sum(GY*GY',2)); 
GY = GY./repmat(NrG,1,3); 
dpar.nor = reshape(GY,[],3); 

er = [sin(u).*cos(v) sin(u).*sin(v) cos(u)];
rdotN = dotp(er(:),dpar.nor(:));
max(abs(rdotN))

dKq = Kernel_Eval(Xtrg,X,dpar)*q;
ShdKq = Sh_Mod_Kernel_Eval_off(qh,'SpMat',0,out,rho,u,v,dpar.nor);
fprintf('rel error for Sp outside (Nr): %e\n',norm(dKq-ShdKq)/nrmq)

plot(real(dKq)); hold on; plot(real(ShKq),'r')
figure; 
GdotN = dotp(GY(:),dpar.nor(:)); 
plot(abs(dKq./ShdKq),'-or')
figure; 
plot(abs(ShdKq./GdotN),'-o')

LSeg  = @(n) (1./(2*n+1)); 
SYext = @(r,n) LSeg(n)*r.^(-n-1); 
ShdKq2 = SYext(rho,1).*NrG.*GdotN; 
plot(abs(ShdKq-ShdKq2))

fprintf('\n Inside the sphere \n')
out=false; 
Xtrg = 0.2*X+repmat([0.1 -0.1 0.1],np,1); 
% Find rho, u,v
[v,u,rho] = cart2sph(Xtrg(:,1),Xtrg(:,2),Xtrg(:,3)); 

v(v<0)=v(v<0)+2*pi; u=pi/2-u; 

dpar.nor = Xtrg./repmat(rho,1,3); 

% Check with regular Kernel Eval
Kq = Kernel_Eval(Xtrg,X,par)*q;
dKq = Kernel_Eval(Xtrg,X,dpar)*q;

ShKq = Sh_Mod_Kernel_Eval_off(q,'SMat',1,out,rho,u,v);
ShdKq = Sh_Mod_Kernel_Eval_off(q,'SpMat',1,out,rho,u,v);
ShdKq2 = Sh_Mod_Kernel_Eval_off(q,'SpMat',1,out,rho,u,v,dpar.nor);

fprintf('Check w/ Nr = er: %e\n',norm(ShdKq2-ShdKq)/nrmq)

fprintf('rel error for S outside: %e\n',norm(Kq-ShKq)/nrmq)
fprintf('rel error for Sp outside: %e\n',norm(dKq-ShdKq)/nrmq)

dpar.nor = Nr; 
dKq = Kernel_Eval(Xtrg,X,dpar)*q;
ShdKq = Sh_Mod_Kernel_Eval_off(q,'SpMat',1,out,rho,u,v,Nr);
fprintf('rel error for Sp outside (Nr): %e\n',norm(dKq-ShdKq)/nrmq)

end
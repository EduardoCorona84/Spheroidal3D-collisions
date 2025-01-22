function [Y,Sc,DMV,D] = Shg_MatVec_RBS(Sc,V,DMV,C,p,n3,flag_pot,type)

if nargin<8
    type='dense'; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Surface parameters: 

% Number of points on each surface
np = 2*p*(p+1); nc = 1;  
if ~iscell(Sc) 
    Shape = Sc;
    Sc = cell(p,1);  
    % Model surface (same shape)
    Sc{p} = SurfaceSph(shape_gallery(p,Shape)); 
else
    if ischar(Sc{1})
       Shape = Sc;
       Sc = cell(p,n3); 
       for j=1:n3
           Sc{p,j} = SurfaceSph(shape_gallery(p,Shape{j}));
       end
       
    end

    nc = size(Sc,2);    
end

% Smooth Quadrature Weights (GL x Trapezoidal)
[~, gwt]=g_grid(p+1);
wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
wt = wt(:);

a = 0.5; 
shOrder  = @(x) (sqrt(1+2*size(x,1)/3)-1)/2;    
evc = @(c) c{1}; 

buildMV = isempty(DMV); 

X = cell(nc,1); W=X; Nr=X; D = X; D2=D; 
if buildMV
    DMV = X; 
end

% Surface info
for j=1:nc
% Points on Sc
X{j} = reshape(Sc{p,j}.cart.to_array,[],3); 
% Area element
W{j} = Sc{p,j}.geoProp.W; 
W{j} = W{j}.*wt; 
% Normal vector
Nr{j} = reshape(Sc{p,j}.geoProp.nor.to_array,[],3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If DMV is empty, generate block-diag MatVec using spherical harmonics
if buildMV
% Traction kernel 
%TMV = @(x) a*x+evc(kerneldS({x},Sc{shOrder(x)})); 
if strcmp(flag_pot,'DT_Stk_3D')
DMV{j} = @(x) a*x+evc(kerneldS({x},Sc{shOrder(x),j},[],1)); 
%[Ck,Bk,Dk] = LOCAL_BuildC(W{j},X{j},[],np);      
     
if strcmp(type,'dense')
   D{j} = DMV{j}(eye(3*np)); %+Bk'*Ck;       
   DMV{j} = @(x) D{j}*x;   
else
   D{j} = [];  
end

elseif strcmp(flag_pot,'SL_Stk_3D')
DMV{j} = @(x) evc(kernelS({x},Sc{shOrder(x),j},[],1));  
     
if strcmp(type,'dense')
   D{j} = DMV{j}(eye(3*np)); 
   DMV{j} = @(x) D{j}*x;
else
   D{j} = [];  
end    
    
end

% Kernel Eval (to subtract block-diag contribution from FMM)
%fprintf('\n Kernel Eval \n') 
params = Kernel_Eval_parameters(flag_pot,0,1,1,1,1e-8,2,400,0);
params.dim = 3; params.mu=1;   
Xv = reshape(repmat(X{j},1,3)',3,[])'; 
Nrv = reshape(repmat(Nr{j},1,3)',3,[])';    
Wv = repmat(W{j},1,3)'; Wv = Wv(:);    
params.X = Xv; params.nor = Nrv; params.W2 = Wv.'; 
params.ci = repmat((1:3)',np,1); params.cj = params.ci; 

K = Kernel_Eval(Xv,Xv,params);   

if strcmp('flag_pot','DT_Stk_3D')
if strcmp(type,'mfree')
   D{j}   = @(x) a*x+evc(kerneldS({x},Sc{shOrder(x),j})); % + Bk'*(Ck*x); 
   DMV{j} = @(x) a*x+evc(kerneldS({x},Sc{shOrder(x),j})) - K*x; % + Bk'*(Ck*x)-K*x;     
else 
   D2{j} = D{j}-K;     
   DMV{j} = @(x) D2{j}*x;    
end
else
if strcmp(type,'mfree')
   D{j}   = @(x) evc(kernelS({x},Sc{shOrder(x),j})); 
   DMV{j} = @(x) evc(kernelS({x},Sc{shOrder(x),j}))-K*x;     
else 
   D2{j} = D{j}-K;   
   DMV{j} = @(x) D2{j}*x;    
end    
end

if nc==1
    D = D{1}; 
    DMV = DMV{1}; 
end   
else
   D{j} = [];
   
end 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(V) 
    Y = V;  
else
% Build 3D array of n^3 objects
%Centers (re-order w/ morton-Z for TT)
% Build 3D array of n^3 objects (Change to Z ordering)
L = log2(n3); 

if isempty(C)
I = dec2bin((0:n3-1).',L); 
Xi = zeros(n3,3); 
for d=1:3
    Xi(:,d) = bin2dec(I(:,d:3:L)); 
end
    
% Centers (morton-Z order for TT)
dst=4; 
C = 2+dst*Xi;   
end

Cg = reshape(repmat(C.',np,1),3,[]).';
% Translate each copy and stack
if ~iscell(DMV)
   Xg = repmat(X{1},n3,1) + Cg;
   Wg = repmat(W{1},n3,1); 
   Nrg = repmat(Nr{1},n3,1);  
else
   Xg = zeros(np*n3,3); Wg=zeros(np*n3,1); Nrg=Xg; Ydg = zeros(3*np*n3,size(V,2));    
   for j=1:nc
      indx=(1:np)+np*(j-1); 
      ind =(1:3*np)+3*np*(j-1);
      Xg(indx,:) = X{j}+Cg(indx,:);  
      Wg(indx,:) = W{j}; 
      Nrg(indx,:) = Nr{j}; 
      
      Ydg(ind,:) = DMV{j}(V(ind,:)); 
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
% FMM dS/dn_x Single layer (?)
fprintf('\n Kernel Eval with STFMM3D \n')

src = X.';

v = zeros(3,np); 
% Add quadrature weights to density
v(1)=params.W2(1)*1;

% Single layer (gradient)
ifsingle=1; 
Ntot = np;   
sigma_sl = v; 
% No Dipoles (Double Layer) 
ifdouble=0;
sigma_dl = zeros(3,Ntot); 
sigma_dv = zeros(3,Ntot); 

ifpot  = 1; 
ifgrad = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Off-d Smooth integrals (FMM)
%fprintf('\n FMM apply \n') 
% Test fast matrix apply (lfmm3d)
% Single Layer (potential and gradient) 
iprec=3; %(9 digits. 4 for 12 digits)        
U=stfmm3dpart(iprec,Ntot,src,ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,ifpot,ifgrad);

PotU = U.pot; %3xN

% Double Layer (potential) 
ifsingle=0; ifdouble=1; sigma_dl=sigma_sl; sigma_dv = Nr.'; ifgrad=0; 
UD=stfmm3dpart(iprec,Ntot,src,ifsingle,sigma_sl,ifdouble,sigma_dl,sigma_dv,ifpot,ifgrad);

PotD = UD.pot; %3xN

% Vectorized form of this: 
%for i=1:Ntot
%   Ug(3*(i-1)+1:3*i) = GradU(:,:,i)*Nr(i,:)';    
%end
GradU = permute((1/4/pi)*U.grad,[2,1,3]); %3x3xN, GU(i,:,:)=U_{x_i}    
Ug = sum(reshape(reshape(GradU,9,[]).*repmat(Nr.',3,1),3,[])); %GU dot N(x)
Ug = Ug(:);      

Up = (1/4/pi)*PotU(:); 
Dp = (-1/4/pi)*PotD(:);      

display(norm(KdS(:,1)-Ug))
display(norm(KS(:,1)-Up))
display(norm(KD(:,1)-Dp))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

Ntot=3*np*n3; 

% We perform off-diagonal kernel eval densely for now
paroffd = Kernel_Eval_parameters(flag_pot,0,1,1,1,1e-8,2,400,0);
paroffd.dim = 3; paroffd.mu=1;     
Xv = reshape(repmat(Xg,1,3)',3,[])'; 
Nrv = reshape(repmat(Nrg,1,3)',3,[])';    
Wv = repmat(Wg,1,3)'; Wv = Wv(:);    
paroffd.X = Xv; paroffd.nor = Nrv; paroffd.W2 = Wv.'; 
paroffd.ci = repmat((1:3)',n3*np,1); paroffd.cj = paroffd.ci; 

Y = Kernel_Eval(Xv,Xv,paroffd)*V;  

% Apply: 
%fprintf('\n spharm diag (aI-D) + off diag FMM  \n') 
if ~iscell(DMV)
    V = reshape(V,3*np,[]);     
    Y = real(reshape(DMV(V),Ntot,[])+Y);      
else
    Y = real(Ydg+Y);    
end

end
%display(toc);      
end

function [C,B,D] = LOCAL_BuildC(W,X,Xc,np)

if ~isempty(Xc)
% Center X
X = X - repmat(Xc,3*np,1); 
end

% C*sigma = [int{sigma} ; int{X \times sigma}]
C = zeros(6,3*np);  
indx =1:3:3*np; indy =2:3:3*np; indz=3:3:3*np;  
% First three vectors are just integrals of f for each coordinate
C(1,indx) = W; 
C(2,indy) = W; 
C(3,indz) = W;
% Last three are integrals of T = (X-X_c) (x) f
C(4,indy) = -W.*X(:,3); C(4,indz) = W.*X(:,2); 
C(5,indz) = -W.*X(:,1); C(5,indx) = W.*X(:,3);
C(6,indx) = -W.*X(:,2); C(6,indy) = W.*X(:,1); 

% Basis for nullspace
% First 3: F_i / |Sc|
D = C;  
D(1,indx) = ones(size(W)); 
D(2,indy) = ones(size(W)); 
D(3,indz) = ones(size(W));
% Last three are (X-X_c) (x) 1
D(4,indy) = -X(:,3); D(4,indz) = X(:,2); 
D(5,indz) = -X(:,1); D(5,indx) = X(:,3);
D(6,indx) = -X(:,2); D(6,indy) = X(:,1); 

% We divide by the integral of ((X-X_c) (x) 1)^2 over Sc. 
W2 = zeros(1,3*np); 
W2(indx)=W; W2(indy)=W; W2(indz)=W; 

B=D; 
for i=1:3
   B(i,:) = B(i,:)./sum(W);    
end
for i=4:6
   B(i,:) = B(i,:)./sum(W2.*B(i,:).^2);       
end 
end

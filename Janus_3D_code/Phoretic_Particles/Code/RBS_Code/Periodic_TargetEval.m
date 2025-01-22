function [Y,params] = Periodic_TargetEval(p,C,r,a,b,flag_pot,dense,neigh,Xtrg,Nrtrg,Q)
%{
Given Gamma the union of boundaries of spheres with radii r(i) and centers
C(:,i) in [a,b]^3, we first produce Gamma_p, by adding 26 periodic copies
in neighboring boxes of the same diameter. This code then evaluates the 
integral equation indicated by flag_pot with density Q: 

Y = S_{Gamma_p}[Q](x) = int_{Gamma_p} K(x,y)Q(y)dSy    ;  x in Xtrg

by discretizing Q using spherical harmonics of degree p. 

INPUTS: 
p - (int) - spherical harmonic degree
C - (double n_c x 3) - centers in original box [a,b]^3
r -  (double n_c x 1) - sphere radii
a,b - (doubles) box start and finish
flag_pot - (string) - the following cases are supported: 
 
'SL_L_3D'     single layer, Laplace
'dSL_L_3D'    normal derivative of single layer, Laplace
'DL_L_3D'     double layer, Laplace
'SL_Stk_3D'   single layer, Stokes
'DL_Stk_3D'   double layer, Stokes
'dSL_Stk_3D'  normal derivative of single layer, Stokes
'TSL_Stk_3D'  traction kernel of single layer, Stokes

dense - (bool) dense apply vs FMM apply for far interactions
Xtrg - (double Ntrg x 3) - target points
Ntrg - (double Ntrg x 3) - normal vector at targets (if required, otherwise empty)
Q - (double n_c*N_deg x 1) array of density or densities 

%}

% Get parameter array, containing centers, radii and geometric data for neighbor copies
params = LOCAL_setparams(p,C,r,a,b,flag_pot,dense,neigh,Xtrg,Nrtrg); 

if isempty(Q)
    Y = VSh_MatVec_RB_trg('Mat',Xtrg,Nrtrg,params);  
else
   % Eval routine (w/ spherical harmonics for near and self-eval)
    Y = VSh_MatVec_RB_trg(Q,Xtrg,Nrtrg,params);  
end


end

function params = LOCAL_setparams(p,C,r,a,b,flag_pot,dense,neigh,Xtrg,Nrtrg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if neigh
%Add 26 periodic copies for neighboring boxes to [a,b]^3
diam = b-a; 
nc = size(C,1); 
[xx,yy,zz] = meshgrid([-diam 0 diam]); 
Cadd = repmat([xx(:) yy(:) zz(:)].',nc,1); 

C = repmat(C,27,1)+reshape(Cadd,3,[]).'; %new center array in [a-diam,b+diam]^3
r = repmat(r,27,1); %new raddi array size 27*nc x 1
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%kernel dimension (3 or 1) 
kerd = 3; 
if strcmp(flag_pot(4),'L') || strcmp(flag_pot(5),'L')
    kerd=1; 
end

% Points on each swimmer
np = 2*p*(p+1); 
% Number of rigid bodies
nc=size(C,1); 
% Model Surfaces
Sc = SurfaceSph(shape_gallery(p,'')); 
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
    'keval','points','dim',3,'mu',1,'Xp',Xg,'X',Xv,'nor',Nrv,'W2',Wv.',...
    'p',p,'kerd',kerd,'n3',nc,'rd',r(:),'dense',dense,'C',C,'out',1,...
    'mdist',4,'a',0); 

if kerd>1
params.ci = repmat((1:kerd)',size(Xtrg,1),1); 
params.cj = repmat((1:kerd)',nc*np,1); 
end

% Find neighbors for near eval
params = LOCAL_findneighbors(params); 

end

function params = LOCAL_findneighbors(params)
% Compute distances and neighboring spheres
distC = LOCAL_CenterDistance(params.C); 
nc = size(params.C,1); 

neigh = cell(params.n3,1);  
for i=1:nc
   neigh{i} = find(distC(i,:)<params.rd(i)*params.mdist); 
end

params.neigh=neigh; 
end

function den = LOCAL_CenterDistance(C)
[Y_g1,  X_g1  ] = meshgrid(C(:,1), C(:,1));
[Y_g2,  X_g2  ] = meshgrid(C(:,2), C(:,2));
[Y_g3,  X_g3  ] = meshgrid(C(:,3), C(:,3));
d1 = (X_g1 - Y_g1); d2 = (X_g2 - Y_g2); d3 = (X_g3 - Y_g3); 
den = sqrt(d1.^2 + d2.^2 + d3.^2); 
end
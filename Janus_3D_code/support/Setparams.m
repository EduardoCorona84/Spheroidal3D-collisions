function params = Setparams(p,C,r,flag_pot,Shape,dense,doAna,mdist,out,lambda)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%kernel dimension (3 or 1) 


if(strcmp(flag_pot,'SL_Stk_3D'))
kerd = 3;
else
    kerd = 1;
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
    'mdist',mdist,'a',0,'doAna',doAna,'lambda',lambda); 

if kerd>1
params.ci = repmat((1:kerd)',nc*np,1); params.cj = params.ci; 
end

end


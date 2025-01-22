function params = RBS_set_params(p,C,r,flag_pot,Shape,Sc,kerd,dense,doAna,mdist,eps,out) 

np=2*p*(p+1); 
Nb = kerd*np; % DOF per particle   
n3 = size(C,1); % no of particles
N = kerd*np*n3; %DOF total
r = r(:);

% Smooth Quadrature Weights (GL x Trapezoidal)
[~, gwt]=g_grid(p+1);
wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
wt = wt(:);

X = cell(n3,1); W=X; Nr=X;
Cg = reshape(repmat(C.',np,1),3,[]).';
rg = repmat(reshape(repmat(r.',np,1),1,[]).',1,3);
Xg = zeros(np*n3,3); Wg=zeros(np*n3,1); Nrg=Xg; Xrp = Xg;  

if size(Sc,2)>1
    tau=X; 
    
for j=1:n3
    % Points on Sc
    X{j} = reshape(Sc{p,j}.cart.to_array,[],3); 
    % Area element
    W{j} = Sc{p,j}.geoProp.W; 
    W{j} = W{j}.*wt; 
    % Normal vector
    Nr{j} = reshape(Sc{p,j}.geoProp.nor.to_array,[],3);
    
    % moment of inertia matrix
    stau = sum(sum(repmat(W{j},1,3).*X{j}.^2)); 
    WX = repmat(W{j},1,3).*X{j}; 
    WXtau = reshape(repmat(WX,3,1),[],9); 
    Xtau=repmat(X{j},1,3); 
    tau{j} = stau*eye(3) - reshape(sum(WXtau.*Xtau),3,3);

    indx=(1:np)+np*(j-1);
    Xrp(indx,:) = r(j)*X{j};
    Xg(indx,:) = Xrp(indx,:)+Cg(indx,:);  
    Wg(indx,:) = W{j}; 
    Nrg(indx,:) = Nr{j}; 
end
else
   X{1} = reshape(Sc{p}.cart.to_array,[],3); 
   W{1} = Sc{p}.geoProp.W; W{1} = W{1}.*wt; 
   Nr{1} = reshape(Sc{p}.geoProp.nor.to_array,[],3);
   
   % moment of inertia matrix
   stau = sum(sum(repmat(W{1},1,3).*X{1}.^2)); 
   WX = repmat(W{1},1,3).*X{1}; 
   WXtau = reshape(repmat(WX,3,1),[],9); 
   Xtau=repmat(X{1},1,3); 
   tau = stau*eye(3) - reshape(sum(WXtau.*Xtau),3,3); 
   
   Xrp = rg.*repmat(X{1},n3,1);
   Xg = Xrp + Cg;
   Wg = (rg(:,1).^2).*repmat(W{1},n3,1); 
   Nrg = repmat(Nr{1},n3,1); 
   W = W{1}; 
end

Xv = reshape(repmat(Xg,1,kerd)',3,[])'; 
Nrv = reshape(repmat(Nrg,1,kerd)',3,[])';    
Wv = repmat(Wg,1,kerd)'; Wv = Wv(:);  
 
diam = repmat(r,1,n3)+repmat(r.',n3,1); %diam(i,j) = r_i + r_j
mxrd = max(repmat(r,1,n3),repmat(r.',n3,1)); %max(r_i,r_j)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Form and fill params struct
params = struct('flag_pot',flag_pot,'kh',0,'transinv',1,'sym',1,'proxy',0,...
    'keval','points','dim',3,'mu',1,'Xrp',Xrp,'Xp',Xg,'Nrp',Nrg,'X',Xv,'nor',Nrv,'targnor',Nrv,'Wg',Wg,'W2',Wv.',...
    'p',p,'np',np,'kerd',kerd,'n3',n3,'Nb',Nb,'N',N,'rd',r,'diam',diam,'mxrd',mxrd,...
    'dense',dense,'C',C,'out',out,'mdist',mdist,'eps',eps,'a',0,'doAna',doAna); 

params.Shape=Shape; 
params.Sc = Sc; 
params.W = W; 
params.tau = tau; 

if kerd>1
    ctmp = repmat((1:kerd)',n3*np,1);
    params.ci = ctmp; 
    params.cj = params.ci; 
end
end


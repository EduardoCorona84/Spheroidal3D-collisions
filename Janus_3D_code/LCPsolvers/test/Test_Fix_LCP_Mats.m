%% Init params 
p = 3;
lambda=0.1; rd=1; n=2; Cdst=2.3;  ep=.3; Nt=500; dt=.1; tdisc='euler'; saveLCPs=true;  
initMode='vesicle';  tol=1e-4; mdist=3;  denseMV=true;  denseforce=1; gamma=1; 
boundary_label =  @(X,y) 0.5*X*y'./sqrt(sum(X.^2,2)).^2 + 1/2;
Fparams = struct('Nt',Nt,'dt',dt,'comp',1,'type','JanusAmp',...
    'lambda',lambda,'gamma',gamma,'denseMV',denseMV,...
    'typeMV','Vsh','tdisc',tdisc, ...
    'boundary_label',boundary_label,'denseforce',denseforce);
[C, init_dir] = init_lattice(2, 2.3); 
nC=size(Ct,1); 
% body parameters
n3 = size(C,1); 
rd=rd*ones(n3,1);
Fparams.parbd = struct('Shape','','n3',n3,'rd',rd,'diam',2*rd,'p',p,'mdist',mdist,'mxrd',rd(1),'eps',ep,'out',1);
Fparams.parbd.Ct = C;
Fparams.init_dir = init_dir;
% LCP solver parameters
Fparams.lcpOpts = struct(...
    'solver','proxquasinewton',...
    'max_iter',100,...
    'tol_rel',1e-6,...
    'tol_abs',1e-5, ...
    'stepSize',struct(...
        'init','uniform',...
        'kappa','uniform',...
        'eta','opt')...
);
% linear solver parameters
Fparams.parslv = struct('solver','gmres','tol',tol,'maxit',200,'rst',4,'prtype','bkdiag','prec',[],'prLCP',false); 
Fparams = RBS_Initialize_params(Fparams);
%% Initial config
Ct = C;
np = 2*p.*(p+1); 
%% Build the velocities the way Corona did for testing
Sc = SurfaceSph(shape_gallery(p,''));
Xp = reshape(Sc.cart.to_array,[],3);
Xg = repmat(Xp,nC,1) + reshape(repmat(Ct',np,1),3,[])';
%% W whatever that is... I think it is something about surface area of the sphere
[~, gwt]=g_grid(p+1);
wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
wt = wt(:);
WT = Sc.geoProp.W; WT= WT.*wt;
Wg = repmat(WT(:),nC,1);
%% Build Sparse representations of these matrices
IC=[];
JC=[];
IT=[];
JT=[];
VT=[];
VA =[];
VC =[];
for i=1:n3
idx = (1:np)+np*(i-1); 
X = Xg(idx,:);
W = Wg(idx); 
oW = ones(size(W));
sumW = sum(W); 
tau1 = sum(W.*X(:,3).^2)+sum(W.*X(:,2).^2); 
tau2 = sum(W.*X(:,1).^2)+sum(W.*X(:,3).^2); 
tau3 = sum(W.*X(:,2).^2)+sum(W.*X(:,3).^2); 
Wav = (1/sum(W))*W; 

VC = [VC; W; W; W ; ...
    -W.*X(:,3); -W.*X(:,1); -W.*X(:,2); ...
    W.*X(:,2);  W.*X(:,3); W.*X(:,1)];
VT = [VT; (1/sum(W)); (1/sum(W)); (1/sum(W)); ...
    (1/tau1); (1/tau2); (1/tau3)];
VA = [VA; Wav; Wav; Wav ; ...
    -(1/tau1)*W.*X(:,3); -(1/tau2)*W.*X(:,1); -(1/tau3)*W.*X(:,2);...
    (1/tau1)*W.*X(:,2);  (1/tau2)*W.*X(:,3); (1/tau3)*W.*X(:,1)];


% norm(VT.*VC - VA) / norm(VA)

indx =(1:3:3*np)+3*np*(i-1); 
indy =(2:3:3*np)+3*np*(i-1); 
indz =(3:3:3*np)+3*np*(i-1);

% C*sigma = [int{sigma} ; int{X \times sigma}]  
IC = [IC; reshape(repmat(6*(i-1)+(1:6),np,1),[],1) ; reshape(repmat(6*(i-1)+(4:6),np,1),[],1)];
JC = [JC; indx'; indy'; indz'; indy'; indz'; indx'; indz'; indx'; indy'];
IT = [IT; (6*(i-1)+(1:6))'];
JT = [JT; (6*(i-1)+(1:6))']; 

end
N = 3*np*n3;
T = diag(VT);
C = sparse(IC,JC,VC,6*n3,N); 
A = sparse(IC,JC,VA,6*n3,N); 

norm(A-T*C,'fro') / norm(A, 'fro')


%%
[A,b,F,Abad] = getLCPfromC(Ct, [], p, Fparams);

Acor = F'*T*(F' \ Abad);
norm(A - Acor) / norm(A)
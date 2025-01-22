function [Y,Itrg,neigh] = VSh_MatVec_RB_trg_old(V,Xtrg,Nrtrg,params)
%{
This code evaluates the integral equation indicated by flag_pot in the params 
struct with density V: 

Y = S_{Gamma}[Q](x) = int_{Gamma} K(x,y)V(y)dSy    ;  x in Xtrg

Where Gamma is the union of boundaries of spheres with radii r(i) and centers
C(:,i), by discretizing Q using spherical harmonics of degree p. 

INPUTS: 
params - parameter struct with fields such as: 
    p - (int) - spherical harmonic degree
    n3 - (int) - number of objects 
    kerd - (int) - kernel dimension 
    dense - (bool) - Dense vs FMM far interactions
    out - (bool) exterior or interior problem
    C - (double n_c x 3) - centers in original box [a,b]^3
    rd -  (double n_c x 1) - sphere radii
    flag_pot - (string) - the following cases are supported: 
        'SL_L_3D'     single layer, Laplace
        'dSL_L_3D'    normal derivative of single layer, Laplace
        'DL_L_3D'     double layer, Laplace
        'SL_Stk_3D'   single layer, Stokes
        'DL_Stk_3D'   double layer, Stokes
        'dSL_Stk_3D'  normal derivative of single layer, Stokes
        'TSL_Stk_3D'  traction kernel of single layer, Stokes

    dense - (bool) dense apply vs FMM apply for far interactions
V - (double n_c*N_deg x 1) array of density or densities 
L - (sparse array) Extra matrix for completion flow / nullspace correction
Xtrg - (Ntrg x 3) Target point array
Nrtrg - (Ntrg x 3) Normal vector at targets (if needed, empty otherwise) 

%}

%TODO: add option to evaluate multiple vectors (nV = size(V,2);)
persistent ShTg ShTr ShY Sc

p  = params.p;      % Spharm degree p 
n3 = params.n3;     % Number of objects
kerd = params.kerd; % Kernel dimension
dense=params.dense; % Dense vs FMM off diagonal 
C = params.C;       % object centers
out = params.out;   % outside vs inside sphere
mdist = params.mdist; %near field is dist < mdist*rd(i)
if isfield(params,'MRot')
   MRot = params.MRot;
   rot = true; 
   % Xt{k} = X_0*MRot{k}
else
   MRot = cell(1,n3);  
   rot = false; 
end
% Radius of spheres
if isfield(params,'rd')
   % make this n3 x 1
   rd=params.rd;
   rda = rd.^(strcmp(params.flag_pot(1:3),'SL_')); 
   rdif1=true; 
else
   rd=ones(n3,1); rda=rd; rdif1=false; 
end

pot=params.flag_pot;
nortrg=true; 
Nrtrgv = reshape(repmat(Nrtrg,1,kerd)',3,[])';    

switch pot(1:3)
case 'SL_'
   pMat = 'SMat';  
   nortrg=false; 
case 'dSL'
   pMat = 'SpMat'; 
case 'TSL' 
   pMat='TSMat'; 
   %Coefficient matrices (add option to load only once) 
    if isempty(ShTg) || size(ShTg{1},1)/3 ~= (p+1)^2
        fprintf('\n Computing Coeff matrices for T \n')
        [ShTg,ShTr,ShY] = Surfgrad_coeffs(ceil(1.5*p),p,1,'VW','VW',0); 
    end
case 'DL_'
   pMat='DMat';
   nortrg=false; 
case 'dDL_'
   pMat='DpMat';
   nortrg=false; 
case 'TDL'
   pMat='TDMat'; 
   %Coefficient matrices (add option to load only once) 
    if isempty(ShTg) || size(ShTg{1},1)/3 ~= (p+1)^2
        fprintf('\n Computing Coeff matrices for T \n')
        [ShTg,ShTr,ShY] = Surfgrad_coeffs(ceil(1.5*p),p,1,'VW','VW',0); 
    end
    nortrg=false; 
end

np = 2*p*(p+1); 
Nb = kerd*np; 
N = Nb*n3;  
X = params.Xp;
Xv = params.X; 
% Target pts and Normal vector for tensor kernels
Xtrgv = reshape(repmat(Xtrg,1,kerd)',3,[])';

% Determine near field / neighbors
%distC = LOCAL_CenterDistance(Xtrg,C); 
[Itrg,neigh] = LOCAL_find_neigh(Xtrg,C,rd,mdist); 

% # of target points
Ntrg = kerd*size(Xtrg,1);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Off-d + diag contributions

% If V is not empty, compute A*V. Otherwise, construct A densely. 
if ~isempty(V) && isnumeric(V)  
    V = reshape(V,Nb,[]);
    
    % Compute spherical harmonic coefficients
    if kerd==3
       V2 = [V(1:3:end,:);V(2:3:end,:);V(3:3:end,:)];  
       if rot
           V2 = reshape(V2,np,[]); 
           for k=1:n3
              cols = (1:3)+3*(n3-1); 
              V2 = V2(:,cols)*MRot{k}';  
           end
           V2 = reshape(V2,Nb,[]); 
       end
       Vh = VshAna(V2,'VW'); 
    else
       Vh = shAna(V);  
    end
    
    V = reshape(V,N,[]); 
    Y = zeros(Ntrg,size(V,2)); 
    
    for nbox=1:n3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Near interactions
        I_box = (1:Nb)+Nb*(nbox-1); 
        % Build neighbor index
        neigh{nbox} = reshape(neigh{nbox},1,[]); 
        
        if sum(neigh{nbox})>0
        I_neigh   = repmat(neigh{nbox},kerd,1); 
        I_neigh = I_neigh(:); 
        Nnear = sum(neigh{nbox}); 
        
         % near target points relative to box center
        Xbtrg = (1/rd(nbox))*(Xtrg(neigh{nbox},:) - repmat(C(nbox,:),Nnear,1)); 
        Nrtrgb = Nrtrg(neigh{nbox},:); 
        Nrtrgb = reshape(repmat(Nrtrgb,1,kerd)',3,[])';    
        
        % rotate back to reference
        if rot
           Xbtrg = Xbtrg*MRot{nbox}'; 
        end
        
        % spherical coordinates
        [th,phi,rho] = cart2sph(Xbtrg(:,1),Xbtrg(:,2),Xbtrg(:,3));  
        th(th<0)=th(th<0)+2*pi; phi=pi/2-phi; 
              
        Vh_ngh = Vh(:,nbox); 
        
        % Compute K[X_ngh,Y_box] using spherical harmonics
        if kerd==1
            Ynear = Sh_Kernel_Eval_off(Vh_ngh,pMat,0,out,rho,phi,th,Nrtrgb);
        elseif kerd==3
            Ynear = Vsh_Kernel_Eval_off(Vh_ngh,pMat,0,out,rho,phi,th,Nrtrgb);
        end
        
        if size(Ynear,2)>1
           Ynear = sum(Ynear,2);  
        end
        
        if rot && kerd==3
            %Rotate back
            Ynear = [Ynear(1:3:end,:);Ynear(2:3:end,:);Ynear(3:3:end,:)];  
            Ynear = reshape(Ynear,Nnear,[])*MRot{nbox};       
            Ynear = reshape(reshape(Ynear,[],3).',[],1); 
        end

        if rdif1
            Ynear = rda(nbox)*Ynear; 
        end
        
        else
           Ynear = []; 
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add far-away interactions 
        parnear = params; 
        parnear.nor = params.nor(I_box,:);  
        
        %boolean index for far away points
        Itrgv   = repmat(Itrg.',kerd,1); 
        Itrgv = Itrgv(:); 
        if sum(neigh{nbox})>0
            faridxv = ~I_neigh & Itrgv; 
        else
            faridxv = Itrgv;
        end
        
        % Far interactions (dense via Kernel Eval, non sense via FMM)
        if dense
            if nortrg
               parnear.nor = Nrtrgv(faridxv,:);  
            end
            
            parnear.W2 = params.W2(I_box); 
            if isfield(parnear,'ci')
                parnear.ci = repmat((1:kerd)',sum(faridxv)/kerd,1);
            end
            if isfield(parnear,'cj')
            parnear.cj = params.cj(I_box);
            end
            
            if sum(neigh{nbox})>0
                Y(I_neigh,:) = Y(I_neigh,:) + Ynear; 
            end
            
            % if dense, we add direct kernel evaluation
            if sum(faridxv)>0
                Y(faridxv,:) =  Y(faridxv,:) + Kernel_Eval(Xtrgv(faridxv,:),Xv(I_box,:),parnear)*V(I_box,:); 
            end
        else
            if nortrg
               parnear.nor = Nrtrgv(I_neigh,:);  
            end
            parnear.W2 = params.W2(I_box); 
            if isfield(parnear,'ci')
                parnear.ci = repmat((1:kerd)',sum(I_neigh)/3,1);
            end
            if isfield(parnear,'cj')
            parnear.cj = params.cj(I_box);
            end
            
            % if not dense, we subtract neighbor kernel evaluation (to
            % cancel that contribution from FMM apply)
            Kngh = Kernel_Eval(Xtrgv(I_neigh,:),Xv(I_box,:),parnear); 
            Y(I_neigh,:) = Y(I_neigh,:) + Ynear - Kngh*V(I_box,:); 
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Add FMM / corrections 
    
    if ~dense
       % Setup and add FMM 
       if nortrg
          NrFMM = Nrtrg;  
       else
          NrFMM = params.nor(1:kerd:end,:); 
       end
       W = params.W2.'; 
       XFMM = Xv(1:kerd:end,:); 
       Y = Y + LOCAL_FMM_Eval(V,W,kerd,pot,Xtrg,XFMM,NrFMM);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(V,'Mat')
    % Build all-to-all interaction matrix densely
    if nortrg
       params.nor = Nrtrgv;  
    end
    
    % Start with Kernel Eval (correct for far-interactions)
    Y = Kernel_Eval(Xtrgv,Xv,params); 
    
    % Self and Near eval correction using spherical harmonics 
    for nbox=1:n3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Near interactions
        I_box = (1:Nb)+Nb*(nbox-1); 
        % Build neighbor index
        neigh{nbox} = reshape(neigh{nbox},1,[]); 
        
        if sum(neigh{nbox})>0
        I_neigh   = repmat(neigh{nbox},kerd,1); 
        I_neigh = I_neigh(:); 
        Nnear = sum(neigh{nbox}); 
        
         % near target points relative to box center
        Xbtrg = (1/rd(nbox))*(Xtrg(neigh{nbox},:) - repmat(C(nbox,:),Nnear,1)); 
        Nrtrgb = Nrtrg(neigh{nbox},:); 
        Nrtrgb = reshape(repmat(Nrtrgb,1,kerd)',3,[])';    
        
        % rotate back to reference
        if rot
           Xbtrg = Xbtrg*MRot{nbox}'; 
        end
        
        % spherical coordinates
        [th,phi,rho] = cart2sph(Xbtrg(:,1),Xbtrg(:,2),Xbtrg(:,3));  
        th(th<0)=th(th<0)+2*pi; phi=pi/2-phi; 
        
        if kerd==1 
            Ynear=Sh_Kernel_Eval_off([p 1],pMat,0,out,rho,phi,th,Nrtrgb);
        elseif kerd==3
            Ynear=Vsh_Kernel_Eval_off([p 1],pMat,0,out,rho,phi,th,Nrtrgb); 
        end
        
        Y(I_neigh,I_box) = rda(nbox)*Ynear; 
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    Y = @(V) VSh_MatVec_RB_trg(V,L,params);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y = LOCAL_FMM_Eval(Q,W,kerd,pot,Xtrg,Xsrc,Nr)

% source and target points variables
target = Xtrg.';     
ntarget = size(Xtrg,1);
if nargin<7
    Xsrc=Xtrg; 
end
source  = Xsrc.'; 
nsource = size(Xsrc,1);  

if strcmp(pot(1:2),'SL') || strcmp(pot(2:3),'SL')
    % SL and dSL/TSL
    sigma_sl = reshape(W.*Q,kerd,nsource); 
    ifsingle=1; 
    ifdouble=0; 
    sigma_dl=zeros(kerd,nsource); 
    sigma_dv=zeros(3,nsource); 
    ifpot=1;  ifpottarg=0;
    ifgradtarg=0;
else
    % DL and dDL/TDL
    sigma_dl = reshape(W.*Q,kerd,nsource); 
    ifsingle=0; 
    ifdouble=1; 
    sigma_sl=zeros(kerd,nsource); 
    sigma_dv=Nr.'; 
    ifpot=1;  ifpottarg=0;
    ifgradtarg=0;
end

if ~strcmp(pot(1:2),'SL') && ~strcmp(pot(1:2),'DL') 
    ifgrad=1;
else
    ifgrad=0; 
end

% precision for FMM, roughly 3*iprec digits of acc
iprec=2; 

if kerd==1
% Laplace particle FMM 
U=lfmm3dpart(iprec,nsource,source,ifsingle,sigma_sl,ifdouble,sigma_dl,...
    sigma_dv,ifpot,ifgrad,ntarget,target,ifpottarg,ifgradtarg);  
else
% Stokes particle FMM 
U=stfmm3dpart(iprec,nsource,source,ifsingle,sigma_sl,ifdouble,sigma_dl,...
    sigma_dv,ifpot,ifgrad,ntarget,target,ifpottarg,ifgradtarg);
end

% Evaluate depending on pot
switch pot
    case 'SL_L_3D'
        % Single layer potential at targets
        Y    = (1/4/pi)*U.pot.'; 
    case 'dSL_L_3D'
        % compute du/dNrtrg
        GSF = -(1/4/pi)*U.fld; % Gradient, size 3 x ntarget
        Y = sum(GSF.*Nr.'); Y=Y(:); 
    case 'DL_L_3D'
         % Double layer potential at targets
        Y    = (1/4/pi)*U.pot.'; 
    case 'SL_Stk_3D'
        % Single layer potential at targets
        Y    = (1/4/pi)*U.pot.'; 
        Y = real(reshape(Y.',[],1));
    case 'DL_Stk_3D'
        % Double layer potential at targets (check ct 1/4/pi)
        Y    = (1/4/pi)*U.pot.'; 
        Y = real(reshape(Y.',[],1));
    case 'TSL_Stk_3D'
        % Pressure
        SFpre = (1/4/pi)*U.pre; 
        % Gradient and Gradient transposed
        GSF   = (1/4/pi)*U.grad; 
        GTSF  = permute(GSF,[2 1 3]);  

        % Compute -pNr+Gu*Nr+Gut*Nr
        PNF = repmat(SFpre.',1,3).*Nr; 
        NrT = zeros(3,3,ntarget); NrT(:,1,:)=Nr.'; NrT(:,2,:)=Nr.'; NrT(:,3,:)=Nr.';
        GuN = reshape(sum(GSF.*NrT),[3 ntarget])+reshape(sum(GTSF.*NrT),[3 ntarget]); 

        Y  = -PNF.'+GuN; 
        Y = reshape(Y,[],1);   
    case 'TDL_Stk_3D'
        % Pressure
        SFpre = (1/4/pi)*U.pre; 
        % Gradient and Gradient transposed
        GSF   = (1/4/pi)*U.grad; 
        GTSF  = permute(GSF,[2 1 3]);  

        % Compute -pNr+Gu*Nr+Gut*Nr
        PNF = repmat(SFpre.',1,3).*Nr; 
        NrT = zeros(3,3,ntarget); NrT(:,1,:)=Nr.'; NrT(:,2,:)=Nr.'; NrT(:,3,:)=Nr.';
        GuN = reshape(sum(GSF.*NrT),[3 ntarget])+reshape(sum(GTSF.*NrT),[3 ntarget]); 

        Y  = -PNF.'+GuN; 
        Y = reshape(Y,[],1); 
    otherwise
        Y = zeros(ntarget,size(Q,2)); 
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function den = LOCAL_CenterDistance(X,C)

[Y_g1,  X_g1  ] = meshgrid(C(:,1), X(:,1));
[Y_g2,  X_g2  ] = meshgrid(C(:,2), X(:,2));
[Y_g3,  X_g3  ] = meshgrid(C(:,3), X(:,3));
d1 = (X_g1 - Y_g1); d2 = (X_g2 - Y_g2); d3 = (X_g3 - Y_g3); 
den = sqrt(d1.^2 + d2.^2 + d3.^2); 

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Itrg,neigh] = LOCAL_find_neigh(Xtrg,C,rd,mdst)

Ntrg = size(Xtrg,1); 
n3 = size(C,1); 
Prod = Ntrg*n3; 
big = Prod > 5*10^6; 
if ~big 
   distC = LOCAL_CenterDistance(Xtrg,C); 
end

neigh = cell(n3,1); 
inside = neigh; 
Itrg = true(Ntrg,1); 
eps=1e-6; 
display(mdst)

for i=1:n3
    if big
        distCXi = LOCAL_CenterDistance(Xtrg,C(i,:)); 
    else
        distCXi = distC(:,i);
    end
    
    inside{i} = distCXi<rd(i)-eps;
    Itrg = Itrg & ~inside{i}; 
    neigh{i} = distCXi<(mdst-1)*rd(i); 
end

for i=1:n3
    neigh{i} = neigh{i} & Itrg; 
end

end
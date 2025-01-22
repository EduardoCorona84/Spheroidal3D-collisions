function Y = VSh_MatVec_RB(V,L,params)
%{
This code evaluates the integral equation indicated by flag_pot in the params 
struct with density V: 

Y = S_{Gamma}[Q](x) = int_{Gamma} K(x,y)V(y)dSy    ;  x in Gamma

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

%}

persistent ShTg ShTr ShY Sc

p  = params.p;      % Spharm degree p 
n3 = params.n3;     % Number of objects
kerd = params.kerd; % Kernel dimension
dense=params.dense; % Dense vs FMM far int 
C = params.C;       % object centers
out = params.out;   % outside vs inside sphere
% Rotation of spheres
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
   rd=params.rd;
   rda = rd.^(strcmp(params.flag_pot(1:3),'SL_')); 
   rdif1=true; 
else
  rd=1; rda=1; rdif1=false; 
end

if isfield(params,'neigh')
    neigh=params.neigh; % Neighbor list (cell(n3,1))
else
   distC = LOCAL_CenterDistance(C);  
   neigh = cell(n3,1);  
   for i=1:n3
      neigh{i} = find(distC(i,:)<rd(i)*params.mdist); 
   end
   params.neigh=neigh;  
end

pot=params.flag_pot; 
nortrg=true; 

switch pot(1:3)
case 'SL_'
   pMat = 'SMat';  
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
sp = (p+1)^2; 
Nb = kerd*np; 
N = Nb*n3;  

X = params.Xp;
Xv = params.X; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Off-d + diag contributions

% If V is not empty, compute A*V. Otherwise, construct A densely. 
if ~isempty(V) && isnumeric(V)    
    V = reshape(V,Nb,[]);
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
    Y = zeros(N,1); 
    
    for nbox=1:n3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Near interactions
        I_box = (1:Nb)+Nb*(nbox-1); 
        num_ngh = length(neigh{nbox}); 
        % Build neighbor sphere index
        neigh{nbox} = reshape(neigh{nbox},1,[]); 
        I_neigh   = repmat((1:Nb)',1,num_ngh)+Nb*(repmat(neigh{nbox},Nb,1)-1);  
        I_neigh = I_neigh(:); 
        
        Xbv = Xv(I_box,:); %box pts
        nghidx = reshape(repmat(neigh{nbox},np,1),[],1); 
         % target points relative to source neighbor box centers
        rngh = rd(neigh{nbox}); rngha = rda(neigh{nbox}); 
        rmult = (1./repmat(reshape(repmat(rngh(:).',np,1),[],1),1,3)); 
        Xtrg = rmult.*(repmat(Xbv(1:kerd:end,:),num_ngh,1) - C(nghidx,:));  
        
        % rotate back to reference
        if rot
           Xtrg = Xtrg*MRot{nbox}'; 
        end
        
        % spherical coordinates
        [th,phi,rho] = cart2sph(Xtrg(:,1),Xtrg(:,2),Xtrg(:,3));  
        th(th<0)=th(th<0)+2*pi; phi=pi/2-phi; 
        
        %boolean index for far away points
        faridxv = true(N,1); 
        faridxv(I_neigh)=false; 
        Nrv = params.nor(I_box,:); 
        Nrtrg = repmat(Nrv,num_ngh,1);
        
        Vh_ngh = Vh(:,neigh{nbox}); 
        
        if rdif1
            Vh_ngh = repmat(rngha(:).',kerd*sp,1).*Vh_ngh;
        end
        
        if kerd==1
            if num_ngh>1
                slf = find(neigh{nbox}==nbox); 
                ind_off  = [1:np*(slf-1) (np*slf+1):np*num_ngh].';
                Vh_slf = Vh_ngh(:,slf); 
                Vh_ngh=[Vh_ngh(:,1:slf-1) Vh_ngh(:,slf+1:end)];  
                
                Ynear = Sh_Kernel_Eval_off(Vh_slf,pMat,0,out,1,[],[],[]) ...
                    + Sh_Kernel_Eval_off(Vh_ngh,pMat,0,out,rho(ind_off),phi(ind_off),th(ind_off),Nrtrg(ind_off,:));
            else
                Ynear = Sh_Kernel_Eval_off(Vh_ngh,pMat,0,out,1,[],[],[]);
            end
        elseif kerd==3
            if num_ngh>1
                slf = find(neigh{nbox}==nbox); 
                indv_off = [1:Nb*(slf-1) (Nb*slf+1):Nb*num_ngh].';
                ind_off  = [1:np*(slf-1) (np*slf+1):np*num_ngh].';
                Vh_slf = Vh_ngh(:,slf); 
                Vh_ngh=[Vh_ngh(:,1:slf-1) Vh_ngh(:,slf+1:end)];  
                
                Ynear = Vsh_Kernel_Eval_off(Vh_slf,pMat,0,out,1,[],[],[]) ...
                    + Vsh_Kernel_Eval_off(Vh_ngh,pMat,0,out,rho(ind_off),phi(ind_off),th(ind_off),Nrtrg(indv_off,:));
            else
                Ynear = Vsh_Kernel_Eval_off(Vh_ngh,pMat,0,out,1,[],[],[]);
            end
        end
        
        if size(Ynear,2)>1
           Ynear = sum(Ynear,2);  
        end 
        
        if rot && kerd==3
            %Rotate back
            Ynear = [Ynear(1:3:end,:);Ynear(2:3:end,:);Ynear(3:3:end,:)];  
            Ynear = reshape(Ynear,np,[])*MRot{nbox};       
            Ynear = reshape(reshape(Ynear,[],3).',[],1); 
        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add far-away interactions 
        parnear = params;
        parnear.nor = params.nor(I_box,:);
        
        if isfield(parnear,'ci')
        parnear.ci = params.ci(I_box); 
        end
        
        if dense
            if ~nortrg
                parnear.nor = params.nor(faridxv,:); 
            end
            parnear.W2 = params.W2(faridxv); 
            if isfield(parnear,'cj')
            parnear.cj = params.cj(faridxv);
            end
            
            % if dense, we add direct kernel evaluation
            if sum(faridxv)>0
                Y(I_box,:) = Ynear + Kernel_Eval(Xbv,Xv(faridxv,:),parnear)*V(faridxv,:); 
            else
                Y(I_box,:) = Ynear; 
            end
        else
            if ~nortrg
                parnear.nor = params.nor(~faridxv,:); 
            end
            parnear.W2 = params.W2(~faridxv); 
            if isfield(parnear,'cj')
            parnear.cj = params.cj(~faridxv);
            end
            
            % if not dense, we subtract neighbor kernel evaluation (to
            % cancel that contribution from FMM apply)
            Y(I_box,:) = Ynear - Kernel_Eval(Xbv,Xv(~faridxv,:),parnear)*V(~faridxv,:); 
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Add FMM / corrections 
    
    if ~dense
       % Setup and add FMM 
       Nr = params.nor(1:kerd:end,:); 
       W = params.W2.'; 
       Y = Y + LOCAL_FMM_Eval(V,W,kerd,pot,X,X,Nr); 
    end
    
    if ~strcmp(pot(1:3),'SL_') 
    if isfield(params,'a')
       if strcmp(pot(2:3),'SL') 
       if isfield(params,'eta') && strcmp(pot(1:3),'dSL') 
           Y = params.eta*Y + (params.a+0.5)*V; 
       else
           Y = Y + (params.a+0.5)*V; 
       end
       else
           Y = Y + (params.a-0.5)*V; 
       end
    end
    end
    
    if ~isempty(L)
    % Add block-diag nullspace correction
    Y = Y + L*V; 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif strcmp(V,'Mat')
    % Build all-to-all interaction matrix densely
    
    % Start with Kernel Eval (correct for far-interactions)
    Y = Kernel_Eval(Xv,Xv,params); 
    
    % Self and Near eval correction using spherical harmonics 
    for nbox=1:n3
        % Box index 
        I_box = (1:Nb)+Nb*(nbox-1); 
        Xbv = Xv(I_box,:); %box pts
        
        % Build neighbor sphere index and points 
        num_ngh = length(neigh{nbox}); 
        neigh{nbox} = reshape(neigh{nbox},1,[]); 
        I_neigh   = repmat((1:Nb)',1,num_ngh)+Nb*(repmat(neigh{nbox},Nb,1)-1);  
        I_neigh = I_neigh(:);
        
        nghidx = reshape(repmat(neigh{nbox},np,1),[],1); 
         % target points relative to source neighbor box centers
        rngh = rd(neigh{nbox}); rngha = rda(neigh{nbox}); 
        rmult = (1./repmat(reshape(repmat(rngh(:).',np,1),[],1),1,3)); 
        Xtrg = rmult.*(repmat(Xbv(1:kerd:end,:),num_ngh,1) - C(nghidx,:));  
        
        % rotate back to reference
        if rot
           Xtrg = Xtrg*MRot{nbox}'; 
        end
        
        % spherical coordinates
        [th,phi,rho] = cart2sph(Xtrg(:,1),Xtrg(:,2),Xtrg(:,3));  
        th(th<0)=th(th<0)+2*pi; phi=pi/2-phi; 
        
        %boolean index for far away points
        Nrv = params.nor(I_box,:); 
        Nrtrg = repmat(Nrv,num_ngh,1);
        
        Rdiag = repmat(rngha(:).',kerd*sp,1); 
        
        if kerd==1 
            Ynear=Sh_Kernel_Eval_off([p num_ngh Rdiag(:).'],pMat,0,out,rho,phi,th,Nrtrg);
        elseif kerd==3
            Ynear=Vsh_Kernel_Eval_off([p num_ngh Rdiag(:).'],pMat,0,out,rho,phi,th,Nrtrg); 
        end
        
        Y(I_box,I_neigh) = Ynear; 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Add diag and block-diagonal corrections        
        if ~strcmp(pot(1:3),'SL_') 
            if isfield(params,'a')
                if strcmp(pot(2:3),'SL') 
                    if isfield(params,'eta') && strcmp(pot(1:3),'dSL') 
                        Y = params.eta*Y + (params.a+0.5)*eye(N); 
                    else
                        Y = Y + (params.a+0.5)*eye(N); 
                    end
                else
                    Y = Y + (params.a-0.5)*eye(N); 
                end
            end
        end
    
        if ~isempty(L)
            % Add block-diag nullspace correction
            Y = Y + L; 
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    Y = @(V) VSh_MatVec_RB(V,L,params);
end

end

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

function den = LOCAL_CenterDistance(C)

[Y_g1,  X_g1  ] = meshgrid(C(:,1), C(:,1));
[Y_g2,  X_g2  ] = meshgrid(C(:,2), C(:,2));
[Y_g3,  X_g3  ] = meshgrid(C(:,3), C(:,3));
d1 = (X_g1 - Y_g1); d2 = (X_g2 - Y_g2); d3 = (X_g3 - Y_g3); 
den = sqrt(d1.^2 + d2.^2 + d3.^2); 

end

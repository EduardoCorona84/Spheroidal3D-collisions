function Y = VSh_MatVec_RB_near(V,L,params)
%{
This code evaluates the integral equation indicated by flag_pot in the params 
struct with density V: 

Y = S_{Gamma}[Q](x) = int_{Gamma} K(x,y)V(y)dSy    ;  x in Gamma

Where Gamma is the union of boundaries of spheres with radii r(i) and centers
C(:,i), by discretizing Q using spherical harmonics of degree p. 

This code calculates ONLY interactions between neighboring spheres. 

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
dense=params.dense; % Dense vs FMM off diagonal 
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

if isfield(params,'doAna')
   doAna = params.doAna;  
else
   doAna = 1;  
end

% Radius of spheres
if isfield(params,'rd')
   if length(params.rd)==1
       params.rd = repmat(params.rd,1,n3); 
   end
   
   rd=params.rd;
   rda = rd.^(strcmp(params.flag_pot(1:3),'SL_')); 
   rdif1=true; 
else
  rd=1; rda=1; rdif1=false; 
end

if isfield(params,'neigh')
   neigh=params.neigh; % Neighbor list (cell(n3,1)) 
elseif isfield(params,'nghmat')
   neigh = cell(n3,1);  
   for i=1:n3
      neigh{i} = find(params.nghmat(:,i)==1); 
   end
   params.neigh=neigh;  
else
   distC = LOCAL_CenterDistance(C);
   nghmat=false(n3,n3); neigh = cell(n3,1);
   for j=1:n3
      nghmat(:,j) = distC(:,j)<rd(j)*params.mdist;  
      neigh{j} = find(nghmat(:,j)==1);
   end
   params.neigh=neigh;  
end

pot=params.flag_pot; 
nortrg=true; 

switch pot(1:3)
case 'SL_'
   pMat = 'SMat';  
case 'SDL'
   pMat = 'SDMat'; 
   nortrg = false; 
case 'dSL'
   pMat = 'SpMat'; 
case 'TSL' 
   pMat='TSMat'; 
   %Coefficient matrices (add option to load only once) 
    if isempty(ShTg) || size(ShTg{1},1)/3 ~= (p+1)^2
        fprintf('\n Computing Coeff matrices for T \n')
        [ShTg,ShTr,ShY] = get_Surfgrad_coeffs(1.5*p,p,1);
        %[ShTg,ShTr,ShY] = Surfgrad_coeffs(ceil(1.5*p),p,1,'VW','VW',0); 
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
        [ShTg,ShTr,ShY] = get_Surfgrad_coeffs(1.5*p,p,1);
        %[ShTg,ShTr,ShY] = Surfgrad_coeffs(ceil(1.5*p),p,1,'VW','VW',0); 
    end
    nortrg=false; 
end

np = 2*p*(p+1); 
Nb = kerd*np; 
N = Nb*n3;  
X = params.Xp;
Xv = params.X; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Off-d + diag contributions

% If V is not empty, compute A*V. Otherwise, construct A densely. 
if ~isempty(V) && isnumeric(V)   
    V = reshape(V,Nb,[]);
    
    if doAna
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
    else
       Vh = params.Vh;  
    end
    
    V = reshape(V,N,[]); 
    Y = zeros(N,size(V,2));
    
    for nbox=1:n3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Near interactions
        I_box = (1:Nb)+Nb*(nbox-1); 
        num_ngh = length(neigh{nbox}); 
        % Build neighbor sphere index
        neigh{nbox} = reshape(neigh{nbox},1,[]);         
        I_nghv   = repmat((1:Nb)',1,num_ngh)+Nb*(repmat(neigh{nbox},Nb,1)-1);  
        I_nghv = I_nghv(:); 
        
        Xbv = Xv(I_nghv,:); %box pts
        %nghidx = reshape(repmat(neigh{nbox},np,1),[],1); 
         % target points relative to source neighbor box centers
        Xtrg = (1/rd(nbox))*(Xbv(1:kerd:end,:) - repmat(C(nbox,:),num_ngh*np,1)); 
        Nrtrg = params.nor(I_nghv,:); 
        
        % rotate back to reference
        if rot
           Xtrg = Xtrg*MRot{nbox}'; 
        end
        
        % spherical coordinates
        [th,phi,rho] = cart2sph(Xtrg(:,1),Xtrg(:,2),Xtrg(:,3));  
        th(th<0)=th(th<0)+2*pi; phi=pi/2-phi; 
        
        %boolean index for far away points
        faridxv = true(N,1); 
        faridxv(I_nghv)=false;         
        Vh_ngh = Vh(:,nbox); 
        
        if kerd==1
            if num_ngh>1
                slf = find(neigh{nbox}==nbox); 
                ind_off  = [1:np*(slf-1) (np*slf+1):np*num_ngh].';
                
                Ynear = zeros(np*num_ngh,1); 
                Ynear(np*(slf-1)+1:np*slf) = Sh_Kernel_Eval_off(Vh_ngh,pMat,0,out,1,[],[],[]); 
                Ynear(ind_off) = Sh_Kernel_Eval_off(Vh_ngh,pMat,0,out,rho(ind_off),phi(ind_off),th(ind_off),Nrtrg(ind_off,:));
            else
                Ynear = Sh_Kernel_Eval_off(Vh_ngh,pMat,0,out,1,[],[],[]);
            end
        elseif kerd==3
            if num_ngh>1
                slf = find(neigh{nbox}==nbox); 
                indv_off = [1:Nb*(slf-1) (Nb*slf+1):Nb*num_ngh].';
                ind_off  = [1:np*(slf-1) (np*slf+1):np*num_ngh].';
                
                Ynear = zeros(Nb*num_ngh,1);
                Ynear(Nb*(slf-1)+1:Nb*slf) = Vsh_Kernel_Eval_off(Vh_ngh,pMat,0,out,1,[],[],[]); 
                Ynear(indv_off) = Vsh_Kernel_Eval_off(Vh_ngh,pMat,0,out,rho(ind_off),phi(ind_off),th(ind_off),Nrtrg(indv_off,:));
            else
                Ynear = Vsh_Kernel_Eval_off(Vh_ngh,pMat,0,out,1,[],[],[]);
                %Ynear = Vsh_Kernel_Eval_off(Vh_ngh,pMat,0,out,rho,phi,th,Nrtrg);
            end
        end
        
        if size(Ynear,2)>1
           Ynear = sum(Ynear,2);  
        end
        
        if rot && kerd==3
            %Rotate back
            Ynear = [Ynear(1:3:end,:);Ynear(2:3:end,:);Ynear(3:3:end,:)];  
            Ynear = reshape(Ynear,[],3)*MRot{nbox};       
            Ynear = reshape(reshape(Ynear,[],3).',[],1); 
        end

        if rdif1
            Ynear = rda(nbox)*Ynear; 
        end
        
        Y(I_nghv,:) = Y(I_nghv,:) + Ynear;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Add diagonal corrections 
    
    sgo = out - ~out; 
    ct = 0.5*sgo; 

    if ~strcmp(pot(1:3),'SL_') 
    if isfield(params,'a')
       if strcmp(pot(2:3),'SL') 
       if isfield(params,'eta') && strcmp(pot(1:3),'dSL') 
           Y = params.eta*Y + (params.a+ct)*V; 
       else
           Y = Y + (params.a+ct)*V; 
       end
       else
           Y = Y + (params.a-ct)*V; 
       end
    end
    end
    
    if ~isempty(L)
    % Add block-diag nullspace correction
    Y = Y + L*V; 
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif strcmp(V,'Mat')
    % Build all-to-all neighbor interaction matrix as dense / spase matrix
    
    % Start with Kernel Eval (correct for far-interactions)
    if ~dense
       nnzKp = nnz(params.nghmat)*Nb^2; 
       IK = zeros(nnzKp,1); JK = IK; VK = IK; 
       nnzK=0; 
    else
       Y = zeros(N,N);  
    end
    
    % Self and Near eval correction using spherical harmonics 
    for nbox=1:n3
        % Near interactions
        I_box = (1:Nb)+Nb*(nbox-1); 
        num_ngh = length(neigh{nbox}); 
        % Build neighbor sphere index
        neigh{nbox} = reshape(neigh{nbox},1,[]);         
        I_nghv   = repmat((1:Nb)',1,num_ngh)+Nb*(repmat(neigh{nbox},Nb,1)-1);  
        I_nghv = I_nghv(:); 
        
        Xbv = Xv(I_nghv,:); %box pts
        %nghidx = reshape(repmat(neigh{nbox},np,1),[],1); 
         % target points relative to source neighbor box centers
        Xtrg = (1/rd(nbox))*(Xbv(1:kerd:end,:) - repmat(C(nbox,:),num_ngh*np,1)); 
        Nrtrg = params.nor(I_nghv,:); 
        
        % rotate back to reference
        if rot
           Xtrg = real(Xtrg*MRot{nbox}'); 
        end
      
        Xtrg = real(Xtrg); 
        
        % spherical coordinates
        [th,phi,rho] = cart2sph(Xtrg(:,1),Xtrg(:,2),Xtrg(:,3));  
        th(th<0)=th(th<0)+2*pi; phi=pi/2-phi;
        
        if kerd==1 
            Ynear=Sh_Kernel_Eval_off([p 1],pMat,0,out,rho,phi,th,Nrtrg);
        elseif kerd==3
            Ynear=Vsh_Kernel_Eval_off([p 1],pMat,0,out,rho,phi,th,Nrtrg); 
        end
        
        if rot
            %Rotate back
            Ynear = permute(reshape(Ynear,[3 np*num_ngh Nb]),[2 3 1]); 
            Ynear = reshape(Ynear,[],3)*MRot{nbox};   
            Ynear = permute(reshape(Ynear,[np*num_ngh Nb 3]),[3 1 2]); 
            Ynear = reshape(Ynear,[],Nb); 
        end
        
        if dense
            Y(I_nghv,I_box) = rda(nbox)*Ynear; 
        else
            [IIb,IIn] = meshgrid(I_box,I_nghv); 
            nY = numel(Ynear); 
            IK(nnzK+1:nnzK+nY)=IIn; 
            JK(nnzK+1:nnzK+nY)=IIb;
            VK(nnzK+1:nnzK+nY)=rda(nbox)*Ynear(:);
            nnzK = nnzK+nY; 
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % sparse mat build
    if ~dense
       Y = sparse(IK,JK,VK,N,N);  
       IdN = speye(N); 
    else
       IdN = eye(N); 
    end
    
    % Add diag and block-diagonal corrections
    sgo = out - ~out; 
    ct = 0.5*sgo;
        
    if ~strcmp(pot(1:3),'SL_') 
        if isfield(params,'a')
            if strcmp(pot(2:3),'SL') 
                if isfield(params,'eta') && strcmp(pot(1:3),'dSL') 
                    fprintf('\n eta = %1.2f',params.eta); 
                    Y = params.eta*Y + (params.a+ct)*IdN; 
                else
                    Y = Y + (params.a+ct)*IdN;
                end
            else
                Y = Y + (params.a-ct)*IdN;
            end
        end
    end
    
    if ~isempty(L)
        % Add block-diag nullspace correction
        Y = Y + L; 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    Y = @(V) VSh_MatVec_RB_near(V,L,params);
end

end

function den = LOCAL_CenterDistance(C)

[Y_g1,  X_g1  ] = meshgrid(C(:,1), C(:,1));
[Y_g2,  X_g2  ] = meshgrid(C(:,2), C(:,2));
[Y_g3,  X_g3  ] = meshgrid(C(:,3), C(:,3));
d1 = (X_g1 - Y_g1); d2 = (X_g2 - Y_g2); d3 = (X_g3 - Y_g3); 
den = sqrt(d1.^2 + d2.^2 + d3.^2); 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ShTg,ShTr,ShY] = get_Surfgrad_coeffs(px,p,sprs)

if sprs
  printMsg('* Reading the Traction kernel coefficient matrices for p=%d: ',p);
  [TrSz,rFlag] = readData('TractionCoeffsSpSize',p,1);

  if(~rFlag)
    printMsg('Failed\n');
    printMsg('* Generating the Traction kernel coefficient matrices for p=%d\n',p);

    [ShTg,ShTr,ShY] = Surfgrad_coeffs(px,p,1,'VW','VW',0);
    
    TrData = zeros(0,5); 
    for j=1:3
       [ig,jg,vg] = find(ShTg{j}); ng = [length(vg) ; zeros(length(vg)-1,1)]; 
       [ir,jr,vr] = find(ShTr{j}); nr = [length(vr) ; zeros(length(vr)-1,1)];
       [iy,jy,vy] = find(ShY{j}); ny = [length(vy) ; zeros(length(vy)-1,1)];
       TrData = [TrData ; [ig ; ir ; iy] [jg ; jr ; jy] [real(vg) ; real(vr) ; real(vy)] [imag(vg) ; imag(vr) ; imag(vy)] [ng ; nr ; ny]];  
    end
    
    writeData('TractionCoeffsSp',p, TrData);
    writeData('TractionCoeffsSpSize',p, size(TrData,1));
    printMsg('* Stored generated Traction coefficient matrices for p=%d\n', p);
  else
    sp = (p+1)^2; 
    [TrData,rFlag] = readData('TractionCoeffsSp',p,[TrSz 5]);  
    ShTg = cell(3,1); ShTr = ShTg; ShY = ShTg;  
    lv = [find(TrData(:,5)>0) ; TrSz];
    sp3 = 3*sp; 
    
    for j=1:3
        ind = lv(1+3*(j-1)):(lv(2+3*(j-1))-1); 
        ShTg{j} = sparse(TrData(ind,1),TrData(ind,2),complex(TrData(ind,3),TrData(ind,4)),sp3,sp3); 
        ind = lv(2+3*(j-1)):(lv(3+3*(j-1))-1); 
        ShTr{j} = sparse(TrData(ind,1),TrData(ind,2),complex(TrData(ind,3),TrData(ind,4)),sp3,sp3);
        ind = lv(3+3*(j-1)):(lv(4+3*(j-1))-1); 
        ShY{j} = sparse(TrData(ind,1),TrData(ind,2),complex(TrData(ind,3),TrData(ind,4)),sp3,sp3);
    end
    
    printMsg('Successful\n');
  end
  
else
    printMsg('* Reading the Traction kernel coefficient matrices for p=%d: ',p);
    sp = (p+1)^2; 
  [TrData,rFlag] = readData('TractionCoeffsDn',p,[9*sp 18*sp]);

  if(~rFlag)
    printMsg('Failed\n');
    printMsg('* Generating the Traction kernel coefficient matrices for p=%d\n',p);

    [ShTg,ShTr,ShY] = Surfgrad_coeffs(px,p,1,'VW','VW',0);
    
    TrData = [full(real(cell2mat(ShTg))) full(imag(cell2mat(ShTg))) ...
        full(real(cell2mat(ShTr))) full(imag(cell2mat(ShTr))) ...
        full(real(cell2mat(ShY))) full(imag(cell2mat(ShY)))]; 
    
    writeData('TractionCoeffsDn',p, TrData);
    printMsg('* Stored generated Traction coefficient matrices for p=%d\n', p);
  else 
    ShTg = cell(3,1); ShTr = ShTg; ShY = ShTg;  
    
    for j=1:3
        ind = (1+3*sp*(j-1)):3*sp*j;  
        ShTg{j} = complex(TrData(ind,1:3*sp),TrData(ind,3*sp+1:6*sp));   
        ShTr{j} = complex(TrData(ind,6*sp+1:9*sp),TrData(ind,9*sp+1:12*sp));
        ShY{j} = complex(TrData(ind,12*sp+1:15*sp),TrData(ind,15*sp+1:18*sp));
    end
    
    printMsg('Successful\n');
  end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

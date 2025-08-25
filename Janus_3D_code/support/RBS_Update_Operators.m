function [Kernels,Nullsp,Fparams,timings] = RBS_Update_Operators(Xt,Ct,Mt,nrmW,Kernels,Fparams,timings,i)
%global denseMV type typeMV timings rd mdist out; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sdim=3; ldim=1; 
parslv = Fparams.parslv; 
prec=parslv.prec; prev=parslv.prev; prtype=parslv.prtype; 
typeMV = Fparams.typeMV; denseMV = Fparams.denseMV; 

tic; 

np = Fparams.parbd.np; n3 = size(Ct,1); p = Fparams.parbd.p; 
dodense=Fparams.denseforce;
Nb = Fparams.parbd.Nb; 
mdist = Fparams.parbd.mdist; out = Fparams.parbd.out; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fparams.parbd update

Shape = Fparams.parbd.Shape; 

if strcmp(typeMV,'Vsh') || strcmp(Shape,'')
    Sc = Fparams.parbd.Sc;  
else 
Sc = cell(p,n3);
for k=1:n3
    Sc{p,k} = SurfaceSph(vec3d(Xt((1:np)+np*(k-1),:)));
end    
end

Fparams.parbd = RBS_set_params(p,Ct,Fparams.parbd.rd,...
        'TSL_Stk_3D',Shape,Sc,sdim,Fparams.denseMV,Fparams.parbd.doAna,mdist,Fparams.parbd.eps,out);     
if(isfield(Fparams,'lambda'))
    Fparams.parmod = RBS_set_params(p,Ct,Fparams.parbd.rd,...
        'SL_LMOD_3D',Shape,Sc,ldim,Fparams.denseMV,Fparams.parbd.doAna,mdist,Fparams.parbd.eps,out);
    Fparams.parmod.lambda=Fparams.lambda;
end

if isfield(Fparams, 'lofi')
    lofi_p = Fparams.lofi.p;
    lofi_Sc = Fparams.lofi.Sc;
    lofi_rd = Fparams.lofi.rd;
    Fparams.lofi = RBS_set_params(lofi_p,Ct,lofi_rd,...
        'TSL_Stk_3D',Shape,lofi_Sc,sdim,Fparams.denseMV,...
        Fparams.lofi.doAna,mdist,Fparams.lofi.eps,out); 
end

Xt = Fparams.parbd.Xrp;

fprintf('\n Time for surface update: %e',toc) 
if i>0
timings.operator.surf(i) = toc;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Kernels / MatVec info

[Ck,Bk,Dk,Lk] = Build_AuxMats(Fparams.parbd.Wg,Xt,[],np,n3); 
Nullsp.C = Ck; Nullsp.B = Bk; Nullsp.D = Dk; Nullsp.L = Lk;
if isfield(Fparams, 'lofi')
    lofi_Xt = Fparams.lofi.Xrp;
    [lofi_Ck,lofi_Bk,lofi_Dk,lofi_Lk] = Build_AuxMats(Fparams.lofi.Wg,lofi_Xt,[],Fparams.lofi.np,n3); 
    Nullsp.lofi_C = lofi_Ck; Nullsp.lofi_B = lofi_Bk; Nullsp.lofi_D = lofi_Dk; Nullsp.lofi_L = lofi_Lk;
end

% Initial block diag build
if i==0
   % Single layer traction
   [~,Kernels.TSSDd] = RBS_Build_DMV('dense',Fparams,0.5,'TSL_Stk_3D'); 
   Kernels.TSSD0 = Kernels.TSSDd;
   Kernels.ITSSD0 = inv(Kernels.TSSDd+Lk(1:Nb,1:Nb));
   % Single layer kernel
   [~,Kernels.SSDd] = RBS_Build_DMV('dense',Fparams,0,'SL_Stk_3D');
   Kernels.SSD0 = Kernels.SSDd;
      
   if strcmp(prtype,'TT')
      Kernels.ITSSDd = Kernels.TSSDd+Lk(1:Nb,1:Nb);
      Kernels.ITSSD0 = Kernels.ITSSDd;
   elseif strcmp(prtype,'bkdiag')   
      Kernels.ITSSDd = Kernels.ITSSD0; 
   end
   
   if strcmp(Fparams.type,'MHD')
       % Laplace kernels

       % Single layer
       [~,Kernels.SLD0] = RBS_Build_DMV('dense',Fparams,0,'SL_L_3D');

       % KL = (1/2 I + eta S')
       [~,Kernels.KLD0] = RBS_Build_DMV('dense',Fparams,0.5,'dSL_L_3D');
   end
end
     
if ~strcmp(typeMV,'Vsh') && n3>1 && i==1
    Kernels.TSSDd = cell(n3,1); Kernels.SSDd = cell(n3,1); 
    if denseMV==0 
       Kernels.ITSSDd = cell(n3,1);  
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diagonal block update (rotation) for typeMV=Rbs 
tic; 
if ~strcmp(typeMV,'Vsh') && i>0
for k=1:n3
   if n3>1
   if nrmW(k)>1e-10
       Kernels.TSSDd{k} = Rotate_Operator(Kernels.TSSD0,Mt{k},np); 
       Kernels.SSDd{k} = Rotate_Operator(Kernels.SSD0,Mt{k},np);
        
       if Fparams.denseMV==0
          Kernels.ITSSDd{k} = Rotate_Operator(Kernels.ITSSD0,Mt{k},np);  
       end
   elseif i==1
       Kernels.TSSDd{k} = Kernels.TSSD0; 
       Kernels.SSDd{k} = Kernels.SSD0;
        
       if Fparams.denseMV==0
           Kernels.ITSSDd{k} = Kernels.ITSSD0;
       end 
   end
   else
       if nrmW(k)>1e-10 
           Kernels.TSSDd = Rotate_Operator(Kernels.TSSD0,Mt{k},np); 
           Kernels.SSDd = Rotate_Operator(Kernels.SSD0,Mt{k},np);   
       end
   end
end
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set preconditioner
if n3>1 && Fparams.denseMV==0
    if strcmp(prtype,'TT')        
        [parslv.prec,parslv.prev] = RBS_setprec(Kernels.ITSSDd,n3,prtype,p,Fparams.parbd,prev); 
    elseif strcmp(prtype,'bkdiag')
        parslv.prec = RBS_setprec(Kernels.ITSSDd,n3,prtype);   
    end
end

Fparams.parslv = parslv; 
Fparams.parslv.Pr = @(x) reshape(VshProj(reshape(x,Nb,[]),'VW'),n3*Nb,[]); 

fprintf('\n Time for diagonal block and prec rotation: %e',toc) 
if i>0
timings.operator.diag(i) = toc; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic; 


%V,L,params,kerd,flag_pot,DMV
% Stokes kernels 
Kernels.TD = RBS_MatVec([],Lk,typeMV,Fparams.parbd,sdim,0.5,'TSL_Stk_3D',Kernels.TSSDd);
Kernels.SD = RBS_MatVec([],[],typeMV,Fparams.parbd,sdim,0,'SL_Stk_3D',Kernels.SSDd); 
if isfield(Fparams, 'lofi')
    if ~strcmpi(typeMV, 'vsh')
        warning('Lofi implementation expects no use of the DMV possible failure');
    end
    Kernels.lofi_TD = RBS_MatVec([],lofi_Lk,typeMV,Fparams.lofi,sdim,0.5,'TSL_Stk_3D',[]);
    Kernels.lofi_SD = RBS_MatVec([],[],typeMV,Fparams.lofi,sdim,0,'SL_Stk_3D',[]); 
end
Fparams.parmod.dense = 1;
if strcmp(Fparams.type,'MHD')
    % Laplace kernels
    Kernels.SLD = RBS_MatVec([],[],typeMV,Fparams.parbd,ldim,0,'SL_L_3D',Kernels.SLD0); 
    Kernels.KLD = RBS_MatVec([],[],typeMV,Fparams.parbd,ldim,0.5,'dSL_L_3D',Kernels.KLD0);
elseif strcmp(Fparams.type,'JanusEHD')
    Kernels.SLD0=[]; Kernels.DLD0 = []; 
    %compute kernels
   
    if (dodense == 1)
    
    Kernels.SLMODD = RBS_MatVec([],[],typeMV,Fparams.parmod,ldim,0,'SL_LMOD_3D',Kernels.SLD0);
    Kernels.DLMODD = RBS_MatVec([],[],typeMV,Fparams.parmod,ldim,0,'DL_LMOD_3D',Kernels.DLD0);
    Kernels.dSLMODD = RBS_MatVec([],[],typeMV,Fparams.parmod,ldim,0,'dSL_LMOD_3D',Kernels.DLD0);
    Kernels.dDLMODD = RBS_MatVec([],[],typeMV,Fparams.parmod,ldim,0,'dDL_LMOD_3D',Kernels.DLD0);
    
    Kernels.SLD = RBS_MatVec([],[],typeMV,Fparams.parmod,ldim,0,'SL_L_3D',Kernels.SLD0);
    Kernels.DLD = RBS_MatVec([],[],typeMV,Fparams.parmod,ldim,0,'DL_L_3D',Kernels.DLD0);
    Kernels.dSLD = RBS_MatVec([],[],typeMV,Fparams.parmod,ldim,0,'dSL_L_3D',Kernels.SLD0);
    Kernels.dDLD = RBS_MatVec([],[],typeMV,Fparams.parmod,ldim,0,'dDL_L_3D',Kernels.DLD0);
    
    else
    Kernels.SLMODD = @(V) RBS_MatVec(V,[],typeMV,Fparams.parmod,ldim,0,'SL_LMOD_3D',Kernels.SLD0);
    Kernels.DLMODD= @(V) RBS_MatVec(V,[],typeMV,Fparams.parmod,ldim,0,'DL_LMOD_3D',Kernels.DLD0);
    Kernels.dSLMODD = @(V) RBS_MatVec(V,[],typeMV,Fparams.parmod,ldim,0,'dSL_LMOD_3D',Kernels.DLD0);
    Kernels.dDLMODD = @(V) RBS_MatVec(V,[],typeMV,Fparams.parmod,ldim,0,'dDL_LMOD_3D',Kernels.DLD0);
    
    Kernels.SLD = @(V) RBS_MatVec(V,[],typeMV,Fparams.parmod,ldim,0,'SL_L_3D',Kernels.SLD0);
    Kernels.DLD = @(V) RBS_MatVec(V,[],typeMV,Fparams.parmod,ldim,0,'DL_L_3D',Kernels.DLD0);
    Kernels.dSLD = @(V) RBS_MatVec(V,[],typeMV,Fparams.parmod,ldim,0,'dSL_L_3D',Kernels.SLD0);
    Kernels.dDLD = @(V) RBS_MatVec(V,[],typeMV,Fparams.parmod,ldim,0,'dDL_L_3D',Kernels.DLD0);
    
    
    end
    
    
    
    
    Epsilon=zeros(n3*np,1);
    ind=false(np,n3);
    point_charges=zeros(size(Fparams.init_charges));
    
    %compute Epsilon
    init_dir=Fparams.init_dir;
    for i=1:n3
        Xsphere=Fparams.parmod.Xrp(np*(i-1)+1:np*i,:);
        Xbackrot=Xsphere*Mt{i};
        ind(:,i)=Xbackrot*init_dir(i,:)'>0;        
        point_charges(:,:,i)=Fparams.init_charges(:,:,i)*Mt{i}';
    end
        Fparams.point_charges=point_charges;
        ind=reshape(ind,[],1);
        Epsilon(ind)=Fparams.epsilon_north;
        Epsilon(~ind)=Fparams.epsilon_south;
        Fparams.Epsilon=Epsilon;
    %update point charges
elseif strcmp(Fparams.type,'JanusAmp')
    
     Kernels.SLD0=[]; Kernels.DLD0 = []; 
    if (dodense ==1)
    
    %compute kernels
    Kernels.SLMODD = RBS_MatVec([],[],typeMV,Fparams.parmod,ldim,0,'SL_LMOD_3D',Kernels.SLD0);
    Kernels.DLMODD = RBS_MatVec([],[],typeMV,Fparams.parmod,ldim,0.5,'DL_LMOD_3D',Kernels.DLD0);
    Kernels.dSLMODD = RBS_MatVec([],[],typeMV,Fparams.parmod,ldim,0,'dSL_LMOD_3D',Kernels.DLD0);
    Kernels.dDLMODD = RBS_MatVec([],[],typeMV,Fparams.parmod,ldim,0,'dDL_LMOD_3D',Kernels.DLD0);
    
    else
    Kernels.SLMODD = @(V) RBS_MatVec(V,[],typeMV,Fparams.parmod,ldim,0,'SL_LMOD_3D',Kernels.SLD0);
    Kernels.DLMODD = @(V) RBS_MatVec(V,[],typeMV,Fparams.parmod,ldim,0.5,'DL_LMOD_3D',Kernels.DLD0);
    Kernels.dSLMODD =@(V) RBS_MatVec(V,[],typeMV,Fparams.parmod,ldim,0,'dSL_LMOD_3D',Kernels.DLD0);
    Kernels.dDLMODD =@(V) RBS_MatVec([],[],typeMV,Fparams.parmod,ldim,0,'dDL_LMOD_3D',Kernels.DLD0)*V;
    end
    
    
    %compute hydrophilic label at each point
    %f: Sphere-->[0,1] vectorized function handle 
    flabel=Fparams.boundary_label;
    surfacelabel=zeros(np,n3);
    init_dir=Fparams.init_dir;
    for i=1:n3
       Xsphere=Fparams.parmod.Xrp(np*(i-1)+1:np*i,:);
       %Center=Fparams.parmod.C(i,:);
       %Xsphere=Xsphere-repmat(Center,np,1);
       Xbackrot=Xsphere*Mt{i}; 
       surfacelabel(:,i)=flabel(Xbackrot,init_dir(i,:));
       %fprintf('\n max surf label = %e, min surf label = %e ',max(surfacelabel(:,i)),min(surfacelabel(:,i)));
       %plot(surfacelabel(:,i)); 
    end
    %pause; 
    Fparams.SurfaceLabel=reshape(surfacelabel,[],1);


end

fprintf('\n Time for kernel eval update: %e',toc) 
if i>0
timings.operator.offd(i) = toc;
end

end

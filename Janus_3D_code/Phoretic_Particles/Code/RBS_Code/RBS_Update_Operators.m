function [Kernels,Nullsp,Fparams,timings] = RBS_Update_Operators(Xt,Ct,Mt,nrmW,Kernels,Fparams,timings,i)
%global denseMV type typeMV timings rd mdist out; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sdim=3; ldim=1; 
parslv = Fparams.parslv; 
prec=parslv.prec; prev=parslv.prev; prtype=parslv.prtype; 
typeMV = Fparams.typeMV; denseMV = Fparams.denseMV; 

tic; 

np = Fparams.parbd.np; n3 = size(Ct,1); p = Fparams.parbd.p;
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
    
if strcmp(Fparams.type,'MHD')
    Fparams.parbdL = RBS_set_params(p,Ct,Fparams.parbd.rd,...
        'SL_L_3D',Shape,Sc,ldim,Fparams.denseMV,Fparams.parbd.doAna,mdist,Fparams.parbd.eps,out);
    Fparams.parbdL.eta = Fparams.eta; 
end

Xt = Fparams.parbd.Xrp;

fprintf('\n Time for surface update: %e',toc) 
if i>0
timings.operator.surf(i) = toc;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Kernels / MatVec info

[Ck,Bk,Dk,Ak,Lk] = Build_AuxMats2(Fparams.parbd.Wg,Xt,[],np,n3); 
Nullsp.C = Ck; Nullsp.B = Bk; Nullsp.D = Dk; Nullsp.A = Ak;Nullsp.L = Lk;

% Initial block diag build
if i==0
   Sc0 = cell(p,1);
   Sc0{p} = SurfaceSph(shape_gallery(p,Shape)); 
   parbd0 = RBS_set_params(p,[0 0 0],1,...
        'TSL_Stk_3D',Shape,Sc0,sdim,Fparams.denseMV,Fparams.parbd.doAna,mdist,Fparams.parbd.eps,out); 
   Fpar0 = Fparams; Fpar0.parbd=parbd0; 
   % Single layer traction
   [~,Kernels.TSSDd] = RBS_Build_DMV('dense',Fpar0,0.5,'TSL_Stk_3D'); 
   
   Kernels.TSSD0 = Kernels.TSSDd;
   Kernels.ITSSD0 = inv(Kernels.TSSDd+Lk(1:Nb,1:Nb));
   % Single layer kernel
   [~,Kernels.SSDd] = RBS_Build_DMV('dense',Fpar0,0,'SL_Stk_3D');
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
       [~,Kernels.SLD0] = RBS_Build_DMV('dense',Fpar0,0,'SL_L_3D');

       % KL = (1/2 I + eta S')
       [~,Kernels.KLD0] = RBS_Build_DMV('dense',Fpar0,0.5,'dSL_L_3D');
   end
end

bld_bkdiag = (i>0) && ((~strcmp(typeMV,'Vsh')) || (strcmp(typeMV,'Vsh') && ~strcmp(prtype(1:2),'Vs')));
     
if bld_bkdiag && i==1 %~strcmp(typeMV,'Vsh') && n3>1 && i==1
    Kernels.TSSDd = cell(n3,1); Kernels.SSDd = cell(n3,1); 
    if denseMV==0 
       Kernels.ITSSDd = cell(n3,1);  
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diagonal block update (rotation) for typeMV=Rbs 
tic; 
if bld_bkdiag
    
for k=1:n3
    rk = Fparams.parbd.rd(k); 
    
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
   
   if abs(rk-1)>0
       Kernels.SSDd = rk*Kernels.SSDd;
   end
end
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set preconditioner
if n3>1 %&& Fparams.denseMV==0
    if strcmp(prtype,'TT')        
        fprintf('\n TT preconditioner build'); 
        [parslv.prec,parslv.prev] = RBS_setprec(Kernels.ITSSDd,n3,prtype,p,Fparams.parbd,prev); 
    elseif strcmp(prtype,'bkdiag')
        fprintf('\n Block diagonal preconditioner build');
        parslv.prec = RBS_setprec(Kernels.ITSSDd,n3,prtype);  
    elseif strcmp(prtype(1:3),'Vsh')
        fprintf(['\n Sparse / TT ' prtype ' preconditioner build']);
        [parslv.prec,parslv.prev] = RBS_setprec(Lk,n3,prtype,p,Fparams.parbd,prev);
        %{
        if Fparams.denseMV==0
           %test
           Fpar=Fparams; Fpar.parbd.dense=1; 
           [precd,prevd] = RBS_setprec(Lk,n3,prtype,p,Fpar.parbd,prev);
           EK = abs(full(parslv.prev{1})-full(prevd{1})); 
           EI = abs(full(parslv.prev{2})-full(prevd{2})); 
           fprintf('\n Test preconditioner vs dense, |TTK-TTKd|=%1.2e,|TTI-TTId|=%1.2e'...
               ,norm(EK)/norm(full(prevd{1})),norm(EI)/norm(full(prevd{2})));
           pause; 
           parslv.prec=precd; parslv.prev=prevd; 
        end
        %}
    else
        parslv.prec=[]; 
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

%{
if ~Fparams.parbd.dense
   %MatVec tests
   parbdd = Fparams.parbd; parbdd.dense=1; 
   TDM = RBS_MatVec([],Lk,typeMV,parbdd,sdim,0.5,'TSL_Stk_3D',Kernels.TSSDd);
   SDM = RBS_MatVec([],[],typeMV,parbdd,sdim,0,'SL_Stk_3D',Kernels.TSSDd);
   
   V = rand(size(TDM,2),1); V=V./norm(V); 
   
   fprintf('\n T apply error = %1.2e',norm(TDM*V-Kernels.TD(V))); 
   fprintf('\n S apply error = %1.2e',norm(SDM*V-Kernels.SD(V)));
   pause; 
   Id = eye(size(TDM,2)); TDM2=zeros(size(Id)); 
   for i=1:size(Id,2)
      TDM2(:,i) = Kernels.TD(Id(:,i));  
   end
   fprintf('\n T mat error = %1.2e',norm(TDM-TDM2)/norm(TDM)); 
   pause; 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Target eval test
   Xtrg = [0 0 20]; Nrtrg = [0 0 1]; parbdd.dense=0; 
   parbdd.flag_pot = 'TSL_Stk_3D'; parbdd.kerd=3; parbdd.a=0.5;
   Y1 = VSh_MatVec_RB_trg(V,Xtrg,Nrtrg,parbdd); 
   parbdd.dense=1; 
   Y2 = VSh_MatVec_RB_trg(V,Xtrg,Nrtrg,parbdd);
   display(norm(Y1))
   display(norm(Y2))
   fprintf('\n T target apply error = %1.2e',norm(Y1-Y2));
   pause; 
   %Kernels.TD = @(V) TDM*V; 
end
%}
    
if strcmp(Fparams.type,'MHD')
% Laplace kernels
Kernels.SLD = RBS_MatVec([],[],typeMV,Fparams.parbdL,ldim,0,'SL_L_3D',Kernels.SLD0); 
Kernels.KLD = RBS_MatVec([],[],typeMV,Fparams.parbdL,ldim,0.5,'dSL_L_3D',Kernels.KLD0); 
end

fprintf('\n Time for kernel eval update: %e',toc) 
if i>0
timings.operator.offd(i) = toc;
end

end

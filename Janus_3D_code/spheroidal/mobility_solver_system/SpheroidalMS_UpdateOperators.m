function [Kernels,Nullsp,Fparams,timings] = SpheroidalMS_UpdateOperators(Xt,Ct,Mt,nrmW,Kernels,Fparams,timings,i)
%{
For each timestep of the mobility solver, we need to update the matvecs associated with
each layer potential (and also update the nullspace operator L_k).

Inputs:
Xt : discretization points of bodies at timestep t=i
Ct : centers of bodies at timestep t=i
Mt : rotation matrices for each body at timestep t=i
nrmW : norm of rotational velocity for each body
Kernels : matvec matrices for operators
Fparams : parameter struct associated with mobility problem
timings : timings struct for mobility problem (for debugging)
i : current timestep in mobility problem (starts at i=0)

Outputs:
Kernels : BIE operators (either matrix/matrix-free depending on params passed)
    Kernels.TSSD0
        seems to be used for old collision resolution
    Kernels.ITSSD0
Nullsp : nullspace completion terms
Fparams : simulation parameters; see spheroidal_mobility.m.
timings : timings struct for debugging purposes

%%%%%%%%%%%%%%%%
CODE ANNOTATIONS
%%%%%%%%%%%%%%%%
Xrp -> X rotated points
prtype -> block diagonal or tensor train (so format of data)?
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize variables to be used later
sdim=3; 
parslv = Fparams.parslv; 
prec=parslv.prec;
prev=parslv.prev;
prtype=parslv.prtype; 
typeMV = Fparams.typeMV;
denseMV = Fparams.denseMV; 

tic; 

np = Fparams.parbd.np; n3 = size(Ct,1); p = Fparams.parbd.p;
Nb = Fparams.parbd.Nb; 
mdist = Fparams.parbd.mdist; out = Fparams.parbd.out; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Need to update the point clouds of each body after applying rotations/translations

Shape = Fparams.parbd.Shape; 

if strcmp(typeMV,'SSph') || strcmp(Shape, '')
    Sc = Fparams.parbd.Sc;  
else % For non-spherical bodies...
    Sc = cell(p,n3);
    for k=1:n3
        Sc{p,k} = SurfaceSph(vec3d(Xt((1:np)+np*(k-1),:)));
    end    
end

% Update Fparams.parbd
Fparams.parbd = SpheroidalMS_set_params(p,Ct,Fparams.parbd.rd,...
        'TSL_Stk_3D',Shape,Sc,sdim,Fparams.denseMV,Fparams.parbd.doAna,mdist,Fparams.parbd.eps,out); 

Xt = Fparams.parbd.Xrp; % X rotated points

surface_update_time = toc;
fprintf('\n Time for surface update: %e',surface_update_time)
if i>0
    timings.operator.surf(i) = surface_update_time;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Kernels / MatVec info

% Here, Lk is the actual nullspace operator (i.e. its matvec).
[Ck, Bk, Dk, Lk] = Build_SpheroidalAuxMats(Fparams.parbd.Wg,Xt,[],np,n3); 
Nullsp.C = Ck; Nullsp.B = Bk; Nullsp.D = Dk; Nullsp.L = Lk;

% Initial block diag build at initial timestep
if i==0
    Sc0 = cell(p,1);
    Sc0{p} = SurfaceSph(shape_gallery(p,Shape)); 
    parbd0 = RBS_set_params(p, [0 0 0], 1,...
        'TSL_Stk_3D',Shape,Sc0,sdim,Fparams.denseMV,Fparams.parbd.doAna,mdist,Fparams.parbd.eps,out); 
    Fpar0 = Fparams; Fpar0.parbd=parbd0;

    % TSL: 0.5*I + TSL
    [~,Kernels.TSSDd] = SpheroidalMS_Build_DMV('dense',Fpar0,0.5,'TSL_Stk_3D'); 
    
    Kernels.TSSD0 = Kernels.TSSDd;
    Kernels.ITSSD0 = inv(Kernels.TSSDd+Lk(1:Nb,1:Nb));

    % SLP: SLP
    [~,Kernels.SSDd] = SpheroidalMS_Build_DMV('dense',Fpar0,0,'SL_Stk_3D');
    Kernels.SSD0 = Kernels.SSDd;
        
    if strcmp(prtype,'TT')
        Kernels.ITSSDd = Kernels.TSSDd+Lk(1:Nb,1:Nb);
        Kernels.ITSSD0 = Kernels.ITSSDd;
    elseif strcmp(prtype,'bkdiag')   
        Kernels.ITSSDd = Kernels.ITSSD0; 
    end
end

bld_bkdiag = (i>0) && ((~strcmp(typeMV,'Vsh')) || (strcmp(typeMV,'Vsh') && ~strcmp(prtype(1:2),'Vs')));
if bld_bkdiag && i==1
    Kernels.TSSDd = cell(n3,1); Kernels.SSDd = cell(n3,1); 
    if denseMV==0 
        Kernels.ITSSDd = cell(n3,1);  
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Diagonal block update (rotation) for typeMV=Rbs
% For spheroidal mobility solve, probably not necessary.
tic; 
if bld_bkdiag      
    for k=1:n3
        rk = Fparams.parbd.rd(k); 
        
    if n3>1
        if nrmW(k)>1e-10 % If the body is rotating
            Kernels.TSSDd{k} = Rotate_Operator(Kernels.TSSD0,Mt{k},np); 
            Kernels.SSDd{k} = Rotate_Operator(Kernels.SSD0,Mt{k},np);
                
            if Fparams.denseMV==0 % Invert inverse of TSL for FMM ?
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
if n3>1
    if strcmp(prtype,'TT')        
        fprintf('\n TT preconditioner build'); 
        [parslv.prec,parslv.prev] = RBS_setprec(Kernels.ITSSDd,n3,prtype,p,Fparams.parbd,prev); 
    elseif strcmp(prtype,'bkdiag')
        fprintf('\n Block diagonal preconditioner build');
        parslv.prec = RBS_setprec(Kernels.ITSSDd,n3,prtype);  
    elseif strcmp(prtype(1:3),'Vsh')
        fprintf(['\n Sparse / TT ' prtype ' preconditioner build']);
        [parslv.prec,parslv.prev] = RBS_setprec(Lk,n3,prtype,p,Fparams.parbd,prev);
    else
        parslv.prec=[]; 
    end
end

Fparams.parslv = parslv; 

fprintf('\n Time for diagonal block and prec rotation: %e',toc) 
if i>0
    timings.operator.diag(i) = toc; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update Stokes kernels
tic;

% Stokes kernels 
Kernels.TD = SpheroidalMS_MatVec([],Lk,typeMV,Fparams.parbd,sdim,0.5,'TSL_Stk_3D',Kernels.TSSDd);
Kernels.SD = SpheroidalMS_MatVec([],[],typeMV,Fparams.parbd,sdim,0,'SL_Stk_3D',Kernels.SSDd); 

fprintf('\n Time for kernel eval update: %e',toc) 
if i>0
    timings.operator.offd(i) = toc;
end

end
    
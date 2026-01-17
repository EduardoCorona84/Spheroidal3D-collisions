function [Kernels,Nullsp,Fparams,timings] = SpheroidalMS_UpdateOperators(Xt,Ct,Mt,nrmW,Kernels,Fparams,timings,i)
%{
For each timestep of the mobility solver, we need to update the matvecs associated with
each layer potential (and also update the nullspace operator L_k).

Inputs:
Xt : discretization points of bodies at timestep t=i
    probably not needed...
    TODO: clean up this argument.
Ct : centers of bodies at timestep t=i
Mt : rotation matrices for each body at timestep t=i
nrmW : norm of rotational velocity for each body
    probably not needed; only used for RBS method
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

% Body parameters
equ_radii = Fparams.parbd.equ_radii;
polar_radii = Fparams.parbd.polar_radii;
eps = Fparams.parbd.eps;

np = Fparams.parbd.np; n3 = size(Ct,1); p = Fparams.parbd.p;
Nb = Fparams.parbd.Nb; 
mdist = Fparams.parbd.mdist;
out = Fparams.parbd.out;
doAna = Fparams.parbd.doAna;
kerd = Fparams.parbd.kerd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Need to update the point clouds of each body after applying rotations/translations

% Update Fparams.parbd
Fparams.parbd = SpheroidalMS_set_params( ...
    equ_radii   = equ_radii, ...
    polar_radii = polar_radii, ...
    p           = p, ...
    C           = Ct, ...
    eps         = eps, ...
    mdist       = mdist, ...
    out         = out, ...
    doAna       = doAna, ...
    flag_pot    = "TSL_Stk_3D", ...
    kerd        = kerd, ...
    dense       = denseMV);

Xt = Fparams.parbd.Xrp; % X rotated points

surface_update_time = toc;
fprintf('Time for surface update: %e\n',surface_update_time)
if i>0
    timings.operator.surf(i) = surface_update_time;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Kernels / MatVec info

% Here, Lk is the actual nullspace operator (i.e. its matvec).
[Ck, Bk, Dk, Lk] = Build_SpheroidalAuxMats(Fparams.parbd.Wg,Xt,[],np,n3); 
Nullsp.C = Ck; Nullsp.B = Bk; Nullsp.D = Dk; Nullsp.L = Lk;

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

fprintf('Time for diagonal block and prec rotation: %e\n',toc) 
if i>0
    timings.operator.diag(i) = toc; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update Stokes kernels
tic;

% Stokes kernels 
Kernels.TD = SpheroidalMS_MatVec([],Lk,typeMV,Fparams.parbd,sdim,0.5,'TSL_Stk_3D');
Kernels.SD = SpheroidalMS_MatVec([],[],typeMV,Fparams.parbd,sdim,0,'SL_Stk_3D'); 

fprintf('Time for kernel eval update: %e\n',toc) 
if i>0
    timings.operator.offd(i) = toc;
end

end
    
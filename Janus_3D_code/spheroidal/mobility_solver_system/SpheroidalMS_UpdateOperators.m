function [Kernels,Nullsp,Fparams,timings] = SpheroidalMS_UpdateOperators(Xt,Ct,Mt,nrmW,Kernels,Fparams,timings,i)
%{
For each timestep of the mobility solver, we need to update the matvecs associated with
each layer potential (and also update the nullspace operator L_k).

Inputs:
    Xt : discretization points of bodies at timestep t=i
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
        Kernels.TDself0
    Nullsp : nullspace completion terms
    Fparams : simulation parameters; see spheroidal_mobility.m.
    timings : timings struct for debugging purposes
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize variables to be used later
sdim=3; ldim=1;
parslv = Fparams.parslv; 
if isfield(parslv,'precond_type')
    precond_type = parslv.precond_type;
elseif isfield(parslv,'prtype')
    % Backward compatibility with legacy naming used in RBS.
    precond_type = parslv.prtype;
else
    precond_type = '';
end
typeMV = Fparams.typeMV;
denseMV = Fparams.denseMV;

tic;

% Body parameters
equ_radii = Fparams.parbd.equ_radii;
polar_radii = Fparams.parbd.polar_radii;
collision_eps = Fparams.parbd.collision_eps;
np = Fparams.parbd.np; n3 = size(Ct,1); p = Fparams.parbd.p;
Nb = Fparams.parbd.Nb; 
mdist = Fparams.parbd.mdist;
doAna = Fparams.parbd.doAna;
kerd = Fparams.parbd.kerd;
bodydist = Fparams.parbd.bodydist;
if isfield(Fparams.parbd, 'tsl_dealiasing')
    tsl_dealiasing = Fparams.parbd.tsl_dealiasing;
else
    tsl_dealiasing = true;
end
if isfield(Fparams.parbd, 'tsl_dealiasing_pad')
    tsl_dealiasing_pad = Fparams.parbd.tsl_dealiasing_pad;
else
    tsl_dealiasing_pad = 4;
end
if isfield(Fparams.parbd, 'tsl_backend')
    tsl_backend = Fparams.parbd.tsl_backend;
else
    tsl_backend = 'spheroidal';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Need to update the point clouds of each body after applying rotations/translations

% Update Fparams.parbd
Fparams.parbd = SpheroidalMS_set_params( ...
    equ_radii   = equ_radii, ...
    polar_radii = polar_radii, ...
    p           = p, ...
    C           = Ct, ...
    collision_eps = collision_eps, ...
    mdist       = mdist, ...
    doAna       = doAna, ...
    flag_pot    = "TSL_Stk_3D", ...
    kerd        = kerd, ...
    dense       = denseMV, ...
    bodydist    = bodydist, ...
    MRot        = Mt, ...
    tsl_backend = tsl_backend, ...
    tsl_dealiasing = tsl_dealiasing, ...
    tsl_dealiasing_pad = tsl_dealiasing_pad ...
    );

Xt = Fparams.parbd.Xrp; % X rotated points

surface_update_time = toc;
fprintf('Time for surface update: %e\n',surface_update_time)
if i>0
    timings.operator.surf(i) = surface_update_time;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Kernels / MatVec info

% Here, Lk is the actual nullspace operator (i.e. its matvec).
[Ck, Bk, Dk, Ak, Lk] = Build_SpheroidalAuxMats(Fparams.parbd.Wg,Xt,[],np,n3); 
Nullsp.C = Ck; Nullsp.B = Bk; Nullsp.D = Dk; Nullsp.A = Ak; Nullsp.L = Lk;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update Stokes kernels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kernel_timer = tic;
tic
Kernels.SD = SpheroidalMS_MatVec([],[],typeMV,Fparams.parbd,sdim,0,'SL_Stk_3D'); 
fprintf('Time for SLP update: %e\n',toc);

tic
Kernels.TD = SpheroidalMS_MatVec([],Lk,typeMV,Fparams.parbd,sdim,0.5,'TSL_Stk_3D');
fprintf('Time for TSL update: %e\n',toc);

% For the initial timestep, the cache must be built for the modified Laplace, which may
% take a long time!
if strcmp(Fparams.type,'MHD')
    % Laplace kernels
    Kernels.SLD = SpheroidalMS_MatVec([],[],typeMV,Fparams.parbd,ldim,0,'SL_L_3D'); 
    Kernels.KLD = SpheroidalMS_MatVec([],[],typeMV,Fparams.parbd,ldim,0.5,'dSL_L_3D');
elseif strcmp(Fparams.type, 'JanusAmp')
    tic

    lap_params = Fparams.parbd;
    lap_params.lambda = Fparams.lambda;
    lap_params.dense = true;
    if isfield(Fparams, 'sphwv_mex_opts')
        lap_params.sphwv_mex_opts = Fparams.sphwv_mex_opts;
    end

    spheroidal_modified_radial_workspace('clear');

    % Kernels.SLMODD = SpheroidalMS_MatVec([],[],typeMV,lap_params,ldim,0,'SL_LMOD_3D');
    Kernels.DLMODD = SpheroidalMS_MatVec([],[],typeMV,lap_params,ldim,0.5,'DL_LMOD_3D');
    % Kernels.dSLMODD = SpheroidalMS_MatVec([],[],typeMV,lap_params,ldim,0,'dSL_LMOD_3D');
    Kernels.dDLMODD = SpheroidalMS_MatVec([],[],typeMV,lap_params,ldim,0,'dDL_LMOD_3D');

    % Compute hydrophilic label at each point in the local body frame.
    flabel = Fparams.boundary_label;
    surfacelabel = zeros(np,n3);
    init_dir = Fparams.init_dir;
    for body_idx=1:n3
       xind = (1:np) + np*(body_idx-1);
       Xbody = Fparams.parbd.Xrp(xind,:);
       Xbackrot = Xbody*Mt{body_idx};
       body_meta = struct( ...
           'shape_type', char(string(Fparams.parbd.shape_type(body_idx))), ...
           'equ_radius', Fparams.parbd.equ_radii(body_idx), ...
           'polar_radius', Fparams.parbd.polar_radii(body_idx) ...
           );
       surfacelabel(:,body_idx) = LOCAL_eval_boundary_label(flabel, Xbackrot, init_dir(body_idx,:), body_meta);
    end
    Fparams.SurfaceLabel = reshape(surfacelabel,[],1);

    fprintf('Time for JanusAmp update: %e\n', toc);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set/update preconditioner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
if n3>1 && ~denseMV && strcmp(precond_type,'bkdiag')
    if i==0 || ~isfield(Kernels,'TDself0')
        fprintf('Building reference TD self blocks for block diagonal preconditioner.\n');
        TDself0 = cell(n3,1);
        prev_inv = cell(n3,1);
    else
        TDself0 = Kernels.TDself0;
        if isfield(Kernels,'ITSSDd') && iscell(Kernels.ITSSDd) && numel(Kernels.ITSSDd)==n3
            prev_inv = Kernels.ITSSDd;
        else
            prev_inv = cell(n3,1);
        end
    end

    [ITSSDd, TDself0] = LOCAL_build_td_inverse_blocks(Fparams.parbd, Lk, Nb, n3, TDself0, nrmW, prev_inv);
    if i>0
        fprintf('Updating block diagonal preconditioner in current row frame.\n');
    end

    if i==0 || ~isfield(Kernels,'ITSSD0')
        Kernels.ITSSD0 = ITSSDd;
    end
    Kernels.TDself0 = TDself0;
    Kernels.ITSSDd = ITSSDd;
    [parslv.prec, ~] = SpheroidalMS_setprec(ITSSDd, n3, precond_type, p, Fparams.parbd, []);
    parslv.prev = [];
else
    parslv.prec = [];
    parslv.prev = [];
end
Fparams.parslv = parslv;

diag_update_time = toc;
fprintf('Time for diagonal block and preconditioner setup: %e\n',diag_update_time) 
if i>0
    timings.operator.diag(i) = diag_update_time; 
end

kernel_eval_time = toc(kernel_timer);
fprintf('Time for kernel eval update: %e\n',kernel_eval_time) 
if i>0
    timings.operator.offd(i) = kernel_eval_time;
end

end %% END MAIN FUNCTION

function vals = LOCAL_eval_boundary_label(flabel, X_body, init_dir, body_meta)
    label_arity = nargin(flabel);
    if label_arity == 2
        vals = flabel(X_body, init_dir);
    elseif label_arity == 3 || label_arity < 0
        vals = flabel(X_body, init_dir, body_meta);
    else
        error('Fparams.boundary_label must accept 2 args (X,init_dir) or 3 args (X,init_dir,meta).');
    end

    vals = vals(:);
end

function [ITSSDd, TDself0] = LOCAL_build_td_inverse_blocks(parbd, Lk, Nb, n3, TDself0, nrmW, prev_inv)
    ROT_TOL = 1e-10;
    Iblock = eye(Nb);
    ITSSDd = cell(n3,1);
    reuse_block = false(n3,1);

    if nargin>=6 && isnumeric(nrmW) && numel(nrmW)==n3
        reuse_block = abs(nrmW(:))<=ROT_TOL;
    end

    if nargin<5 || ~iscell(TDself0) || numel(TDself0)~=n3
        TDself0 = cell(n3,1);
    end

    if nargin>=7 && iscell(prev_inv) && numel(prev_inv)==n3
        if nargin>=6 && isnumeric(nrmW) && numel(nrmW)==n3
            reuse_block = abs(nrmW(:))<=ROT_TOL;
        end
    else
        prev_inv = cell(n3,1);
    end

    for k=1:n3
        if reuse_block(k) && isnumeric(prev_inv{k}) && isequal(size(prev_inv{k}), [Nb Nb])
            ITSSDd{k} = prev_inv{k};
            continue;
        end

        [Akk, TDself0{k}] = SpheroidalMS_BuildTDSelfBlock(parbd, Lk, k, TDself0{k});
        ITSSDd{k} = Akk\Iblock;
    end
end

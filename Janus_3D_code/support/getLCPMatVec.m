function [A,Abad] = getLCPMatVec(Fparams, F, Kernels, Nullsp, p, gmresTol, Ct, denseMV, debug)
% override p
if exist('p','var') && ~isempty(p)
    Fparams.parbd.p = p;
end
% override gmresTol
if exist("gmresTol","var") && ~isempty(gmresTol)
    Fparams.parslv.tol = gmresTol;
end
% override denseMV
if exist("denseMV","var") && ~isempty(denseMV)
    Fparams.parbd.C = Ct;
end
% override denseMV
if exist("denseMV","var") && ~isempty(denseMV)
    Fparams.denseMV = denseMV;
end
% check for the mvp
if ~exist("debug","var") || isempty(debug)
    debug = false;
end

if ~exist('Kernels','var') || ~exist('Nullsp','var')  || exist('p','var')
    if ~strcmpi(Fparams.typeMV, 'vsh')
        warning('DMV needs to be set and is not implemented');
    end
    typeMV = 'Vsh'; 
    flag_pot ='SL_Stk_3D';  
    kerd = Fparams.parbd.kerd; 
    Ct = Fparams.parbd.C;
    rd = Fparams.parbd.rd;
    p = Fparams.parbd.p;
    doAna = Fparams.parbd.doAna;
    mdist = Fparams.parbd.mdist;
    matVecEps = Fparams.parbd.eps;
    Shape = Fparams.parbd.Shape;
    out = Fparams.parbd.out;
    denseMV = Fparams.denseMV;
    n3 = size(Ct,1);
    Sc = cell(p,1);
    Sc{p} = SurfaceSph(shape_gallery(p,Shape));
    Fparams.parbd = RBS_set_params(p,Ct,rd,...
        flag_pot,Shape,Sc,kerd,denseMV, ... 
        doAna,mdist,matVecEps,out);
    np = Fparams.parbd.np; Xt = Fparams.parbd.Xrp; Wg = Fparams.parbd.Wg;
    % get the aux mats
    [Ck,Bk,~,Ak,Lk] = Build_AuxMats2(Wg,Xt,[],np,n3);
    % get the kernels
    a=0; DMV =[];
    SD = RBS_MatVec([],[],typeMV,Fparams.parbd,kerd,a,flag_pot,DMV);
    flag_pot='TSL_Stk_3D'; a = 0.5;
    TD = RBS_MatVec([],Lk,typeMV,Fparams.parbd,kerd,a,flag_pot,DMV);
else 
    % If the kernels are already provide then use them
    Bk = Nullsp.B; Ck = Nullsp.C; Ak = Nullsp.A; Lk = Nullsp.L; 
    SD = Kernels.SD; TD = Kernels.TD;
end

if Fparams.denseMV
    A_basic = real((F'*Ak)*SD*(TD\(Bk'*(Ck*(Bk'*F)))));
    if nargout >1
        Abad = real((F'*Ck)*SD*(TD\((Bk'*(Ck*(Bk'*F))))));
    end
else
    if strcmpi(Fparams.parslv.prLCP, 'bkdiag')
        warning('Not tested')
        S0 = @(x) reshape(Kernels.SSD0*(repmat(rd.',Nb,size(x,2)).*reshape(x,Nb,n3*size(x,2))),[],size(x,2)); 
        IT0 = @(x) reshape(Kernels.ITSSD0*reshape(x,Nb,n3*size(x,2)),[],size(x,2));
        A_basic = @(x) real(F.'*(Ck*(S0(-IT0(Lapp(TD,Bf(x))+Lk*Bf(x))+Bf(x)))));
    else
        % TODO implement preconditioner
        parslv = Fparams.parslv;
        parslv.prec = []; 
        BkF = @(x) (Bk.')*(F*x);
        TD = @(x) Lapp(TD, x);
        A_basic = @(x) real(F.'*(Ak*Lapp(SD,Lslv(TD,-Lapp(TD,BkF(x))+Lk*BkF(x),parslv)+BkF(x))));
    end
    if debug 
        typeMV = 'Vsh'; flag_pot ='SL_Stk_3D'; a=0; kerd = Fparams.parbd.kerd; DMV =[];
        SD = RBS_MatVec('Mat',[],typeMV,Fparams.parbd,kerd,a,flag_pot,DMV);
        flag_pot='TSL_Stk_3D'; a = 0.5;
        TD = RBS_MatVec('Mat',Lk,typeMV,Fparams.parbd,kerd,a,flag_pot,DMV);
        n = size(F,2);
        Amat = real((F'*Ak)*SD*(TD\(Bk'*(Ck*(Bk'*F)))));
        AA = eye(n);
        for i = 1:n
            AA(:,i) = A_basic(AA(:,i));
        end
        fprintf('relErr in MVP %.4g\n',norm(Amat-AA) / norm(Amat))
    end
    if nargout >1
        Abad = @(x) real(F.'*(Ck*Lapp(SD,Lslv(TD,-Lapp(TD,BkF(x))+Lk*BkF(x),parslv)+BkF(x))));
    end
    A = @(x) A_verbose(x, A_basic, true);
end

function y = A_verbose(x, A_basic, verbose)
    start = tic;
    y = A_basic(x);
    stop = toc(start);
    if verbose 
        fprintf('-- A[x] MVP time %.3g sec \n', stop)
    end
end
function [A,Abad] = getLCPMatVec(Fparams, F, Kernels, Nullsp, p, gmresTol, Ct, denseMV)
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
% If the null space has already been computed then use it
if ~exist("Nullsp","var") ||isempty(Nullsp)
    np = Fparams.parbd.np; Xt = Fparams.parbd.Xrp; Wg = Fparams.parbd.Wg;
    [Ck,Bk,~,Ak,Lk] = Build_AuxMats2(Wg,Xt,[],np,n3);
else
    Bk = Nullsp.B; Ck = Nullsp.C; Ak = Nullsp.A; Lk = Nullsp.L; 
end

parslv = Fparams.parslv;
% If the kernels are already provide then use them
if ~exist('Kernels','var') || isempty('Kernels') || (exist('p','var') && ~isempty(p))
    if ~strcmpi(Fparams.typeMV, 'vsh')
        warning('DMV needs to be set and is not implemented');
    end
    typeMV = 'Vsh'; flag_pot ='SL_Stk_3D'; a=0; kerd = Fparams.parbd.kerd; DMV =[];
    SD = RBS_MatVec([],[],typeMV,params,kerd,a,flag_pot,DMV);
    flag_pot='TSL_Stk_3D'; a = 0.5;
    TD = RBS_MatVec([],Lk,typeMV,Fparams.parbd,kerd,a,flag_pot,DMV);
else 
    SD = Kernels.SD; TD = Kernels.TD;
end

if Fparams.denseMV
    A = real((F'*Ak)*SD*(TD\(Bk'*(Ck*(Bk'*F)))));
    if nargout >1top 
        Abad = real((F'*Ck)*SD*(TD\((Bk'*(Ck*(Bk'*F))))));
    end
else
    if strcmpi(Fparams.parslv.prLCP, 'bkdiag')
        warning('Not tested')
        S0 = @(x) reshape(Kernels.SSD0*(repmat(rd.',Nb,size(x,2)).*reshape(x,Nb,n3*size(x,2))),[],size(x,2)); 
        IT0 = @(x) reshape(Kernels.ITSSD0*reshape(x,Nb,n3*size(x,2)),[],size(x,2));
        A = @(x) real(F.'*(Ck*(S0(-IT0(Lapp(TD,Bf(x))+Lk*Bf(x))+Bf(x)))));
    else
        % TODO implement preconditioner
        parslv.prec = []; 
        BkF = @(x) (Bk.')*(F*x);
        A = @(x) real(F.'*(Ak*Lapp(SD,Lslv(TD,-Lapp(TD,BkF(x))+Lk*BkF(x),parslv)+BkF(x))));
    end
    if nargout >1
        Abad = @(x) real(F.'*(Ck*Lapp(SD,Lslv(TD,-Lapp(TD,BkF(x))+Lk*BkF(x),parslv)+BkF(x))));
    end
end
function [A,Abad] = getMatVec(Fparams, F, C, Kernels, Nullsp, p, gmresTol, denseMV)
% override p
if exists('p','var') && ~isempty(p)
    Fparams.parsh.p = p;
end
% override gmresTol
if exist("gmresTol","var") && ~isempty(gmresTol)
    Fparams.parslv.tol = gmresTol;
end
% override denseMV
if exist("denseMV","var") && ~isempty(denseMV)
    Fparams.denseMV = denseMV;
end

sdim = 3;
if ~exist('Kernels','var') || isempty('Kernels') || (exists('p','var') && ~isempty(p))
    Vshparams = params.parsh;
    Vshparams.prec = []; Vshparams.flag_pot ='SL_Stk_3D'; Vshparams.a=0;
    SD = VSh_MatVec_RB2('Mat',[],Vshparams); 
    Vshparams.a=.5; Vshparams.flag_pot='TSL_Stk_3D'; Vshparams.a = 0.5;
    TD = VSh_MatVec_RB2('Mat',[],Vshparams); 
else 
    SD = Kernels.SD; TD = Kernels.TD;
end

if ~exist("Nullsp","var") ||isempty(Nullsp)
    B = Nullsp.B; C = Nullsp.C; A = Nullsp.A; L = Nullsp.L; 
else
    np = TD.np;Xt = TD.Xrp; Wg = TD.Wg;
    [C,B,~,A,L] = Build_AuxMats2(Wg,Xt,[],np,n3);
end

if denseMV
    SMat = VSh_MatVec_RB2('Mat',[],SD); 
    TMat = VSh_MatVec_RB2('Mat',L,TD); 
    M = (A*SMat)*(((TMat\B')*C)*B');
    % Mbad = (C*SMat)*(((TMat\B')*C)*B');
    A = F.'*M*F;
    % Abad = F.'*Mbad*F;
else
    % TODO: Preconditioner?
    prec = [];
    % Traction of Single Layer stokes
    TD = RBS_MatVec([],Lk,'Vsh',...% typeMV,
        TD,sdim,0.5,'TSL_Stk_3D',prec);
    SD = RBS_MatVec([],[],'Vsh',...% typeMV,
        SD,sdim,0,'SL_Stk_3D',prec);
    Bf = @(x) (B.')*(F*x);
    parslv.prec = [];
    A = @(x) real(F.'*(C*Lapp(SD,Lslv(TD,-Lapp(TD,Bf(x))+L*Bf(x),parslv)+Bf(x))));
    % Abad = @(x) real(F.'*(A*Lapp(SD,Lslv(TD,-Lapp(TD,Bf(x))+L*Bf(x),parslv)+Bf(x))));
end
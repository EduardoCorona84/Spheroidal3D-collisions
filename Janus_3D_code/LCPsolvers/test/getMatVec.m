function [A,Abad] = getMatVec(Fparams, F, C, p, gmresTol, denseMV)

parslv = Fparams.parslv;
sdim = 3;
rd = Fparams.parbd.rd;
doAna = Fparams.parbd.doAna;
mdist = Fparams.parbd.mdist;
matVecEps = Fparams.parbd.eps;
Shape = Fparams.parbd.Shape;
out = Fparams.parbd.out;
kerd=3;
n3 = size(C,1);
Sc = cell(p,1);% TODO CHANGE
Sc{p} = SurfaceSph(shape_gallery(p,Shape));

parS = RBS_set_params(p,C,rd,...
    'SL_Stk_3D',Shape,Sc,kerd,denseMV, ... 
    doAna,mdist,matVecEps,out);
parT = RBS_set_params(p,C,rd,...
    'TSL_Stk_3D',Shape,Sc,kerd,denseMV, ... 
    doAna,mdist,matVecEps,out);
np = parT.np;
Xt = parT.Xrp;
Wg = parT.Wg;
[C,B,~,A,L] = Build_AuxMats2(Wg,Xt,[],np,n3);

if denseMV
    SMat = VSh_MatVec_RB2('Mat',[],parS); 
    parT.a = 0.5; 
    TMat = VSh_MatVec_RB2('Mat',L,parT); 
    M = (A*SMat)*(((TMat\B')*C)*B');
    Mbad = (C*SMat)*(((TMat\B')*C)*B');
    A = F.'*M*F;
    Abad = F.'*Mbad*F;
else
    % TODO: Preconditioner?
    prec = [];
    if exist('gmresTol', 'var') && ~isempty(gmresTol)
        parslv.tol = gmresTol;
    end
    Traction of Single Layer stokes
    TD = RBS_MatVec([],Lk,'Vsh',...% typeMV,
        parT,sdim,0.5,'TSL_Stk_3D',prec);
    SD = RBS_MatVec([],[],'Vsh',...% typeMV,
        parS,sdim,0,'SL_Stk_3D',prec);
    Bf = @(x) (B.')*(F*x);
    parslv.prec = [];
    A = @(x) real(F.'*(C*Lapp(SD,Lslv(TD,-Lapp(TD,Bf(x))+L*Bf(x),parslv)+Bf(x))));
end
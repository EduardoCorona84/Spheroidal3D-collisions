function A = getMatVec(Fparams, F, C, p)

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
matVecParams = RBS_set_params(p,C,rd,...
    'TSL_Stk_3D',Shape,Sc,kerd,false, ... %denseMV
    doAna,mdist,matVecEps,out);
np = matVecParams.np;
Xt = matVecParams.Xrp;
Wg = matVecParams.Wg;
[Ck,Bk,~, Lk] = Build_AuxMats(Wg,Xt,[],np,n3);
% TODO: Preconditioner?
prec = [];
% Traction of Single Layer stokes
TD = RBS_MatVec([],Lk,'Vsh',...% typeMV,
    matVecParams,sdim,0.5,'TSL_Stk_3D',prec);
SD = RBS_MatVec([],[],'Vsh',...% typeMV,
    matVecParams,sdim,0,'SL_Stk_3D',prec);
Bf = @(x) (Bk.')*(F*x);
parslv.prec = [];
A = @(x) real(F.'*(Ck*Lapp(SD,Lslv(TD,-Lapp(TD,Bf(x))+Lk*Bf(x),parslv)+Bf(x))));
function [Serr1,Serr2,Serr3,Mc1,Mc2,Mc3] = Test_BIEmatsym(pv,Ct)

if ~exist('pv','var') || isempty(pv)
    pv = 2:7;
end
if ~exist('Ct','var') || isempty(Ct)
    Ct = init_lattice(2,2.5,1,0); % n,cDist,rd,polydisperseRatio
end

nC=size(Ct,1); 
m = length(pv); 
np = @(p) 2*p.*(p+1); 

Serr1=zeros(m,1); Serr2 = Serr1; Serr3 = Serr1; 
Mc1 = cell(m,1); Mc2 = cell(m,1); Mc3 = cell(m,1); 

for k=1:m
    p = pv(k); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setup 
    Sc = SurfaceSph(shape_gallery(p,''));
    Xp = reshape(Sc.cart.to_array,[],3);
    SMat = kernelS([],Sc);
    TMat = kerneldS_RI(Sc,[]);

    [~, gwt]=g_grid(p+1);
    wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
    wt = wt(:);
    WT = Sc.geoProp.W; WT= WT.*wt;
    Nr = reshape(Sc.geoProp.nor.to_array,[],3);
    Xv = reshape(repmat(Xp,1,3)',3,[])';
    Nrv = reshape(repmat(Nr,1,3)',3,[])';
    Wv = repmat(WT,1,3)'; Wv = Wv(:); 

    par = Kernel_Eval_parameters('SL_Stk_3D',0,1,1,1,1e-8,2,400,1);
    par.dim = 3; par.mu=1;

    XpC = repmat(Xp,nC,1) + reshape(repmat(Ct',np(p),1),3,[])'; 
    XvC = repmat(Xv,nC,1) + reshape(repmat(Ct',3*np(p),1),3,[])';
    NrvC = repmat(Nrv,nC,1); 
    WvC = repmat(Wv,nC,1); 

    par.Xp = XpC; par.X = XvC; par.nor = NrvC; par.W2 = WvC';
    par.ci = repmat((1:3)',nC*np(p),1);
    par.cj = repmat((1:3)',nC*np(p),1);
    parT = par; parT.flag_pot='TSL_Stk_3D';

    [Ck,Bk,Dk,Ak,Lk] = Build_AuxMats2(repmat(WT(:),nC,1),XpC,Ct,np(p),nC);

    if 3*nC*np(p)<1000
        AM = struct('C',full(Ck),'B',full(Bk),'D',full(Dk),'A',full(Ak),'L',full(Lk));
    else
        AM = struct('C',Ck,'B',Bk,'D',Dk,'A',Ak,'L',Lk);
    end

    Vshparams = par; 
    Vshparams.p = p; Vshparams.kerd=3; Vshparams.out=1; Vshparams.dense=1; 
    Vshparams.C = Ct; Vshparams.n3=nC; Vshparams.rd=1; Vshparams.mdist=3; 
    Vshparams.a=0; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n Computing M using sing + smooth quad: \n')
    M1 = LOCAL_GenMobilityMat_V1(SMat,TMat,par,parT,AM,nC,np(p)); 
    fprintf('Computing M using spharm + smooth quad (GE solve): \n')
    M2 = LOCAL_GenMobilityMat_V2(AM,Vshparams); 
    fprintf('Computing M using spharm + smooth quad (GMRES solve): \n')
    M3 = LOCAL_GenMobilityMat_V3(AM,Vshparams);

    Mc1{k} = M1; Mc2{k} = M2; Mc3{k} = M3; 

    Serr1(k) = max(max(abs(M1-M1')))./max(abs(M1(:))); 
    Serr2(k) = max(max(abs(M2-M2')))./max(abs(M2(:)));
    Serr3(k) = max(max(abs(M3-M3')))./max(abs(M3(:)));

    display(log10(Serr1(k)));
    display(log10(Serr2(k)));
    display(log10(Serr3(k)));

end

nn = 3*nC*np(pv); 

figure; 
plot(log10(nn),log10(Serr1),'-ob','DisplayName','Sing Quad + KE'); hold on; 
plot(log10(nn),log10(Serr2),'-xr','DisplayName','Spectral + KE mdist=3');
plot(log10(nn),log10(Serr3),'-.g','DisplayName','Spectral for all'); 
hold off; 

end

function M = LOCAL_GenMobilityMat_V1(SMat,TMat,parS,parT,AM,nC,np)
 
Xv = parS.X; 

SKE = Kernel_Eval(Xv,Xv,parS);
TKE = Kernel_Eval(Xv,Xv,parT); 

%ILam = diag(repmat((1./[4*pi 4*pi 4*pi 8*pi/3 8*pi/3 8*pi/3])',nC,1));

for i=1:nC
    ind = (3*np*(i-1)+1):(3*np*i); 
    idl = (6*(i-1)+1):(6*i);

    SKE(ind,ind)=SMat; 
    TKE(ind,ind)=0.5*eye(3*np) + TMat + AM.B(idl,ind)'*AM.C(idl,ind);  %AM.L(ind,ind); 
end

M = (AM.A*SKE)*(((TKE\AM.B')*AM.C)*AM.B');

end

function M = LOCAL_GenMobilityMat_V2(AM,params,mdist)

if nargin>2
    params.mdist = mdist; 
end

%ILam = diag(repmat((1./[4*pi 4*pi 4*pi 8*pi/3 8*pi/3 8*pi/3])',params.n3,1));

SMat = VSh_MatVec_RB2('Mat',[],params); 
paramsT = params; paramsT.flag_pot='TSL_Stk_3D'; paramsT.a = 0.5; 
TMat = VSh_MatVec_RB2('Mat',AM.L,paramsT);

M = (AM.A*SMat)*(((TMat\AM.B')*AM.C)*AM.B');

end

function M = LOCAL_GenMobilityMat_V3(AM,params,mdist)

if nargin>2
    params.mdist = mdist; 
end

%ILam = diag(repmat((1./[4*pi 4*pi 4*pi 8*pi/3 8*pi/3 8*pi/3])',params.n3,1));

SMat = VSh_MatVec_RB2('Mat',[],params); 
paramsT = params; paramsT.flag_pot='TSL_Stk_3D'; paramsT.a = 0.5; 
TMat = VSh_MatVec_RB2('Mat',AM.L,paramsT);

TB = zeros(size(TMat,1),size(AM.B,1)); 
rs=1; TOL=1e-10; MaxIT = 100; 

for i=1:size(AM.B,1)
    [TB(:,i),~,~,iteri] = gmres(TMat,full(AM.B(i,:)'),rs,TOL,MaxIT);
    display(iteri)
end

M = (AM.A*SMat)*((TB*AM.C)*AM.B');

end
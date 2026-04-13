function Yt = RBS_FMM_Eval(Qt,W,kerd,pot,Xtrg,Xsrc,Nr)

d2=size(Qt,2); 
Yt=zeros(kerd*size(Xtrg,1),d2); 

for i=1:d2

Q = Qt(:,i); 

% source and target points variables
target = Xtrg.';     
ntarget = size(Xtrg,1);
if nargin<7
    Xsrc=Xtrg; 
end
source  = Xsrc.'; 
nsource = size(Xsrc,1);  

% right now only SL and dSL/DT are supported
sigma_sl = reshape(W.*Q,kerd,nsource); 
ifsingle=1; 
ifdouble=0; 
sigma_dl=zeros(kerd,nsource); 
sigma_dv=zeros(3,nsource); 
ifpot=1;  ifpottarg=0;
ifgradtarg=0;

if ~strcmp(pot(1:3),'SL_')
    ifgrad=1;
else
    ifgrad=0; 
end

% precision for FMM, roughly 3*iprec digits of acc
iprec=2; 

if kerd==1
% Laplace particle FMM 
U=lfmm3dpart(iprec,nsource,source,ifsingle,sigma_sl,ifdouble,sigma_dl,...
    sigma_dv,ifpot,ifgrad,ntarget,target,ifpottarg,ifgradtarg);  
else
% Stokes particle FMM 
U=stfmm3dpart(iprec,nsource,source,ifsingle,sigma_sl,ifdouble,sigma_dl,...
    sigma_dv,ifpot,ifgrad,ntarget,target,ifpottarg,ifgradtarg);
end

% Evaluate depending on pot
switch pot
    case 'SL_L_3D'
        % Single layer potential at targets
        Y    = (1/4/pi)*U.pot.'; 
    case 'dSL_L_3D'
        % compute du/dNrtrg
        GSF = -(1/4/pi)*U.fld; % Gradient, size 3 x ntarget
        Y = sum(GSF.*Nr.'); Y=Y(:); 
    case 'SL_Stk_3D'
        % Single layer potential at targets
        Y    = (1/4/pi)*U.pot.'; 
        Y = real(reshape(Y.',[],1));
    case 'TSL_Stk_3D'
        % Pressure
        SFpre = (1/4/pi)*U.pre; 
        % Gradient and Gradient transposed
        GSF   = (1/4/pi)*U.grad; 
        GTSF  = permute(GSF,[2 1 3]);  

        % Compute -pNr+Gu*Nr+Gut*Nr
        PNF = repmat(SFpre.',1,3).*Nr; 
        NrT = zeros(3,3,ntarget); NrT(:,1,:)=Nr.'; NrT(:,2,:)=Nr.'; NrT(:,3,:)=Nr.';
        GuN = reshape(sum(GSF.*NrT),[3 ntarget])+reshape(sum(GTSF.*NrT),[3 ntarget]); 

        Y  = -PNF.'+GuN; 
        Y = reshape(Y,[],1);    
    otherwise
        Y = zeros(ntarget,1); 
end

Yt(:,i)=Y; 

end

end
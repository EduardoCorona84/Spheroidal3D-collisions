function Y = SSph_FMM_Eval(Q, W, kerd, pot, Xtrg, Xsrc, Nr)
%{
Wrapper function to call the external FMM library.
Inputs
Q - (double) kerd*N_src × 1 column vector 
    source densities per DOF, ordered as follows:
    [q1(p1); q2(p1); ...; q3(p1); q1(p2); ...].
W - (double)  kerd*N_src × 1 column vector 
    quadrature weights aligned with Q (typically W2 from params)
kerd - (int) dimension of kernel (should be 1 or 3 for now)
pot - (string) list of potentials; see main function for description
Xtrg - (double)
Xsrc - (double)
Nr - (double) N_trg x 3 array
    target normals; required for certain potentials
%}

% source and target points variables
target = Xtrg.';     
ntarget = size(Xtrg,1);
if nargin<7
    Xsrc=Xtrg; 
end
source  = Xsrc.'; 
nsource = size(Xsrc,1);  

% Determine SL/DL and whether gradient (target normals) is needed
if strcmp(pot(1:2),'SL') || strcmp(pot(2:3),'SL')
    % SL and dSL/TSL
    sigma_sl = reshape(W.*Q,kerd,nsource); 
    ifsingle=1; 
    ifdouble=0; 
    sigma_dl=zeros(kerd,nsource); 
    sigma_dv=zeros(3,nsource); 
    ifpot=1;  ifpottarg=0;
    ifgradtarg=0;
else
    % DL and dDL/TDL
    sigma_dl = reshape(W.*Q,kerd,nsource); 
    ifsingle=0; 
    ifdouble=1; 
    sigma_sl=zeros(kerd,nsource); 
    sigma_dv=Nr.'; 
    ifpot=1;  ifpottarg=0;
    ifgradtarg=0;
end

if ~strcmp(pot(1:2),'SL') && ~strcmp(pot(1:2),'DL') 
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
    case 'DL_L_3D'
        % Double layer potential at targets
        Y    = (1/4/pi)*U.pot.'; 
    case 'SL_Stk_3D'
        % Single layer potential at targets
        Y    = (1/4/pi)*U.pot.'; 
        Y = real(reshape(Y.',[],1));
    case 'DL_Stk_3D'
        % Double layer potential at targets (check ct 1/4/pi)
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
    case 'TDL_Stk_3D'
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
        Y = zeros(ntarget,size(Q,2)); 
end
end
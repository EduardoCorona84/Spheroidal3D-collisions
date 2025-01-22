function [Dens_err,Coeff_err]=FarEvaluationTest(Mat,f,center,thetax,thetay)
%------------------------------------------------------------
% 
% Mat=layer potential (SL,DL)
% a=semi-major axis of ellipse
% f=function of (v,phi) surface density to test
% center=center of target spheroid
% thetax, thetay= rotation angles of target spheroid.
% 
%------------------------------------------------------------


% Initial parameters
%------------------------------------------------------------
u0 = 2/sqrt(3); a = sqrt(3)/2; Shape='ellipseZ'; 
params=SpheroidalParameters; params.u0=u0; params.a=a;
pstart=4; pend=16; parr=pstart:pend; % Select range of p-values to test
%------------------------------------------------------------


% Target points outside spheroid
%--------------------------------------------------------------
p_target=8; %order p for the target spheroid
Stns = SurfaceSph(shape_gallery(p_target,Shape)); 
Xtns_tmp = reshape(Stns.cart.to_array,[],3); %unshifted, unrotated

Rx=[1,0,0;0,cos(thetax),-sin(thetax);0,sin(thetax),cos(thetax)];
Ry=[cos(thetay),0,sin(thetay);0,1,0;-sin(thetay),0,cos(thetay)];

% Rotate and shift so that targets lie on a spheroid in the exterior
Xtns=zeros(size(Xtns_tmp));
for j=1:2*p_target*(p_target+1)
    Xtns(j,:)=((Xtns_tmp(j,:)-center)*Rx')*Ry';
end

%check that spheroids are not intersecting:
spheroidal_coord_targets=cart2spheroidal(Xtns,a);
ut=spheroidal_coord_targets(:,1);
%u coordinates should all be greater than u_0
if sum(ut<=u0)>0
    error('error: spheroids are intersecting.')
end
[~,ut_sort_index]=sort(ut);
%--------------------------------------------------------------


% Compute Layer potential using smooth quadrature
%--------------------------------------------------------------------
%--------------------------------------------------------------------
% Kernel Eval at highest order p, to compare with spheroidal layer
% potentials
Sns = SurfaceSph(shape_gallery(pend,Shape)); 
Xns = reshape(Sns.cart.to_array,[],3);

if strcmp(Mat,'SL')
    pot='SL_L_3D';
elseif strcmp(Mat,'DL')
    pot='DL_L_3D';
else
    pot='dSL_L_3D';
end

%Kernel parameters
pns = Kernel_Eval_parameters(pot,0,1,1,1,1e-8,2,400,1);
pns.dim = 3;
[~, gwt]=g_grid(pend+1);
wt = pi/pend*repmat(gwt', 2*pend, 1)./sin(gl_grid(pend));
wt = wt(:);
Wns = Sns.geoProp.W; Wns= Wns.*wt;
Nrns = reshape(Sns.geoProp.nor.to_array,[],3);
pns.X = Xns; pns.nor = Nrns; pns.W2 = Wns.';

[theta,phi]=gl_grid(pend);
sigmaKE=f(cos(theta),phi);

KYns = Kernel_Eval(Xtns,Xns,pns)*sigmaKE;
KYnssh = shAna(KYns);
%--------------------------------------------------------------------
%--------------------------------------------------------------------


% Now compute layer potentials over range of orders p=4:16
%--------------------------------------------------------------------
%--------------------------------------------------------------------
% Matrix to store error btwn KernelEval and Double Layer formula
Dens_err=zeros(length(parr),2*p_target*(p_target+1));
Coeff_err=zeros(length(parr),(p_target+1)^2);

for i=1:length(parr)

    % set current p value 
    p=parr(i);
    
    % Generate spheroidal surface (self)
    Sns = SurfaceSph(shape_gallery(p,Shape)); 
    Xns = reshape(Sns.cart.to_array,[],3);
    [theta,phi]=gl_grid(p);
    sigma=f(cos(theta),phi);

    params.sigma=sigma;
    params.get_shc;
    
    % Compute Layer potential 
    if strcmp(Mat,'SL')
        myL=spheroidalSL(params,Xtns);
        myLsh=shAna(myL);
    elseif strcmp(Mat,'DL')
        myL=spheroidalDL(params,Xtns);
        myLsh=shAna(myL);
    end
    
    %Compare with kernel eval 

%     this_err=(abs(myL-KYns)./abs(KYns))'; %relative error
    Dens_err(i,:)=(abs(myL-KYns))';
    Coeff_err(i,:)=(abs(myLsh-KYnssh))';
    
end
%--------------------------------------------------------------------
%--------------------------------------------------------------------

figure;
imagesc(log10(Dens_err(:,ut_sort_index)))
colorbar
ylabel('p')
xlabel('target point (increasing u)')
yticks(1:length(parr))
yticklabels(pstart:pend)
if strcmp(Mat,'DL')
    title('$D[\sigma]$ Error',Interpreter='latex')
elseif strcmp(Mat,'SL')
    title('$S[\sigma]$ Error',Interpreter='latex')
end

figure;
imagesc(log10(Coeff_err))
colorbar
ylabel('p')
xlabel('Coefficient index')
yticks(1:length(parr))
yticklabels(pstart:pend)
if strcmp(Mat,'DL')
    title('$D[\sigma]$ Spheroidal Harmonic Coefficients Error',Interpreter='latex')
elseif strcmp(Mat,'SL')
    title('$S[\sigma]$ Spheroidal Harmonic Coefficients Error',Interpreter='latex')
end

end

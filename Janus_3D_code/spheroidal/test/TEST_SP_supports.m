% tests for S'[]

% -------------------------------------------------------------------------
% test normal vectors: [test_normal.m]
clear;
pars=SpheroidalParameters; pars.u0 = 2/sqrt(3); pars.a = sqrt(3)/2;
pars.centers=[1,1,1];
Xvec = pars.get_X();
Xt = pars.get_X_targets([2,3,4]);
nor = pars.get_Norm();
pars.plot(); hold on;
quiver3(Xvec(:,1),Xvec(:,2),Xvec(:,3),nor(:,1),nor(:,2),nor(:,3)); hold off;
pause
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% test spheroidal - cartician conversion [TEST_nu_v.m]
%{
- Calculate formula for all n,m and compare;
!- Calculate formula for all u,v,phi, and compare;
- Use SL[] from two points, take normal to be the difference vector, and
approximate S' this way
!- S' formula for oblates
- Good results documented on notes, caption of what is done;
- Isolate code common to SPtest, DLtest_oblate, and SLtest_oblate, but not in prolate convergence tests.
 to see what is the issue -- guesses: solid harmonics, coefficients, 
%}

clear;
% pstart=2; pend=16; 
% parr=pstart:pend; spend=(pend+1)^2;
pend=10; spend=(pend+1)^2;
uu=2/sqrt(3);
a=1./uu;

pars=SpheroidalParameters;
pars.u0 = uu;
pars.a = a;
[up,vp]=gl_grid(pend);
% Yp=zeros(size(up,1),spend);
% for n=0:pend  %loop over terms in spheroidal harmonic expansion
%     Yn=Ynm(n,[],up,vp);
%     fac=1./sqrt(uu.^2-cos(up).^2); 
%     Yp(:,n^2+1:(n+1)^2)=Yn.*fac;
% end
n0=0;m0=0;nmind=n0^2+m0+n0+1;
Yin = Ynm(n0,m0,up,vp);
Yin = Yin ./sqrt(uu.^2-cos(up).^2);

pot='dSL_L_3D';
KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-8,2,400,1);
KEparams.dim = 3;
Xself=prolate_spheroid_shape(pend,uu,a);
Sns = SurfaceSph(Xself);
KEparams.X = Xself;
[~, gwt]=g_grid(pend+1);
wt = pi/pend*repmat(gwt', 2*pend, 1)./sin(gl_grid(pend));
wt = wt(:);
Wns = Sns.geoProp.W; Wns= Wns.*wt;
KEparams.W2 = Wns.';

dists=[2];
Xtrg=prolate_spheroid_shape(pend,uu,a);
Xtrg_dist = Xtrg+dists(1).*ones(size(Xtrg));
nu_sph=[0,1,0];
strg_dist=cart2spheroidal(Xtrg_dist(10,:),a); % for nu conversion
[nu_cart,~] = spheroidalNu2cart(nu_sph,strg_dist,a);
u_r=strg_dist(:,1); v_r=strg_dist(:,2); phi_r=strg_dist(:,3);
KEparams.nor = nu_cart;
dSmat = Kernel_Eval(Xtrg_dist(10,:),Xself,KEparams);
calcDY = dSmat * Yin;

pars.sigma=Yin;
mySP = spheroidalSP(pars,Xtrg_dist(10,:),nu_cart);
display(calcDY-mySP)
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% test components of S'[sigma] [TEST_ev_sp.m]
pend=16;
uu=2/sqrt(3); a=1/uu;
u_r=2+uu; v_r=-0.6; phi_r=pi/7;
Xtrg = [a*sqrt(u_r^2-1)*sqrt(1-v_r^2)*cos(phi_r),a*sqrt(u_r^2-1)*sqrt(1-v_r^2)*sin(phi_r),a*u_r*v_r]; 
strg = cart2spheroidal(Xtrg,a);
[u,v]=gl_grid(pend);
Yp = Ynm(1,1,u,v);
Yp = Yp ./sqrt(uu.^2-cos(u).^2);

nu21_sph = [1,0,0];
nu21_cart = spheroidalNu2cart(nu21_sph,strg,a);

nv2phi=(2*v_r*sqrt(u_r^2-1))/sqrt(v_r^2-u_r^2);
nu11_sph = [0,1,-1*nv2phi];
nu11_cart = spheroidalNu2cart(nu11_sph,strg,a);

% manual calculation of S'[Y11*fac]
% to check if implementation of the formula in notes is correct.
anm=-1/2*sqrt(uu^2-1);
PQ=legendre_otc(8,uu);
P=PQ{1}; 
p11uu=P(4); p21uu=P(8);
PQ2=legendre_otc(8,u_r,1,1,1);
Q=PQ2{2}; dQ=PQ2{4};
q11u=Q(4); q21u=dQ(8);
% % FOR [0,1,0]
% f1 =2*sqrt((1-v_r^2)/(u_r^2-v_r^2))*v_r/(1-v_r^2);
% f2 = sqrt((1-v_r^2)/(u_r^2-v_r^2))/(1-v_r^2);
% Ytemp = Ynm(1,1,acos(v_r),phi_r);
% Ytemp2=Ynm(2,1,acos(v_r),phi_r);
% Spmanual=anm*gnm*f1*fnm*Ytemp - anm*gnm*f2*fnm*Ytemp2;

% FOR [0,0,1]
% f1 = 1j/sqrt((u_r^2-1)*(1-v_r^2));
% Spmanual=anm*gnm*f1*fnm*Ytemp;

% FOR [0,1,-2vsqrt(u^2-1)/isqrt(u^2-v^2)] to isolate last component.
pars=SpheroidalParameters;
pars.u0 = uu;
pars.a = a;
pars.sigma=Yp;
dS11=spheroidalSP(pars,Xtrg,nu11_cart); % last component lambda

Y21_spectral = dS11*(2*sqrt(u_r^2-v_r^2)*sqrt(1-v_r^2))/(sqrt(uu^2-1)*p11uu*q11u);

% Kernel_Eval for SpMat * Y21
pot='dSL_L_3D';
KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-8,2,400,1);
KEparams.dim = 3;
Xself=prolate_spheroid_shape(pend,uu,a);
Sns = SurfaceSph(Xself);
KEparams.X = Xself;
[~, gwt]=g_grid(pend+1);
wt = pi/pend*repmat(gwt', 2*pend, 1)./sin(gl_grid(pend));
wt = wt(:);
Wns = Sns.geoProp.W; Wns= Wns.*wt;
KEparams.W2 = Wns.';
KEparams.nor = nu21_cart';
dSmat = Kernel_Eval(Xtrg,Xself,KEparams);
[u,v]=gl_grid(pend);
Yp21 = Ynm(2,1,u,v);
Yp21 = Yp21 ./sqrt(uu.^2-cos(u).^2);
dS21 = dSmat * Yp21;

Y21_kernel = dS21*(-6*sqrt(u_r^2-v_r^2))/(sqrt(uu^2-1)*sqrt(u_r^2-1)*p21uu*q21u);

display(Y21_spectral);
display(Y21_kernel);
display(Y21_true);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% test second component of S'[sigma] [TEST_bv.m]
clear;
p=8; np=2*p*(p+1); u0=2/sqrt(3); a=1/u0;
params=SpheroidalParameters;
params.u0=u0; params.a=a; params.isReal=true; 
[theta,phi]=gl_grid(p);
v=cos(theta);

Sns = SurfaceSph(shape_gallery(p,'ellipseZ')); 
Ys=[u0*ones(length(v),1) , v, phi];
Y=spheroidal2cart(Ys,a);
NrY=params.get_Norm;

% % CHECK1 Sns and Y should be same -- Good
Nrns = reshape(Sns.geoProp.nor.to_array,[],3);
% display(norm(abs(NrY-Nrns)));
% Ysns = reshape(Sns.cart.to_array,[],3);
% display(norm(abs(Ysns-Y)));

% Source
Xptch=[[0,0,0];[.1,.1,.1];[3,4,5]];
ptch=[1;2;4.5];
% Xptch=[0,0,0]; ptch =[1];

%% Using physics formuli
truesolnSurf=PtChargePotential(ptch,Xptch,Y); 
fluxSurf=PtChargeFlux(ptch,Xptch,Y,NrY);

%% Using Kernel Eval
pot='SL_L_3D';
KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-8,2,400,1);
KEparams.dim = 3;
Xself=reshape(Sns.cart.to_array,[],3);
KEparams.X = Xself;
KEparams.W2 = ones(size(ptch))'; 
SLmat = Kernel_Eval(Xself,Xptch,KEparams);

pot='dSL_L_3D';
KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-8,2,400,1);
KEparams.dim = 3;
Xself=reshape(Sns.cart.to_array,[],3);
KEparams.X = Xself;
KEparams.W2 = ones(size(ptch))'; 
KEparams.nor = Nrns; 
SPmat = Kernel_Eval(Xself,Xptch,KEparams);

surf_dir=SLmat*ptch;
surf_neu=SPmat*ptch;
display(norm(abs(surf_dir-truesolnSurf)))
display(norm(abs(surf_neu-fluxSurf)))

function pcp = PtChargePotential(ptch,Xptch,Y)
M=length(ptch);
np=length(Y);
Rptch=zeros(np,M);
for i=1:M
    Rptch(:,i)=sqrt(sum((Xptch(i,:)-Y).^2,2));
end
pcp=1./(4*pi*Rptch)*ptch; 
end

function flux=PtChargeFlux(ptch,Xptch,Y,NrY)
    M=length(ptch);
    np=length(Y);
    Rptch=zeros(np,M);
    EdotN=zeros(np,M);
    for i=1:M
        r=Xptch(i,:)-Y;
        Rptch(:,i)=sqrt(sum(r.^2,2));
        EdotN(:,i)=dot(1.*NrY,r./Rptch(:,i),2);
    end
    flux=EdotN./(4*pi*Rptch.^2)*ptch;
end

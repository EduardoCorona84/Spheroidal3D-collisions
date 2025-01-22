% check spectra calculated by Kernel_Eval.m

clear;
pend=8;sp=(pend+1)^2;uu=2/sqrt(3);a=1/sqrt(1+uu^2);
Xtrg=oblate_spheroid_shape(pend,uu,a);
Xtrg_dist=Xtrg+2.*ones(size(Xtrg));
Strg=cart2spheroidal(Xtrg_dist,a);
ux=Strg(:,1);v_r=Strg(:,2); phi_r=Strg(:,3);


[up,vp]=gl_grid(pend);
Y53src=Ynm(5,3,up,vp);
Yp=Y53src./sqrt(uu.^2+cos(up).^2);

c53=-2/factorial(8)*1j;
L0=legendre_otc(pend,1j.*uu);P0=L0{1};
L1=legendre_otc(pend,1j.*ux,1);Q1=L1{2};Fr=Q1.';
lambda53=c53*P0(34).*Fr(:,34);
Y53trg=Ynm(5,3,acos(v_r),phi_r);
manual_SL=lambda53.*Y53trg;

% % using spheroidalSL
% pars=SpheroidalParameters;
% pars.u0 = uu;
% pars.a = a;
% pars.oblate=1;
% pars.sigma=Yp;
% SLYp=spheroidalSL(pars,Xtrg_dist);
% display(norm(abs(SLYp-manual_SL)));

% Using Kernel_Eval
pot='SL_L_3D';
KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-8,2,400,1);
KEparams.dim = 3;
Xself=oblate_spheroid_shape(pend,uu,a);
Sns = SurfaceSph(Xself);
KEparams.X = Xself;
[~, gwt]=g_grid(pend+1);
wt = pi/pend*repmat(gwt', 2*pend, 1)./sin(gl_grid(pend));
wt = wt(:);
Wns = Sns.geoProp.W; Wns= Wns.*wt;
KEparams.W2 = Wns.';

SLmat = Kernel_Eval(Xtrg_dist,Xself,KEparams);
calcSY = SLmat * Yp;

display(norm(abs(calcSY-manual_SL)));


%{
%%%%%%%%%% CHECK for spectral convergence with higher p
pend=16;sp=(pend+1)^2;uu=2/sqrt(3);a=1/sqrt(1+uu^2);
[up,vp]=gl_grid(pend);
Y53=Ynm(5,3,up,vp);
Yp=Y53./sqrt(uu.^2+cos(up).^2);
c53=-2/factorial(8)*1j;
L0=legendre_otc(pend,1j.*uu);P0=L0{1};

Xtrg=oblate_spheroid_shape(pend,uu,a);
Xtrg_dist=Xtrg+2.*ones(size(Xtrg));
Strg=cart2spheroidal(Xtrg_dist,a);
ux=Strg(:,1);v_r=Strg(:,2); phi_r=Strg(:,3);

L1=legendre_otc(pend,1j.*ux,1);Q1=L1{2};Fr=Q1.';
lambda53=c53*P0(34).*Fr(:,34);
Y53p=Ynm(5,3,acos(v_r),phi_r);
manual_SL=lambda53.*Y53p;

% % using spheroidalSL
% pars=SpheroidalParameters;
% pars.u0 = uu;
% pars.a = a;
% pars.oblate=1;
% pars.sigma=Yp;
% SLYp=spheroidalSL(pars,Xtrg_dist);
% display(norm(abs(SLYp-manual_SL)));

% Using Kernel_Eval
pot='SL_L_3D';
KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-8,2,400,1);
KEparams.dim = 3;
Xself=oblate_spheroid_shape(pend,uu,a);
Sns = SurfaceSph(Xself);
KEparams.X = Xself;
[~, gwt]=g_grid(pend+1);
wt = pi/pend*repmat(gwt', 2*pend, 1)./sin(gl_grid(pend));
wt = wt(:);
Wns = Sns.geoProp.W; Wns= Wns.*wt;
KEparams.W2 = Wns.';

SLmat = Kernel_Eval(Xtrg_dist,Xself,KEparams);
calcSY = SLmat * Yp;

display(norm(abs(calcSY-manual_SL)));
%}


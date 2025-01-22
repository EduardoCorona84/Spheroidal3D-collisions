% Test solid harmonic calculation for oblate points (evaluated at imaginary
% inputs)

clear;
%% Solid_Harmonic checked for imaginary points OTC
% %{
p=4; uu=2/sqrt(3); a=1/sqrt(uu^2+1);
% u_x = uu + (1.1:0.5:2.5)';
u_x=uu-[0.1;0.3;0.5];
display(u_x);
% F = solid_harmonic(p,uu,u_x,1);
% for n=0:p
%     fprintf("\n n=%d\n", n);
%     display(F(:,(n)^2+1:(n+1)^2));
% end

PQ=legendre_otc(p,1j.*u_x,1,1,1);
P=PQ{1};Q=PQ{2};dP=PQ{3};dQ=PQ{4};
for n=0:p
    fprintf("\n n=%d\n", n);
    display(Q((n)^2+1:(n+1)^2,:));
end
% %}

%% DLspectrum checked for P(iu0), Q(iu0), and dP,dQ.
%{
p=4; uu=2/sqrt(3); a=1/uu;
u_x=uu;
[P,Q,dQ,dP,anm,lambda_surf]=DLspectrum(p,uu);
for n=0:p
    fprintf("\n n=%d", n);
    fprintf("\n dP val")
    display(dP((n)^2+1:(n+1)^2,:));
    fprintf("\n dQ val")
    display(dQ((n)^2+1:(n+1)^2,:));
    fprintf("\n anm val")
    display(anm((n)^2+1:(n+1)^2,:));
    % mm=(-n:n);
    % display(n-mm);
    % display(n+mm);
    % f1=factorial(n-mm); f2=factorial(n+mm);
    % a_manual=f1./f2.*(-1).^(mm+1);
    % a_manual_2=a_manual.*(uu^2+1);
    % fprintf("\n fact(n-m),fact(n+m),anm");
    % display(f1); display(f2); display(a_manual); display(a_manual_2)
    fprintf("\n lambda val")
    display(lambda_surf((n)^2+1:(n+1)^2,:));
    % fprintf("\n manual lambda for first point");
    % l_manual=a_manual_2(1)/2*(P(n^2+1)*dQ(n^2+1)+dP(n^2+1)*Q(n^2+1));
    % display(l_manual);
end
%}


%% SP spectrum vectorization check
% Vectorization CHECKED for all n and m.
%{
nu=[0,1,0];p=16;u0=2/sqrt(3);a=1/u0;sp=(p+1)^2;
Xself=prolate_spheroid_shape(p,u0,a);
Xtrg=prolate_spheroid_shape(p,u0,a);
Xtrg_dist=Xtrg+2.*ones(size(Xtrg));
Sself=cart2spheroidal(Xself,a);
Strg_centered=cart2spheroidal(Xtrg,a);
Strg_away=cart2spheroidal(Xtrg_dist,a);
u_x=Strg_away(:,1); v_x=Strg_away(:,2); phi_x=Strg_away(:,3);

% spectrum inside
nt_r=length(u_x);
Yr=zeros(nt_r*2,sp);
for n=0:p  %loop over terms in spheroidal harmonic expansion
    Yn=Ynm(n,[],acos(v_x)',phi_x);
    Yr(1:nt_r,n^2+1:(n+1)^2)=Yn;
    Yn1 = Ynm(n+1,-n:n,acos(v_x)',phi_x);
    Yr(nt_r+1:end,n^2+1:(n+1)^2) = Yn1;
end
[lambda_nm_prime, lambda_nm,lambda_n1m] = SPspectrum_away(p,u0,u_x,v_x,nu,0);
[F,Fp]=solid_harmonic_prime(p,u0,u_x,0);
Fmul_y11 = lambda_nm_prime.*Fp + lambda_nm.*F;
Fmul_y21= lambda_n1m.*F;
FYr=Fmul_y11.*Yr(1:nt_r,:)+Fmul_y21.*Yr(nt_r+1:end,:);

% manual
PQ=legendre_otc(p,u0);
PQ_r=legendre_otc(p,u_x,1,1,1);
P0=PQ{1}; Q1=PQ_r{2}; Qp1=PQ_r{4}; % Pnm(u0), Qnm(u_x), and Qnm'(u_x)
F=Q1.'; Fp=Qp1.';
Y_trgval=zeros(size(u_x,1),(p+2)^2); % size for one extra Yn+1m
for n=0:p+1  
    Yn=Ynm(n,[],acos(v_x),phi_x);
    Y_trgval(:,n^2+1:(n+1)^2)=Yn;
end
fac_u=zeros(size(u_x,1),sp); 
fac_v1=zeros(size(u_x,1),sp); 
fac_v2=zeros(size(u_x,1),sp);
fac_phi=zeros(size(u_x,1),sp); 
exp_u=zeros(size(u_x,1),sp); 
exp_v=zeros(size(u_x,1),sp); 
exp_phi=zeros(size(u_x,1),sp); 
for n=0:p
    for m=-n:n
        iind=n^2+m+n+1;
        common_fac=factorial(n-m)/factorial(n+m)*(-1)^m*sqrt(u0^2-1);
        gnm=common_fac.*P0(iind);
        normal_u_fac=sqrt((u_x.^2-1)./(u_x.^2-v_x.^2));
        normal_v_fac1=(n+1)./sqrt((u_x.^2-v_x.^2).*(1-v_x.^2)).*v_x;
        normal_v_fac2=-(n-m+1)./sqrt((u_x.^2-v_x.^2).*(1-v_x.^2));
        normal_phi_fac=1j*m./sqrt((u_x.^2-1).*(1-v_x.^2));
        fac_u(:,iind)=gnm.*normal_u_fac;
        fac_v1(:,iind)=gnm.*normal_v_fac1;
        fac_v2(:,iind)=gnm.*normal_v_fac2;
        fac_phi(:,iind)=gnm.*normal_phi_fac;
        exp_u(:,iind)=gnm.*normal_u_fac.*Fp(:,iind).*Y_trgval(:,iind);
        exp_v(:,iind)=gnm.*normal_v_fac1.*F(:,iind).*Y_trgval(:,iind)+gnm.*normal_v_fac2.*F(:,iind).*Y_trgval(:,(n+1)^2+m+(n+1)+1);
        exp_phi(:,iind)=gnm.*normal_phi_fac.*F(:,iind).*Y_trgval(:,iind);
    end
end

% with input
[u,v]=gl_grid(p);
Yin=zeros(size(u,1),sp);
for n=0:p
    Yn=Ynm(n,[],u,v);
    Yin(:,n^2+1:(n+1)^2)=Yn./sqrt(u0^2-cos(u).^2);
end
nu_sph=repmat(nu,size(u,1),1);
[nu_cart,~] = spheroidalNu2cart(nu_sph,Strg_away,a);
pars=SpheroidalParameters;
pars.u0 = u0;
pars.a = a;
pars.sigma=Yin;
spDY=spheroidalSP(pars,Xtrg_dist,nu_cart);

insideDY=FYr; % since Gshc should be identity.
manualDY=exp_v;
fprintf("\n insideSP-manualSP %f",norm(abs(insideDY-manualDY)));
fprintf("\n spheroidalSP-insideSP %f",norm(abs(insideDY-spDY))); 
% fprintf("\n insideSP-spheroidalSP, nm=00: %f",norm(insideDY(:,1)-spDY(:,1)));
% fprintf("\n maunalSP-spheroidalSP, nm=00: %f",norm(manualDY(:,1)-spDY(:,1)));
%}

%% Pnm' recursion
% % % TODO: Test Pnm'(v)=1/(1-v^2)*[(n+1)*v*Pnm(v)-(n-m+1)*P(n+1)m(v)]
% v_x=0.45;phi_x=pi/3.5;n=1;m=0;
% % ynm=Ynm(n,m,v_x,phi_x);yn1m=Ynm(n+1,m,v_x,phi_x);
% ynm=Ynm(n,m,acos(v_x),phi_x);yn1m=Ynm(n+1,m,acos(v_x),phi_x);
% recur=((n+1).*v_x.*ynm-(n-m+1).*yn1m)./(1-v_x.^2);
% display(recur)
% nu=[0,1,0];p=3;u0=2/sqrt(3);a=1/u0;sp=(p+1)^2;
% Xself=prolate_spheroid_shape(p,u0,a);
% Xtrg=prolate_spheroid_shape(p,u0,a);
% Xtrg_dist=Xtrg+2.*ones(size(Xtrg));
% Sself=cart2spheroidal(Xself,a);
% Strg_centered=cart2spheroidal(Xtrg,a);
% Strg_away=cart2spheroidal(Xtrg_dist,a);
% u_x=Strg_away(:,1); v_x=Strg_away(:,2); phi_x=Strg_away(:,3);
% 
% v_x=(-0.9999:0.3:1.0001)'; phi_x=Strg_away(1:size(v_x,1),3);
% nt_r=length(v_x);
% Yr=zeros(nt_r*2,sp);
% for n=0:p 
%     Yn=Ynm(n,[],acos(v_x)',phi_x);
%     Yr(1:nt_r,n^2+1:(n+1)^2)=Yn;
%     Yn1 = Ynm(n+1,-n:n,acos(v_x)',phi_x);
%     Yr(nt_r+1:end,n^2+1:(n+1)^2) = Yn1;
%     pp=legendre(n,v_x);
%     display(pp)
%     for m=-n:n
%         pnm=Yr(1:nt_r,n^2+m+n+1)./exp(1j*m.*phi_x);
%         pn1m=Yr(1+nt_r:end,n^2+m+n+1)./exp(1j*m.*phi_x);
%         display(pnm);
%         display(pn1m);
%         recur=((n+1).*v_x.*pnm-(n-m+1).*pn1m)./(1-v_x.^2);
%         % display(recur);
%     end
% end


%% SP spectrum with Kernel
%{
nu=[0,1,0];p=16;u0=2/sqrt(3);a=1/u0;sp=(p+1)^2;
Xself=prolate_spheroid_shape(p,u0,a);
Xtrg=prolate_spheroid_shape(p,u0,a);
Xtrg_dist=Xtrg+2.*ones(size(Xtrg));
Sself=cart2spheroidal(Xself,a);
Strg_centered=cart2spheroidal(Xtrg,a);
Strg_away=cart2spheroidal(Xtrg_dist,a);
u_x=Strg_away(:,1); v_x=Strg_away(:,2); phi_x=Strg_away(:,3);

[u,v]=gl_grid(p);
Yin=zeros(size(u,1),sp);
for n=0:p
    Yn=Ynm(n,[],u,v);
    % Yin(:,n^2+1:(n+1)^2)=Yn;
    Yin(:,n^2+1:(n+1)^2)=Yn./sqrt(u0^2-cos(u).^2);
end
nu_sph=repmat(nu,size(u,1),1);
[nu_cart,~] = spheroidalNu2cart(nu_sph,Strg_away,a);
pars=SpheroidalParameters;
pars.u0 = u0;
pars.a = a;
pars.sigma=Yin;
spDY=spheroidalSP(pars,Xtrg_dist,nu_cart);

pot='dSL_L_3D';
KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-8,2,400,1);
KEparams.dim = 3;
Sns = SurfaceSph(Xself);
KEparams.X = Xself;
[~, gwt]=g_grid(p+1);
wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
wt = wt(:);
Wns = Sns.geoProp.W; Wns= Wns.*wt;
KEparams.W2 = Wns.';
KEparams.nor=nu_cart;
dSmat=Kernel_Eval(Xtrg_dist,Xself,KEparams);
calcDY=dSmat*Yin;

fprintf("\n Kernel-spheroidalSP %f",norm(abs(calcDY-spDY)));
%}

function [P,Q,dQ,dP,anm,lambda_surf]=DLspectrum(p,u0)
    sp=(p+1)^2;
    ii = (1:sp)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;
	anm=factorial(nn-mm)./factorial(nn+mm).*(-1).^(mm+1).*(u0.^2+1);
    L=legendre_otc(p,1j.*u0,1,1,1);
    P=L{1}; Q=L{2}; dP=L{3}; dQ=L{4};
    lambda_surf=anm./2.*(P.*dQ+dP.*Q);
end

function F=solid_harmonic(p,u0,u_x,oblate)
    F=1;
    if oblate
	    u_x = 1j.*u_x;
    end
    if abs(u_x) - u0 < -1e-14
        PQ=legendre_otc(p,u_x,0);
        P=PQ{1};
        F=P.';
    elseif abs(u_x) - u0 > 1e-14
        PQ=legendre_otc(p,u_x,1);
        Q=PQ{2};
        F=Q.';
        % 
        PQ=legendre_otc(p,u_x,1,1,1);
        P=PQ{1};Q=PQ{2};dP=PQ{3};dQ=PQ{4};
        F=dQ.';
    end
end

function [lambda_nm_prime, lambda_nm,lambda_n1m] = SPspectrum_away(p,u0,u_x,v_x,nu,oblate)
    % nu = [nu_u,nu_v,nu_phi] Nx x 3 normal vector in spheroidal basis
    sp=(p+1)^2;
    ii = (1:sp)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1; 
    anm=factorial(nn-mm)./factorial(nn+mm).*(-1).^mm .*sqrt(u0.^2-1);
    if norm(abs(u_x)-u0)<1e-14
        fprintf("\n On surface, using SPspectrum instead.\n ");
        PQ= legendre_otc(p,u0,1,1,1);
        P=PQ{1}; Q=PQ{2}; dP=PQ{3}; dQ=PQ{4};
        lambda_nm_prime=zeros(size(u_x));
        lambda_n1m=zeros(size(u_x));
        lambda_nm=anm./2.*(P.*dQ+dP.*Q);
        return;
    end
    if abs(u_x)-u0<-1e-14
        PQ = legendre_otc(p,u0,1);
        gnm = PQ{2};
    elseif abs(u_x)-u0>1e-14
        PQ = legendre_otc(p,u0);
        gnm = PQ{1};
    else
        fprintf("\n Something went wrong.\n")
    end
    nu_u = nu(:,1); nu_v = nu(:,2); nu_phi = nu(:,3);

    lambda_nm_prime = anm.'.*gnm.'.*sqrt((u_x.^2-1)./(u_x.^2-v_x.^2)).*nu_u;
    lambda_nm = anm.'.*gnm.'.*((nn'+1).*v_x.*nu_v./sqrt((u_x.^2-v_x.^2).*(1-v_x.^2))+1j.*mm'.*nu_phi./sqrt((u_x.^2-1).*(1-v_x.^2)));
    lambda_n1m = -anm.'.*gnm.'.*(nn'-mm'+1).*nu_v./sqrt((u_x.^2-v_x.^2).*(1-v_x.^2));
end

function [F,Fp]=solid_harmonic_prime(p,u0,u_x,oblate)
    F=ones(1,(p+1)^2);
    Fp=ones(1,(p+1)^2);
    if nargin < 4
        oblate = false;
    end
    if oblate
        u_x = 1j.*u_x;
    end
    if abs(u_x)-u0 < -1e-14 
        PQ=legendre_otc(p,u_x,1,1);
        P=PQ{1}; dP=PQ{3};
        F=P.'; Fp=dP.';
    elseif abs(u_x)-u0 > 1e-14
        PQ=legendre_otc(p,u_x,1,1,1);
        Q=PQ{2}; dQ=PQ{4};
        F=Q.'; Fp=dQ.';
    end
end
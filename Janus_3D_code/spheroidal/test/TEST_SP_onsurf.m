ns=1;
uu=1.1; a=1/uu;
pars=SpheroidalParameters;
pars.u0 = uu;
pars.a = a;
pars.oblate=0; 

% ns=3;
% uu=[1.1,1.2,1.3];
% a=1./uu;
% pars=SpheroidalParameters;
% pars.u0=uu;
% pars.a=a;
% pars.oblate=[0 0 0];
% pars.centers = [0 0 0; 150 60 0; 130 -140 150];
% thetas = [0 pi/10 5*pi/3];
% phis = [0 0 pi/5];
% 
% Ri=zeros(3,3,3);
% for ii=1:ns_test
%     thetai=thetas(ii);
%     phii=phis(ii);
%     Riy=[cos(thetai) 0 sin(thetai); 0 1 0; -sin(thetai) 0 cos(thetai)];
%     Riz=[cos(phii) -sin(phii) 0; sin(phii) cos(phii) 0; 0 0 1];
%     Ri(:,:,ii)=Riz*Riy;
% end
% pars.Rmat=Ri;

p=8; np=2*p*(p+1);
pars.sigma=ones(np,1,ns); pars.get_shc();

nu_x=repmat([1,0,0],np,1,ns);
nu_y=repmat([0,1,0],np,1,ns);
nu_z=repmat([0,0,1],np,1,ns);

[X_self,~]=pars.get_X(); 

% 1) Implementation of [] input.
% TODO: test for multiple spheroids.
[LP1,LP2,LP3] = spheroidalSP(pars,X_self,nu_x,nu_y,nu_z);

[LP4,LP5,LP6] = spheroidalSP(pars,[],nu_x,nu_y,nu_z);

display(norm(abs([LP1;LP2;LP3]-[LP4;LP5;LP6])));
% CHECKED: negligible error for 1 spheroid.


% 2) When nu_x is e_u, should get back the surface outward normal result.
nu_norm_cart=zeros(np,3,ns);
for i=1:ns
    [~,X_self]=pars.get_X(i);
    Strg=cart2spheroidal(X_self,a(i),pars.oblate(i));
    [nu_norm_cart_i,~]=spheroidalNu2cart(nu_x(:,:,i),Strg,a(i),0);
    nu_norm_cart(:,:,i)=nu_norm_cart_i;
end
LP7 = spheroidalSP(pars,[],nu_norm_cart);

LP8 = spheroidalSP(pars);

display(norm(abs(LP7(:,:,1)-LP8(:,:,1))));
% CHECKED: negligible error for 3 spheroids.





% Testing script for Laplace to Stokes, manually calculate SL[sigma] and
% compare.
% sigma=[2,3,5; 2,3,5;...; 2,3,5]
% Xeval=[1,1,1]
clear;

params=SpheroidalParameters; 
params.matvec_eta=10; 
params.u0=1.1;
params.a=1/1.1;
params.oblate=0; 
params.centers=[0 0 0];
p=8; np=2*p*(p+1);

% sigma_x=2.*ones(np,1); 
% sigma_y=3.*ones(np,1); 
% sigma_z=15.*ones(np,1);
% Xeval=[1 1 1];
sigma_x=rand(np,1)+0.5;
sigma_y=rand(np,1)-0.5;
sigma_z=rand(np,1);
Xtrg=prolate_spheroid_shape(p,1.2,1/1.2)+[3,3,1];

Xself=prolate_spheroid_shape(p,params.u0,params.a);
Sns=SurfaceSph(Xself);

for trgind=1:size(Xtrg,1)
    % was having trouble with Rvec for multiple target points. So for
    % testing, loop over different test points. However, no restrictions on
    % test points for actual L2Stk function, so can vectorize.
    
    Xeval=Xtrg(trgind,:);
    Rvec=[Xeval(1)-Xself(:,1) Xeval(2)-Xself(:,2) Xeval(3)-Xself(:,3)];
    r=(vecnorm(Rvec'))';

    fprintf('\n--- Checking Laplace potentials... ---');
    %% first check Laplace
    % %{
    sigma_list=[sigma_x sigma_y sigma_z];
    Imat=eye(3);
    for sigma_ind=1:3
        sig=sigma_list(:,sigma_ind);
        params.sigma=sig;
        params.get_shc;
    
        Sl_sig=1/4/pi.*integrateOverS(Sns,sig./r);
        SL_matvec=spheroidalMatVec(params,'SL',Xeval);
        fprintf("\n SL sigma_i vs matvec inf error %e", max(abs(Sl_sig-SL_matvec)));
    
        for j=1:3
            dSl_sig=-1/4/pi.*integrateOverS(Sns,sig.*Rvec(:,j)./(r.^3));
            dSL_matvec=spheroidalMatVec(params,'SP',Xeval,Imat(j,:));
            fprintf("\n SP sigma_i dxj inf error %e",max(abs(dSl_sig-dSL_matvec)));
        end
    end
    % %}
    fprintf('\n-------------- Done. ---------------------\n');

    %% Next check Stokes
    % %{
    % i=1
    F1=sigma_x./r;
    
    F11x=sigma_x.*Rvec(:,1).*Rvec(:,1)./(r.^3);
    F12y=sigma_y.*Rvec(:,1).*Rvec(:,2)./(r.^3);
    F13z=sigma_z.*Rvec(:,1).*Rvec(:,3)./(r.^3);
    
    SLx=1/8/pi.*integrateOverS(Sns,F1+F11x+F12y+F13z);
    
    % i=2
    F1=sigma_y./r;
    
    F12x=sigma_x.*Rvec(:,2).*Rvec(:,1)./(r.^3);
    F22y=sigma_y.*Rvec(:,2).*Rvec(:,2)./(r.^3);
    F23z=sigma_z.*Rvec(:,2).*Rvec(:,3)./(r.^3);
    
    Sns=SurfaceSph(Xself);
    SLy=1/8/pi.*integrateOverS(Sns,F1+F12x+F22y+F23z);
    
    % i=3
    F1=sigma_z./r;
    
    F13x=sigma_x.*Rvec(:,3).*Rvec(:,1)./(r.^3);
    F23y=sigma_y.*Rvec(:,3).*Rvec(:,2)./(r.^3);
    F33z=sigma_z.*Rvec(:,3).*Rvec(:,3)./(r.^3);
    
    Sns=SurfaceSph(Xself);
    SLz=1/8/pi.*integrateOverS(Sns,F1+F13x+F23y+F33z);
    
    % calculated from laplace operators.
    [theta,phi]=gl_grid(p);
    % v=cos(theta);
    [s1,s2,s3]=L2Stk(Xeval,params,sigma_x,sigma_y,sigma_z,1);
    
    fprintf("\n inf error Stokes in x: %e, in y: %e, in z: %e\n", max(abs(s1-SLx)),max(abs(s2-SLy)),max(abs(s3-SLz)));
    % %}
end
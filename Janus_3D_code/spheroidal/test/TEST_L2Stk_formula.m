% Test Laplace to Stokes formula by forming each with Kernel_Eval.
clear;
p_test=8; np_test=2*p_test*(p_test+1);

% % Set up --- 1 spheroid
% ns_test=1;
% pars_test=SpheroidalParameters;
% pars_test.u0=1.1; pars_test.a=1/1.1;
% pars_test.p=p_test; pars_test.oblate=0;
% pars_test.centers = [0 0 0];
% pars_test.Rmat=eye(3);

% Set up ---- 3 spheroids
ns_test=3;
u0=[1.1 1.2 1.3];
pars_test=SpheroidalParameters;
pars_test.matvec_eta = 2;
pars_test.isReal=0;
pars_test.u0=u0;
pars_test.p=p_test;
alist = 1./u0;
pars_test.a=alist;
pars_test.oblate=[0 0 0];

pars_test.centers = [0 0 0; 150 60 0; 130 -140 150];
thetas = [0 pi/10 5*pi/3];
phis = [0 0 pi/5];

Ri=zeros(3,3,3);
for ii=1:3
    thetai=thetas(ii);
    phii=phis(ii);
    Riy=[cos(thetai) 0 sin(thetai); 0 1 0; -sin(thetai) 0 cos(thetai)];
    Riz=[cos(phii) -sin(phii) 0; sin(phii) cos(phii) 0; 0 0 1];
    Ri(:,:,ii)=Riz*Riy;
end
pars_test.Rmat=Ri;

pars_test.sigma=zeros(np_test,1,ns_test);
pars_test.get_shc();
    
% Point charges
rng(10.5);
nc = 2;
c=pars_test.centers;
Xptch = reshape(repmat(reshape(c',3,1,[]),1,nc),3,[],1)';
% random point charges
scale = repmat(.5.*pars_test.a .*sqrt(pars_test.u0.^2-1),nc,1);
scale = scale(:);
d = repmat(scale,1,3).*(rand(size(Xptch))-.5);

Xptch = Xptch + d; %point charge locations, near centers of spheroids
ptch = (2.*rand(ns_test*nc,1)-1); %charge value, between -1 and 1

% E field from point charges on surface of each spheroid
sigma1=zeros(np_test,1,ns_test); sigma2=sigma1; sigma3=sigma1;
[Xsrc,~]=pars_test.get_X();
for ii=1:ns_test
    Xsrci=Xsrc((ii-1)*np_test+1:ii*np_test,:);
    E=PtChargeE(ptch,Xptch,Xsrci); % nt x 3
    sigma1(:,:,ii)=E(:,1);
    sigma2(:,:,ii)=E(:,2);
    sigma3(:,:,ii)=E(:,3);
end


ydotsig=zeros(size(sigma1));

% src_x=Xsrc(:,1); src_y=Xsrc(:,2); src_z=Xsrc(:,3);
% for jj=1:ns_test
%     % ydotsig(:,:,jj)=sigma1(:,:,jj).*src_x((jj-1)*np_test+1:jj*np_test)+...
%     %     sigma2(:,:,jj).*src_y((jj-1)*np_test+1:jj*np_test)+...
%     %     sigma3(:,:,jj).*src_z((jj-1)*np_test+1:jj*np_test);
% end

% Target points off surface of each spheroid
pars_trg=copy(pars_test);
pars_trg.a=30.5.*pars_test.a; % set target points far enough away to avoid large error from quadrature.
[Xtrg,~]=pars_trg.get_X();

% % check target points ----------------------------------
% [Xsrc,~]=pars_test.get_X();
% for ii=1:ns_test
%     Xsrci=Xsrc((ii-1)*np_test+1:ii*np_test,:);
%     Xtrgi=Xtrg((ii-1)*np_test+1:ii*np_test,:);
%     plotb([Xsrci(:,1);Xsrci(:,2);Xsrci(:,3)]); hold on;
%     scatter3(Xtrgi(:,1),Xtrgi(:,2),Xtrgi(:,3));
% end
% hold off;
% % -------------------------------------------------------

Xtrg_test=pars_test.get_X_targets(Xtrg);

% %%%%%%%%%%%%%%%%%%%%
% Xtrg=[-3,0,-5;100,0,120];
% Xtrg=[(10:110)',(110:-1:10)',13.8.*ones(101,1)];
% Xtrg_test=pars_test.get_X_targets(Xtrg);
% %%%%%%%%%%%%%%%%%%%%

nt=size(Xtrg_test,1);

SLx=zeros(nt,1);SLy=SLx;SLz=SLx;
KELSLx=zeros(nt,1); KELSLy=KELSLx; KELSLz=KELSLx;

SPtemp=zeros(nt,1);

for ii=1:ns_test
    % Xself_ii=Xsrc((ii-1)*np_test+1:ii*np_test,:);
    Xself_centered=prolate_spheroid_shape(p_test,pars_test.u0(ii),pars_test.a(ii));
    Sns=SurfaceSph(Xself_centered);

    ydotsig(:,:,ii)=sigma1(:,:,ii).*Xself_centered(:,1)+...
        sigma2(:,:,ii).*Xself_centered(:,2)+...
        sigma3(:,:,ii).*Xself_centered(:,3);
    
    [~, gwt]=g_grid(p_test+1);
    wt = pi/p_test*repmat(gwt', 2*p_test, 1)./sin(gl_grid(p_test));
    wt = wt(:);
    Wns = Sns.geoProp.W; Wns= Wns.*wt;
    
    KEparams = Kernel_Eval_parameters('SL_Stk_3D',0,1,1,1,1e-8,2,400,1);
    KEparams.dim = 3; KEparams.mu=1;
    Wv = repmat(Wns,1,3)'; Wv=Wv(:);
    Xv = reshape(repmat(Xself_centered,1,3)',3,[])';
    KEparams.X = Xv;
    KEparams.W2 = Wv.';
    KEparams.cj=repmat((1:3)',np_test,1);
    KEparams.ci=repmat((1:3)',nt,1);
    Xtrg_ii=reshape(repmat(Xtrg_test(:,:,ii),1,3)',3,[])';
    LP_Kernel = Kernel_Eval(Xtrg_ii,Xv,KEparams);

    sig = reshape([sigma1(:,:,ii),sigma2(:,:,ii),sigma3(:,:,ii)].',[],1);
    LP = LP_Kernel*sig;
    LP=reshape(LP,3,[]).';
    SLx = SLx + LP(:,1);
    SLy = SLy + LP(:,2);
    SLz = SLz + LP(:,3); % CHECKED

    % Forming SL matrices to make Stk using Kernel_Eval to check KE and
    % formula for Stk.
    KEparams_1=Kernel_Eval_parameters('SL_L_3D',0,1,1,1,1e-8,2,400,1);
    KEparams_1.dim = 3; KEparams_1.mu=1;
    KEparams_1.X=Xself_centered;
    KEparams_1.W2=Wns.';
    SLker=Kernel_Eval(Xtrg_test(:,:,ii),Xself_centered,KEparams_1);

    KEparams_2=Kernel_Eval_parameters('dSL_L_3D',0,1,1,1,1e-8,2,400,1);
    KEparams_2.dim = 3; KEparams_2.mu=1;
    KEparams_2.X=Xself_centered;
    KEparams_2.W2=Wns.';
    Nu_x=repmat([1,0,0],nt,1);
    Nu_centered=(pars_test.Rmat(:,:,ii) * Nu_x')';
    KEparams_2.nor=Nu_centered;
    SPxker=Kernel_Eval(Xtrg_test(:,:,ii),Xself_centered,KEparams_2);

    KEparams_3=Kernel_Eval_parameters('dSL_L_3D',0,1,1,1,1e-8,2,400,1);
    KEparams_3.dim = 3; KEparams_3.mu=1;
    KEparams_3.X=Xself_centered;
    KEparams_3.W2=Wns.';
    Nu_y=repmat([0,1,0],nt,1);
    KEparams_3.nor=(pars_test.Rmat(:,:,ii) * Nu_y')';
    SPyker=Kernel_Eval(Xtrg_test(:,:,ii),Xself_centered,KEparams_3);

    KEparams_4=Kernel_Eval_parameters('dSL_L_3D',0,1,1,1,1e-8,2,400,1);
    KEparams_4.dim = 3; KEparams_4.mu=1;
    KEparams_4.X=Xself_centered;
    KEparams_4.W2=Wns.';
    Nu_z=repmat([0,0,1],nt,1);
    KEparams_4.nor=(pars_test.Rmat(:,:,ii) * Nu_z')';
    SPzker=Kernel_Eval(Xtrg_test(:,:,ii),Xself_centered,KEparams_4);

    SPx = SLker*sigma1(:,:,ii) + Xtrg_test(:,1,ii).*SPxker*sigma1(:,:,ii) + ...
        Xtrg_test(:,2,ii).*SPxker*sigma2(:,:,ii) + ...
        Xtrg_test(:,3,ii).*SPxker*sigma3(:,:,ii) + ...
        SPxker*ydotsig(:,:,ii);
    SPy = SLker*sigma2(:,:,ii) + Xtrg_test(:,1,ii).*SPyker*sigma1(:,:,ii) + ...
        Xtrg_test(:,2,ii).*SPyker*sigma2(:,:,ii) + ...
        Xtrg_test(:,3,ii).*SPyker*sigma3(:,:,ii) + ...
        SPyker*ydotsig(:,:,ii);
    SPz = SLker*sigma3(:,:,ii) + Xtrg_test(:,1,ii).*SPzker*sigma1(:,:,ii) + ...
        Xtrg_test(:,2,ii).*SPzker*sigma2(:,:,ii) + ...
        Xtrg_test(:,3,ii).*SPzker*sigma3(:,:,ii) + ...
        SPzker*ydotsig(:,:,ii);
    KELSLx=KELSLx + 1/2.*SPx; 
    KELSLy=KELSLy + 1/2.*SPy; 
    KELSLz=KELSLz + 1/2.*SPz;

end

display(sqrt(sum((SLx-KELSLx).^2)./size(SLx,1)))
display(sqrt(sum((SLy-KELSLy).^2)./size(SLy,1)))
display(sqrt(sum((SLz-KELSLz).^2)./size(SLz,1)))


function E=PtChargeE(ptch,Xptch,Y)
    M=length(ptch);
    if size(Y,2)~=3
        if size(Y,1)~=3
            error("one of target point's dimention should be 3")
        end
        Y = Y';
    end
    n=size(Y,1);
    Etemp=zeros(3*n,M);
    for ii=1:M
        r=Xptch(ii,:)-Y; % n x 3
        R=sqrt(sum(r.^2,2)); % n x 1
        Etemp(1:n,ii)=r(:,1)./(R.^3); % n x 1
        Etemp(n+1:2*n,ii)=r(:,2)./(R.^3);
        Etemp(2*n+1:3*n,ii)=r(:,3)./(R.^3);
    end
    E=Etemp./(4*pi)*ptch;
    E=reshape(E,n,3);
end
% ----------------------
% Test L2StkMatVec
% ----------------------------------------

pstart=2; pend=16;
parr=(pstart:2:pend)';

% Test 1: all far evaluation and compare with Kernel_Eval;
% Test 2: all near evaluation and compare wiht Kernel_Eval, on points far
% enough away. Potentially incorrect because of spectral method
% %{
% 1 spheroid:
% eta=0.1; % (all far)
% eta = 100; % (all near)
% 3 spheroids:
eta=1000; % (all near)
    
% Set up --- 1 spheroid
% ns_test=1;
% pars_test=SpheroidalParameters;
% pars_test.u0=1.1; pars_test.a=1/1.1;
% pars_test.oblate=0;
% pars_test.centers = [0 0 0];
% pars_test.Rmat=eye(3);
% pars_test.matvec_eta=eta;

% Set up ---- 3 spheroids 
ns_test=3;
u0=[1.1 1.2 1.3];
pars_test=SpheroidalParameters;
pars_test.matvec_eta = eta;
pars_test.isReal=0;
pars_test.u0=u0;
alist = 1./u0; 
alist(end)=1/sqrt(u0(end)^2+1);
pars_test.a=alist;
pars_test.oblate=[0 0 1];
% pars_test.oblate=[0 0 0];

pars_test.centers = [0 0 0; 150 60 0; 130 -140 150];
thetas = [0 pi/10 5*pi/3];
phis = [0 0 pi/5];

Ri=zeros(3,3,3);
for ii=1:ns_test
    thetai=thetas(ii);
    phii=phis(ii);
    Riy=[cos(thetai) 0 sin(thetai); 0 1 0; -sin(thetai) 0 cos(thetai)];
    Riz=[cos(phii) -sin(phii) 0; sin(phii) cos(phii) 0; 0 0 1];
    Ri(:,:,ii)=Riz*Riy;
end
pars_test.Rmat=Ri;

err=zeros(length(parr),1);
for pind=1:length(parr)
    p_test=parr(pind); 
    fprintf("\n Current p: %d\n",p_test);
    np_test=2*p_test*(p_test+1);

    pars_test.sigma=zeros(np_test,1,ns_test); pars_test.get_shc();
        
    % Point charges
    rng(10.5);
    nc = 2;
    c=pars_test.centers;
    Xptch = reshape(repmat(reshape(c',3,1,[]),1,nc),3,[],1)';
    % (pseudo)random point charges
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

    % Target points off surface of each spheroid
    pars_trg=copy(pars_test);
    pars_trg.a=30.5.*pars_test.a; % set target points far enough away to avoid large error from quadrature.
    [Xtrg,~]=pars_trg.get_X();

    % check target points ----------------------------------
    [Xsrc,~]=pars_test.get_X();
    for ii=1:ns_test
        Xsrci=Xsrc((ii-1)*np_test+1:ii*np_test,:);
        Xtrgi=Xtrg((ii-1)*np_test+1:ii*np_test,:);
        plotb([Xsrci(:,1);Xsrci(:,2);Xsrci(:,3)]); hold on;
        scatter3(Xtrgi(:,1),Xtrgi(:,2),Xtrgi(:,3));
    end
    hold off;
    % -------------------------------------------------------

    Xtrg_test=pars_test.get_X_targets(Xtrg);
    nt=size(Xtrg_test,1);

    SLx=zeros(nt,1);SLy=SLx;SLz=SLx;
    
    % tic
    for ii=1:ns_test
        if ~pars_test.oblate(ii)
            Xself_centered=prolate_spheroid_shape(p_test,pars_test.u0(ii),pars_test.a(ii));
        else
            Xself_centered=oblate_spheroid_shape(p_test,pars_test.u0(ii),pars_test.a(ii));
        end
        Sns=SurfaceSph(Xself_centered);

        ydotsig(:,:,ii)=sigma1(:,:,ii).*Xself_centered(:,1)+...
            sigma2(:,:,ii).*Xself_centered(:,2)+...
            sigma3(:,:,ii).*Xself_centered(:,3);
        
        KEparams = Kernel_Eval_parameters('SL_Stk_3D',0,1,1,1,1e-8,2,400,1);
        KEparams.dim = 3; KEparams.mu=1;
        [~, gwt]=g_grid(p_test+1);
        wt = pi/p_test*repmat(gwt', 2*p_test, 1)./sin(gl_grid(p_test));
        wt = wt(:);
        Wns = Sns.geoProp.W; Wns= Wns.*wt;
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
        SLz = SLz + LP(:,3);

    end
    % toc

    % tic
    [LSLx,LSLy,LSLz]=L2StkMatVec(pars_test,sigma1,sigma2,sigma3,Xtrg);
    LSLx=real(sum(LSLx,3));
    LSLy=real(sum(LSLy,3));
    LSLz=real(sum(LSLz,3));
    % toc
    full_diff=[SLx;SLy;SLz] - [LSLx;LSLy;LSLz];
    err(pind)=norm(sqrt(sum(full_diff.^2)./size(SLx,1)./3));
end

figure;
semilogy(parr,err);
xlabel("p");
ylabel("scaled 2-norm")
title("L2Stk (all spectral) error from Kernel\_Eval");
% %}


% Test 3: all near evaluation and compare with highest p, at points too
% close to the surface for Kernel_Eval.
% Test 4: mixed evaluation and compare with highest p, with eta=10.
% %{
eta=10; 
% 
% % Set up --- 1 spheroid
% ns_test=1;
% pars_test=SpheroidalParameters;
% pars_test.u0=1.1; pars_test.a=1/1.1;
% pars_test.oblate=0;
% pars_test.centers = [0 0 0];
% pars_test.Rmat=eye(3);
% pars_test.matvec_eta=eta;

% Set up ---- 3 spheroids
ns_test=3;
u0=[1.1 1.2 1.3];
pars_test=SpheroidalParameters;
pars_test.matvec_eta = 2;
pars_test.isReal=0;
pars_test.u0=u0;
alist = 1./u0; alist(end)=1/sqrt(u0(end)^2+1);
pars_test.a=alist;
pars_test.oblate=[0 0 1];

pars_test.centers = [0 0 0; 5 0 0; 1.6 1.6 1.6];
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

err=zeros(length(parr)-1,1);
for pind=length(parr):-1:1
    p_test=parr(pind); 
    fprintf("\n Current p: %d",p_test);
    np_test=2*p_test*(p_test+1);

    pars_test.sigma=zeros(np_test,1,ns_test); pars_test.get_shc();
        
    % Point charges
    rng(10.5);
    nc = 2;
    c=pars_test.centers;
    Xptch = reshape(repmat(reshape(c',3,1,[]),1,nc),3,[],1)';
    % (pseudo)random point charges
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

    if p_test==pend
        % Target points off surface of each spheroid
        pars_trg=copy(pars_test);
        pars_trg.a=2.5.*pars_test.a; 
        [Xtrg,~]=pars_trg.get_X();
    end


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
    nt=size(Xtrg_test,1);
    [LSLx,LSLy,LSLz]=L2StkMatVec(pars_test,sigma1,sigma2,sigma3,Xtrg);
    LSLx=real(sum(LSLx,3));
    LSLy=real(sum(LSLy,3));
    LSLz=real(sum(LSLz,3));

    if p_test==pend
        LSLx_p=LSLx;
        LSLy_p=LSLy;
        LSLz_p=LSLz;
    else
        full_diff=[LSLx_p;LSLy_p;LSLz_p]-[LSLx;LSLy;LSLz];
        err(pind)=norm(sqrt(sum(full_diff.^2)./size(LSLx,1)./3));
    end
end

figure;
semilogy([parr(1:end-1),parr(1:end-1),parr(1:end-1)],[err_x',err_y',err_z']);
legend("dim 1", "dim 2", "dim 3");
xlabel("p");
ylabel("scaled 2-norm")
title("L2Stk (all spectral) spectral error");
% %}



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
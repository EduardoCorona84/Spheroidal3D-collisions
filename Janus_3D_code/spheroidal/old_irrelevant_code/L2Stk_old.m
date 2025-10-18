function [Stk_x,Stk_y,Stk_z,varargout]=L2Stk_old(Xeval,pars,sigma_x,sigma_y,sigma_z,ns)
%--------------------------------------------------------------------%
% 
% %--------------------------------------------------------------------%

if nargin==0
    [SLx,SLy,SLz,Stk_x,Stk_y,Stk_z] =LOCAL_test_L2Stk();
    varargout{1}=SLx; varargout{2}=SLy; varargout{3}=SLz;
    return;
end


% i=1
pars.sigma=sigma_x; pars.get_shc;
SL1=spheroidalMatVec(pars,'SL',Xeval);
[SP1dx,SP1dy,SP1dz]=spheroidalMatVec_gradSP(pars,Xeval);

% i=2
pars.sigma=sigma_y; pars.get_shc;
SL2=spheroidalMatVec(pars,'SL',Xeval);
[SP2dx,SP2dy,SP2dz]=spheroidalMatVec_gradSP(pars,Xeval);

% i=3
pars.sigma=sigma_z; pars.get_shc;
SL3=spheroidalMatVec(pars,'SL',Xeval);
[SP3dx,SP3dy,SP3dz]=spheroidalMatVec_gradSP(pars,Xeval);

% extra term y dot sigma, for each i=1,2,3
Yself=pars.get_X();
y_src_x=Yself(:,1); y_src_y=Yself(:,2); y_src_z=Yself(:,3);
new_sig=zeros(size(sigma_x));
np=size(sigma_x,1);
for i=1:size(sigma_x,3)
    new_sig(:,:,i)=sigma_x(:,:,i).*y_src_x((i-1)*np+1:i*np)+sigma_y(:,:,i).*y_src_y((i-1)*np+1:i*np)+sigma_z(:,:,i).*y_src_z((i-1)*np+1:i*np);
end
pars.sigma=new_sig;
pars.get_shc;
[Fdx,Fdy,Fdz]=spheroidalMatVec_gradSP(pars,Xeval);

Stk_x=1/2.*(SL1-Xeval(:,1).*SP1dx-Xeval(:,2).*SP2dx-Xeval(:,3).*SP3dx+Fdx);
Stk_y=1/2.*(SL2-Xeval(:,1).*SP1dy-Xeval(:,2).*SP2dy-Xeval(:,3).*SP3dy+Fdy);
Stk_z=1/2.*(SL3-Xeval(:,1).*SP1dz-Xeval(:,2).*SP2dz-Xeval(:,3).*SP3dz+Fdz);


function [SLx,SLy,SLz,LSLx,LSLy,LSLz] =LOCAL_test_L2Stk()
    % Test Stokes operator made from Laplace operators on a system of 3
    % prolate spheroids, compared to Kernel_Eval for reference.

    p_test=8; np_test=2*p_test*(p_test+1);

    % Set up
    u0=[1.1 1.2 1.3];
    pars_test=SpheroidalParameters;
    pars_test.matvec_eta = 200;
    pars_test.isReal=0;
    pars_test.u0=u0;
    pars_test.p=p_test;
    alist = 1./u0;
    pars_test.a=alist;
    pars_test.oblate=[0 0 0];

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

    pars_test.sigma=zeros(np_test,1,3);
        
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
    ptch = (2.*rand(3*nc,1)-1); %charge value, between -1 and 1

    % E field from point charges on surface of each spheroid
    sigma1=zeros(np_test,1,3); sigma2=sigma1; sigma3=sigma1;
    [Xsrc,~]=pars_test.get_X();
    for ii=1:3
        Xsrci=Xsrc((ii-1)*np_test+1:ii*np_test,:);
        E=PtChargeE(ptch,Xptch,Xsrci); % nt x 3
        sigma1(:,:,ii)=E(:,1);
        sigma2(:,:,ii)=E(:,2);
        sigma3(:,:,ii)=E(:,3);
        % sigma_norm=sqrt(sum(E.^2,2));
        % plotb([Xsrci(:,1);Xsrci(:,2);Xsrci(:,3)],sigma_norm); hold on;
    end
    % hold off;

    % Target points off surface of each spheroid
    pars_trg=copy(pars_test);
    pars_trg.a=1.5.*pars_test.a;
    [Xtrg,~]=pars_trg.get_X();

    % % check target points ----------------------------------
    % [Xsrc,~]=pars_test.get_X();
    % for ii=1:3
    %     Xsrci=Xsrc((ii-1)*np_test+1:ii*np_test,:);
    %     Xtrgi=Xtrg((ii-1)*np_test+1:ii*np_test,:);
    %     plotb([Xsrci(:,1);Xsrci(:,2);Xsrci(:,3)]); hold on;
    %     scatter3(Xtrgi(:,1),Xtrgi(:,2),Xtrgi(:,3));
    % end
    % hold off;
    % % -------------------------------------------------------

    Xtrg_test=pars_test.get_X_targets(Xtrg);
    nt=size(Xtrg_test,1);

    SLx=zeros(nt,1);SLy=SLx;SLz=SLx;
    
    for ii=1:3
        Xself_ii=Xsrc((ii-1)*np_test+1:ii*np_test,:);
        Sns=SurfaceSph(Xself_ii);
        
        KEparams = Kernel_Eval_parameters('SL_Stk_3D',0,1,1,1,1e-8,2,400,1);
        KEparams.dim = 3; KEparams.mu=1;
        [~, gwt]=g_grid(p_test+1);
        wt = pi/p_test*repmat(gwt', 2*p_test, 1)./sin(gl_grid(p_test));
        wt = wt(:);
        Wns = Sns.geoProp.W; Wns= Wns.*wt;
        Wv = repmat(Wns,1,3)'; Wv=Wv(:);
        Xv = reshape(repmat(Xself_ii,1,3)',3,[])';
        KEparams.X = Xv;
        KEparams.W2 = Wv.';
        KEparams.cj=repmat((1:3)',np_test,1);
        KEparams.ci=repmat((1:3)',nt,1);
        Xtrg_ii=reshape(repmat(Xtrg_test(:,:,ii),1,3)',3,[])';
        LP_Kernel = Kernel_Eval(Xtrg_ii,Xv,KEparams);
        LP = LP_Kernel*[sigma1(:,:,ii);sigma2(:,:,ii);sigma3(:,:,ii)];

        SLx = SLx + LP(1:nt,:);
        SLy = SLy + LP(nt+1:2*nt,:);
        SLz = SLz + LP(2*nt+1:3*nt,:);
    end

    [LSLx_cell,LSLy_cell,LSLz_cell]=L2Stk_old(Xtrg,pars_test,sigma1,sigma2,sigma3,3);
    LSLx=sum(LSLx_cell,3);
    LSLy=sum(LSLy_cell,3);
    LSLz=sum(LSLz_cell,3);

    display(norm(abs(SLx-LSLx)))
    display(norm(abs(SLy-LSLy)))
    display(norm(abs(SLz-LSLz)))


end

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

end

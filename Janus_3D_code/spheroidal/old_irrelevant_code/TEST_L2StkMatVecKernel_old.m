% ----------------------
% Test L2StkMatVec
% ----------------------------------------

pstart=2; pend=16;
parr=(pstart:2:pend)';

% For one spheroid, compare SL matrix made by L2StkMatVecKernel
% with Kernel_Eval;
% Only self-eval necessary.
% Current code capacity: self-eval on one spheroid available; no self to
% all or interactions. Test code for one spheroid only, interchanging N-D
% array and 2D vectors sometimes.

eta=0.1; 
    
ns_test=1;
pars_test=SpheroidalParameters;
pars_test.u0=1.1; 
% pars_test.a=1/1.1;
% pars_test.oblate=0;
pars_test.a=1/sqrt(1.1^2+1);
pars_test.oblate=1;
pars_test.centers = [0 0 0];
pars_test.Rmat=eye(3);
pars_test.matvec_eta=eta;

err=zeros(1,length(parr));
for pind=1:length(parr)
    p_test=parr(pind); 
    fprintf("\n Current p: %d\n",p_test);
    np_test=2*p_test*(p_test+1);

    pars_test.sigma=zeros(np_test,1,ns_test); pars_test.get_shc();

    % E field from randomly placed point charges. 
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
    
    % target points on surface.
    % tic;
    [LSLx,LSLy,LSLz]=L2StkMatVec(pars_test,sigma1,sigma2,sigma3);
    % toc;
    LSLx=sum(LSLx,3);
    LSLy=sum(LSLy,3);
    LSLz=sum(LSLz,3);

    % Form kernel and solve for density
    % tic;
    StkKernel = L2StkMatVecKernel(pars_test,p_test);
    % toc;
    u_kernel=StkKernel*[sigma1;sigma2;sigma3];

    % TODO: because of kernel space, solution incorrect.
    %       turning to real does not solve the issue.
    % mu = StkKernel \ [LSLx;LSLy;LSLz];

    err(pind)=norm(abs(u_kernel-[LSLx;LSLy;LSLz]));
    % err(pind)=norm(abs(mu-[sigma1;sigma2;sigma3]));
end

display(err)


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
function [Stk_x,Stk_y,Stk_z]=L2Stk(Xeval,pars,sigma_x,sigma_y,sigma_z,ns)
    %{
        An implementation of the Laplace to Stokes single layer potential.

        Inputs
            Xeval       -   target points
            pars        -   parameters needed for to calculate spheroidal Laplace LPs
            sigma_x     -   vector
            sigma_y     -
            sigma_z     -
            ns          -   number of spheroidal bodies
    %}

    if nargin==0
        LOCAL_test_L2Stk();
        return;
    end

    % Self evaluation ----------------------------------------------
    if isempty(Xeval)
        % Setup for directional derivative
        nu_x_spectral=repmat([1,0,0],size(sigma_x,1),size(sigma_x,2),ns);
        nu_y_spectral=repmat([0,1,0],size(sigma_x,1),size(sigma_x,2),ns);
        nu_z_spectral=repmat([0,0,1],size(sigma_x,1),size(sigma_x,2),ns);

        % i=1; d/dx SL[sigma_x], etc.
        pars.sigma=sigma_x; pars.get_shc;
        SL1=spheroidalSL(pars); % SL self always returns N-D array
        [SP1dx,SP1dy,SP1dz] = spheroidalSP(pars,Xeval,nu_x_spectral,nu_y_spectral,nu_z_spectral);
        
        % i=2
        pars.sigma=sigma_y; pars.get_shc;
        SL2=spheroidalSL(pars);
        [SP2dx,SP2dy,SP2dz] = spheroidalSP(pars,Xeval,nu_x_spectral,nu_y_spectral,nu_z_spectral);
        
        % i=3
        pars.sigma=sigma_z; pars.get_shc;
        SL3=spheroidalSL(pars);
        [SP3dx,SP3dy,SP3dz] = spheroidalSP(pars,Xeval,nu_x_spectral,nu_y_spectral,nu_z_spectral);
    % All to targets
    else
        error_msg = "The evaluation points (i.e. the first argument) need to have type " + ...
            "'cell' of size 1 x ns. Each cell should contain the target points " + ...
            "for the associated body, and each cell should be of size nt x 3, where " + ...
            "nt is the number of target points on the body.\n" + ...
            "If you have a nt x 3 x ns matrix, then use squeeze(num2cell(Xeval, [1, 2])).' " + ...
            "to reduce back to a 1 x ns cell.";
        assert(isa(Xeval, "cell"), error_msg);
        assert(size(Xeval,1) == 1 && size(Xeval, 2) == ns, error_msg)
        nu_x_spectral=cell(1,ns); nu_y_spectral=nu_x_spectral; nu_z_spectral=nu_x_spectral; 
        for i=1:ns
            nu_x_spectral{i} = repmat([1,0,0],size(Xeval{i},1),1);
            nu_y_spectral{i} = repmat([0,1,0],size(Xeval{i},1),1);
            nu_z_spectral{i} = repmat([0,0,1],size(Xeval{i},1),1);
        end

        % i=1
        pars.sigma=sigma_x; pars.get_shc;
        SL1=spheroidalSL(pars,Xeval);
        [SP1dx,SP1dy,SP1dz] = spheroidalSP(pars,Xeval,nu_x_spectral,nu_y_spectral,nu_z_spectral);
        
        % i=2
        pars.sigma=sigma_y; pars.get_shc;
        SL2=spheroidalSL(pars,Xeval);
        [SP2dx,SP2dy,SP2dz] = spheroidalSP(pars,Xeval,nu_x_spectral,nu_y_spectral,nu_z_spectral);
        
        % i=3
        pars.sigma=sigma_z; pars.get_shc;
        SL3=spheroidalSL(pars,Xeval);
        [SP3dx,SP3dy,SP3dz] = spheroidalSP(pars,Xeval,nu_x_spectral,nu_y_spectral,nu_z_spectral);
    end

    % extra term y dot sigma, for each i=1,2,3
    new_sig=zeros(size(sigma_x));
    for i=1:ns
        if ~pars.oblate(i)
            Xloc=prolate_spheroid_shape(pars.p,pars.u0(i),pars.a(i));
        else
            Xloc=oblate_spheroid_shape(pars.p,pars.u0(i),pars.a(i));
        end
        new_sig(:,:,i)=sigma_x(:,:,i).*Xloc(:,1)+sigma_y(:,:,i).*Xloc(:,2)+sigma_z(:,:,i).*Xloc(:,3);
    end
    pars.sigma=new_sig;
    pars.get_shc();
    [Fdx,Fdy,Fdz] = spheroidalSP(pars,Xeval,nu_x_spectral,nu_y_spectral,nu_z_spectral);

    if isempty(Xeval)
        [Xloc,~]=pars.get_X();
        np_div=size(Xloc,1)/ns;
        np=round(np_div);
        if abs(np_div-np)>1e-8
            error("\n size of target on surface is not multiple of ns\n")
        end
        if isa(sigma_x,"cell")
            Stk_x=cell(1,ns); Stk_y=Stk_x; Stk_z=Stk_x;
            for i=1:ns
                Xloc_i=Xloc((i-1)*np+1:i*np,:);
                Stk_x{i}=1/2.*(SL1{i}-Xloc_i(:,1).*SP1dx{i}-Xloc_i(:,2).*SP2dx{i}-Xloc_i(:,3).*SP3dx{i}+Fdx{i});
                Stk_y{i}=1/2.*(SL2{i}-Xloc_i(:,1).*SP1dy{i}-Xloc_i(:,2).*SP2dy{i}-Xloc_i(:,3).*SP3dy{i}+Fdy{i});
                Stk_z{i}=1/2.*(SL3{i}-Xloc_i(:,1).*SP1dz{i}-Xloc_i(:,2).*SP2dz{i}-Xloc_i(:,3).*SP3dz{i}+Fdz{i});
            end
        else
            Stk_x=zeros(size(sigma_x,1),size(sigma_x,2),ns); Stk_y=Stk_x; Stk_z=Stk_x;
            for i=1:ns
                Xloc_i=Xloc((i-1)*np+1:i*np,:);
                Stk_x(:,:,i)=1/2.*(SL1(:,:,i)-Xloc_i(:,1).*SP1dx(:,:,i)-Xloc_i(:,2).*SP2dx(:,:,i)-Xloc_i(:,3).*SP3dx(:,:,i)+Fdx(:,:,i));
                Stk_y(:,:,i)=1/2.*(SL2(:,:,i)-Xloc_i(:,1).*SP1dy(:,:,i)-Xloc_i(:,2).*SP2dy(:,:,i)-Xloc_i(:,3).*SP3dy(:,:,i)+Fdy(:,:,i));
                Stk_z(:,:,i)=1/2.*(SL3(:,:,i)-Xloc_i(:,1).*SP1dz(:,:,i)-Xloc_i(:,2).*SP2dz(:,:,i)-Xloc_i(:,3).*SP3dz(:,:,i)+Fdz(:,:,i));
            end
        end
    else
        % return cell/array of contribution from each spheroid separately.
        if isa(Xeval,"cell")
            Stk_x=cell(1,ns); Stk_y=Stk_x; Stk_z=Stk_x;
            for i=1:ns
                Stk_x{i}=1/2.*(SL1{i}-Xeval{i}(:,1).*SP1dx{i}-Xeval{i}(:,2).*SP2dx{i}-Xeval{i}(:,3).*SP3dx{i}+Fdx{i});
                Stk_y{i}=1/2.*(SL2{i}-Xeval{i}(:,1).*SP1dy{i}-Xeval{i}(:,2).*SP2dy{i}-Xeval{i}(:,3).*SP3dy{i}+Fdy{i});
                Stk_z{i}=1/2.*(SL3{i}-Xeval{i}(:,1).*SP1dz{i}-Xeval{i}(:,2).*SP2dz{i}-Xeval{i}(:,3).*SP3dz{i}+Fdz{i});
            end
        else
            Stk_x=zeros(size(Xeval,1),size(sigma_x,2),ns); Stk_y=Stk_x; Stk_z=Stk_x;
            for i=1:ns
                Stk_x(:,:,i)=1/2.*(SL1(:,:,i)-Xeval(:,1,i).*SP1dx(:,:,i)-Xeval(:,2,i).*SP2dx(:,:,i)-Xeval(:,3,i).*SP3dx(:,:,i)+Fdx(:,:,i));
                Stk_y(:,:,i)=1/2.*(SL2(:,:,i)-Xeval(:,1,i).*SP1dy(:,:,i)-Xeval(:,2,i).*SP2dy(:,:,i)-Xeval(:,3,i).*SP3dy(:,:,i)+Fdy(:,:,i));
                Stk_z(:,:,i)=1/2.*(SL3(:,:,i)-Xeval(:,1,i).*SP1dz(:,:,i)-Xeval(:,2,i).*SP2dz(:,:,i)-Xeval(:,3,i).*SP3dz(:,:,i)+Fdz(:,:,i));
            end
        end
    end
end

function LOCAL_test_L2Stk()
    % Test Stokes operator made from Laplace operators on a system of 3
    % prolate spheroids, compared to Kernel_Eval for reference.

    clear;
    pstart=2; pend=12;
    parr=(pstart:2:pend)';
    
    p_test=16; np_test=2*p_test*(p_test+1);
    
    % % Set up --- 1 spheroid
    ns_test=1;
    pars_test=SpheroidalParameters;
    pars_test.u0=1.1; pars_test.a=1/1.1;
    pars_test.p=p_test; pars_test.oblate=0;
    pars_test.centers = [0 0 0];
    pars_test.Rmat=eye(3);

    % Set up ---- 3 spheroids
    % ns_test=3;
    % u0=[1.1 1.2 1.3];
    % pars_test=SpheroidalParameters;
    % pars_test.matvec_eta = 2;
    % pars_test.isReal=0;
    % pars_test.u0=u0;
    % alist = 1./u0;
    % pars_test.a=alist;
    % pars_test.oblate=[0 0 0];

    % pars_test.centers = [0 0 0; 150 60 0; 130 -140 150];
    % thetas = [0 pi/10 5*pi/3];
    % phis = [0 0 pi/5];

    % Ri=zeros(3,3,3);
    % for ii=1:3
    %     thetai=thetas(ii);
    %     phii=phis(ii);
    %     Riy=[cos(thetai) 0 sin(thetai); 0 1 0; -sin(thetai) 0 cos(thetai)];
    %     Riz=[cos(phii) -sin(phii) 0; sin(phii) cos(phii) 0; 0 0 1];
    %     Ri(:,:,ii)=Riz*Riy;
    % end
    % pars_test.Rmat=Ri;

    err_x=zeros(1,length(parr));
    err_y=err_x; err_z=err_x;
    for pind=1:length(parr)
        %p_test=parr(pind); 
        p_test=8;
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
    
        % This is moved into main loop below.
        % src_x=Xsrc(:,1); src_y=Xsrc(:,2); src_z=Xsrc(:,3);
        % for jj=1:ns_test
        %     ydotsig(:,:,jj)=sigma1(:,:,jj).*src_x((jj-1)*np_test+1:jj*np_test)+...
        %         sigma2(:,:,jj).*src_y((jj-1)*np_test+1:jj*np_test)+...
        %         sigma3(:,:,jj).*src_z((jj-1)*np_test+1:jj*np_test);
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
        nt=size(Xtrg_test,1);
    
        SLx=zeros(nt,1);SLy=SLx;SLz=SLx;
        
        % tic
        for ii=1:ns_test
            fprintf("\n Spheroid number %d\n", ii);
            Xself_centered=prolate_spheroid_shape(p_test,pars_test.u0(ii),pars_test.a(ii));
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
        [LSLx,LSLy,LSLz] = L2Stk(squeeze(num2cell(Xtrg_test, [1, 2]))', pars_test, sigma1, sigma2, sigma3, ns_test);
        LSLx = real(sum(LSLx,3));
        LSLy = real(sum(LSLy,3));
        LSLz = real(sum(LSLz,3));
        % toc
    
        % display(norm(sqrt(sum((SLx-LSLx).^2)./size(SLx,1))))
        % display(norm(sqrt(sum((SLy-LSLy).^2)./size(SLy,1))))
        % display(norm(sqrt(sum((SLz-LSLz).^2)./size(SLz,1))))
        err_x(pind)=norm(sqrt(sum((SLx-LSLx).^2)./size(SLx,1)));
        err_y(pind)=norm(sqrt(sum((SLy-LSLy).^2)./size(SLy,1)));
        err_z(pind)=norm(sqrt(sum((SLz-LSLz).^2)./size(SLz,1)));
    end

    figure;
    semilogy([parr,parr,parr],[err_x',err_y',err_z']);
    legend("dim 1", "dim 2", "dim 3");
    xlabel("p");
    ylabel("scaled 2-norm")
    title("L2Stk (all spectral) error from Kernel\_Eval");
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


% Test scripts for basic DY testings. By L. Crowder

close all
% ---------------------------------------------------------------------
p=8;
np=2*p*(p+1);
ns=3;

[D,trueD,DM,DMsigma,Dsigma] = testYbasic(p,10);

[DY,DY2,DYsh,DYsh2,DYsh_err]=test_DY(p,8,10,3);

[D,trueD] = testY(16,10);
% PASS. With 3 particles having nonzero sigma, obtained double layer
% matches with explicit construction of solution using Ynms

[D,Dmatvec] = compare(20,10); 
% PASS. With only one nonzero sigma,
% evaluation at points agrees with doing one DL on that nonzero sigma and
% ignoring the others.

[D,Dmatvec] = compare_self(20,10); 
% PASS. With only one nonzero sigma,
% the DL on its surface is only the self evaluation since all other
% particles make no contribution.


function pcp = PtChargePotential(ptch,Xptch,Y)
M=length(ptch);
% np=length(Y);
np=size(Y,1);
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
        EdotN(:,i)=dot(NrY,r./Rptch(:,i),2);
    end
    flux=EdotN./(4*pi*Rptch.^2)*ptch;
end


function [DY,DY2,DYsh,DYsh2,DYsh_err]=test_DY(p,ptest,eta,ns)
    % Set up spheroid system
    np=2*p*(p+1);
    
    pars=SpheroidalParameters;
    pars.matvec_eta = eta;
    if ns==3
%         pars.u0=[1.1 1.2 1.3];
        pars.u0=[1.1 1.2 1.3];
        pars.a = 1./pars.u0;
        pars.centers = [0 0 0; 4 0 0; 4 4 4];
        pars.thetas = [0 pi/10 5*pi/3];
        pars.phis = [0 0 pi/5];
    elseif ns==1
        pars.u0=1.1;
        pars.a=1/pars.u0;
        pars.centers=zeros(1,3);
        pars.phis=0;
        pars.thetas=0;
    end

    % Construct DL/SL on-surface matrices
    DM = spheroidalMatVecKernel(pars,'DL',p);

    sp=(ptest+1)^2;
    ii = (1:sp)'; 
    nn = floor(sqrt(ii-1));
    mm=ii-nn.^2-nn-1;
    Y=zeros(ns*np,sp);
    for i=1:sp
        Y(:,i) = [Yp(nn(i),mm(i),p); zeros((ns-1)*np,1)]; 
    end
    
    Ymv = permute(reshape(permute(Y,[1 3 2]),np,ns,sp),[1 3 2]);
    pars.sigma = Ymv;

    DY=DM*Y;
    
    DYsh = reshape(shAna(reshape(DY,np,[])),[],sp);
    figure; 
    imagesc(log10(abs((DYsh)))); 
    colorbar;
    title('log10 of matrix product coefficients');
    clim([-16,0]);
    
    
    DY2 = spheroidalMatVec(pars,'DL');
    DY2 = reshape(permute(DY2,[1 3 2]),ns*np,[]);
    DYsh2 = reshape(shAna(reshape(DY2,np,[])),[],sp);

    figure; 
    imagesc(log10(abs((DYsh2)))); 
    colorbar;
    title('log10 of spheroidalMatVec coefficients');
    clim([-16,0]);

    DYsh_err = abs(DYsh-DYsh2);

    figure; 
    imagesc(log10(abs(DYsh_err))); 
    colorbar;
    title('log10 error of coefficients');

end

function [D,Dmatvec]=compare(p,eta)
    
    np=2*p*(p+1);
    pars=SpheroidalParameters;
    pars.matvec_eta = eta;

    pars.u0=[1.1 1.2 1.3];
    pars.a = 1./pars.u0;
    pars.centers = [0 0 0; 2 0 0; 2 2 2];
    pars.thetas = [0 pi/10 5*pi/3];
    pars.phis = [0 0 pi/5];
    pars.sigma = reshape([zeros(np,1) Yp(2,2,p) zeros(np,1)],[],1,3);
    
    target_distances=[1.2 1.2 2];
    % Target points to test at
    Y=pars.get_X;
    td = repmat(target_distances,np,1);
    td = td(:);
    cen = reshape(repmat(pars.centers',np,1),3,[])';
    displacement = Y-cen;
    Xeval = td.*displacement +cen;

    Dmatvec = spheroidalMatVec(pars,'DL',Xeval);
    Xt = pars.get_X_targets(Xeval);
    Dcell = spheroidalDL(pars,{[],Xt(:,:,2),[]});
    D=Dcell{2};

    fprintf('p=%d: infinity-norm err=%.6e\n',p, max(abs(D-Dmatvec)))
    
    figure; 
    hold on;
    plot(1:length(D),D);
    plot(1:length(Dmatvec),Dmatvec,'--');
    legend('spheroidalDL','spheroidalMatVec')
    hold off;
    
    figure;
    semilogy(1:length(D),abs(D-Dmatvec)/abs(D))
    xlabel('evaluation point');
    ylabel('relative error')
end

function [D,Dmatvec]=compare_self(p,eta)
    
    np=2*p*(p+1);
    pars=SpheroidalParameters;
    pars.matvec_eta = eta;
    pars.isReal=1;

    pars.u0=[1.1 1.1 1.1];
    pars.a = 1./pars.u0;
    pars.centers = [0 0 0; 4 0 0; 4 4 4];
    pars.thetas = [0 pi/10 5*pi/3];
    pars.phis = [0 0 pi/5];
    pars.sigma = reshape([zeros(np,1) Yp(2,2,p) zeros(np,1)],[],1,3);
    
    Dmatvec = spheroidalMatVec(pars,'DL');
    D = spheroidalDL(pars);

    D=D(:,1,2);
    Dmatvec=Dmatvec(:,1,2);

    fprintf('p=%d: infinity-norm err=%.6e\n',p, max(abs(D-Dmatvec)))
    
    figure; 
    hold on;
    plot(1:length(D),D);
    plot(1:length(Dmatvec),Dmatvec,'--');
    legend('spheroidalDL','spheroidalMatVec')
    hold off;
    
    figure;
    semilogy(1:length(D),abs(D-Dmatvec)/abs(D))
    xlabel('evaluation point');
    ylabel('relative error')
end

function [D,trueD] = testY(p,eta)
    np=2*p*(p+1);
    ns=3;

    pars=SpheroidalParameters;
    pars.matvec_eta = eta;
    pars.u0=[1.1 1.2 1.3];
    pars.a = 1./pars.u0;
    pars.centers = [0 0 0; 2 0 0; 2 2 2];
    pars.thetas = [0 pi/10 5*pi/3];
    pars.phis = [0 0 pi/5];

    target_distances=[1.1 1.1 2];

    pars.sigma=reshape([Yp(2,0,p) Yp(2,1,p) Yp(2,2,p)],[],1,ns);
    pars.sigma = [pars.sigma pars.sigma];
    [~,~,le]=DLspectrum(p,pars.u0);

    %relevant eigenvalues for the Ynms used above (Y_2^0, Y_2^1, Y_2^2)

    le20 = le(7,1);
    le21 = le(8,2);
    le22 = le(9,3);

    Y=pars.get_X;

    % Target points to test at
    td = repmat(target_distances,np,1);
    td = td(:);
    cen = reshape(repmat(pars.centers',np,1),3,[])';
    displacement = Y-cen;
    Xeval = td.*displacement + cen;

    pars.plot(Xeval);

    sep=pars.separate_targets(Xeval);
    sep = sum(sep,2) ~=0;
    nfar = sum(sep);
    nclose = sum(~sep);
    fprintf('%d targets are near a spheroid\n',nclose);
    fprintf('%d targets are far from a spheroid\n',nfar);

    Xt = pars.get_X_targets(Xeval);
    St = cart2spheroidal(Xt,pars.a);
    ut = permute(St(:,1,:),[1 3 2]);
    vt = permute(St(:,2,:),[1 3 2]);
    phit = permute(St(:,3,:),[1 3 2]);

    Y20 = Ynm(2,0,acos(vt(:,1))',phit(:,1));
    Y21 = Ynm(2,1,acos(vt(:,2))',phit(:,2));
    Y22 = Ynm(2,2,acos(vt(:,3))',phit(:,3));

    PQ=legendre_otc(2,ut(:),1);
    Q=PQ{2};
    F20=reshape(Q(7,:)',np*ns,[]);
    F20 = F20(:,1);
    F21=reshape(Q(8,:)',np*ns,[]);
    F21 = F21(:,2);
    F22=reshape(Q(9,:)',np*ns,[]);
    F22 = F22(:,3);

    trueD = le20.*F20.*Y20 + le21.*F21.*Y21 + le22.*F22.*Y22; 
    disp(size(trueD))

    D = spheroidalMatVec(pars,'DL',Xeval);

    err1 = abs(D(:,1)-trueD)./abs(trueD);
    err2 = abs(D(:,2)-trueD)./abs(trueD);
    fprintf('p=%d: infinity-norm err1=%.6e\n',p, max(err1))
    fprintf('p=%d: infinity-norm err2=%.6e\n',p, max(err2))

    figure; 
    hold on;
    plot(1:length(trueD),real(trueD));
    plot(1:length(D(:,1)),real(D(:,1)),'--');
    legend('True D[Y]','spheroidalMatVec')
    title('Real')
    hold off;

    figure; 
    hold on;
    plot(1:length(trueD),imag(trueD));
    plot(1:length(D(:,1)),imag(D(:,1)),'--');
    legend('True D[Y]','spheroidalMatVec')
    hold off;
    title('Imaginary')

    figure;
    semilogy(1:length(err1),err1);
    title('error 1');

    figure;
    semilogy(1:length(err2),err2);
    title('error 2');

end

function [D,trueD,DM,DMsigma,Dsigma] = testYbasic(p,eta)
    np=2*p*(p+1);
    ns=2;

    pars=SpheroidalParameters;
    pars.matvec_eta = eta;
    pars.u0=[1.1 1.1];
    pars.a = 1./pars.u0;
    pars.centers = [0 0 0; 2 0 0];
    pars.thetas = [0 0];
    pars.phis = [0 pi/2];
    pars.sigma=[repmat(Yp(2,1,p),1,1,ns) repmat(Yp(2,2,p),1,1,ns)];

    Xeval = [1 1 1];

    pars.plot(Xeval);

    
    [~,~,le]=DLspectrum(2,pars.u0);

    %relevant eigenvalues for the Ynm used above Y_2^1
    le21 = le(8,:)';
    le22 = le(9,:)';
    disp(le21);


    % SOMETHING MIGHT BE FUNNY WITH THE DIMENSIONS HERE
    Xt = pars.get_X_targets(Xeval);
    St = cart2spheroidal(Xt,pars.a);
    ut = permute(St(:,1,:),[1 3 2]);
    vt = permute(St(:,2,:),[1 3 2]);
    phit = permute(St(:,3,:),[1 3 2]);

    disp(ut)
    disp(acos(vt))
    disp(phit)

    Y21 = Ynm(2,1,acos(vt),phit');
    disp(Y21)

    Y22 = Ynm(2,2,acos(vt),phit');

    PQ=legendre_otc(2,ut(:),1);
    Q=PQ{2};
    F21=Q(8,:)';
    F22=Q(9,:)';
    disp(F21)

    Dvec = [le21.*F21.*Y21 le22.*F22.*Y22];
    disp(Dvec)
    trueD = sum(Dvec,1); 

    D = spheroidalMatVec(pars,'DL',Xeval);

    D1 = spheroidalDL(pars,{Xeval,Xeval-pars.centers(2,:)});
    disp(D1)
    disp(sum(cell2mat(D1)))

    err = abs(D-trueD)./abs(trueD);
    for i=1:length(err)
        fprintf('p=%d: err(%d)=%.6e\n',p, i,err(i));
    end


%     Test DM matrix
    DM = spheroidalMatVecKernel(pars,'DL',p);
    sigma_vec = reshape(permute(pars.sigma,[1 3 2]),np*ns,[],1);
    DMsigma = DM*sigma_vec;
    Dsigma = spheroidalMatVec(pars,'DL');
    Dsigma = reshape(permute(Dsigma,[1 3 2]),np*ns,[],1);
    err2 = abs(DMsigma-Dsigma)./abs(Dsigma);

    for i=1:size(err2,2)
        fprintf('p=%d: infinity-norm matrix product err(%d)=%.6e\n',p,i, max(err2(:,i)))
    end

    figure; 
    hold on;
    plot(1:length(DMsigma(:,1)),real(DMsigma(:,1)),'+-');
    plot(1:length(Dsigma(:,1)),real(Dsigma(:,1)));
    legend('Matrix product','spheroidalMatVec')
    title('Real')
    hold off;

    figure; 
    hold on;
    plot(1:length(DMsigma(:,1)),imag(DMsigma(:,1)),'+-');
    plot(1:length(Dsigma(:,1)),imag(Dsigma(:,1)));
    legend('Matrix product','spheroidalMatVec')
    title('Imaginary')
    hold off;

    figure; 
    semilogy(1:length(DMsigma(:,1)),abs(imag((DMsigma(:,1)-Dsigma(:,1)))));
    title('Imaginary error')

    figure;
    semilogy(1:length(err2),err2);
end

function [soln,soln_quad,truesoln,sigma_vec] = test_three(p,eta,u0,target_distances,useS,plt,obl,neumann,test_quad)
    % (1) set up system of 3 spheroids, 2 close and 1 far
    % (2) put a few point charges around the center of each spheroid.
    %       - Calculate the potential explicitly.
    %       - Determine potential on the boundary of each spheroid for BCs, f
    % (3) Use MatVec self-evaluation (no X) to construct matrix operator D.
    % (4) Solve (1/2*I + D + Completion)*sigma = f
    %       * Completion term will need some care, since it's multiple particles.
    %       For each particle it will be something like sum_i[ 1/|| x- c_i||
    %       integral(sigma_i) ]
    % (5) Compare sigma with true point charge potential at targets that
    %     are 'target_distances' away from surface of each spheroid.

    % Set up spheroid system
    np=2*p*(p+1);
    ns=3;
    
    pars=SpheroidalParameters;
    pars.matvec_eta = eta;
    pars.isReal=0;
    pars.u0=u0;

    if ~obl
        pars.a = 1./u0;
        pars.oblate = [0 0 0];
    else
        % pars.a = 1./sqrt(u0.^2+1);
        % pars.oblate = [1 1 1];
        pars.oblate=[1,0,1];
        pars.a=~pars.oblate./pars.u0 + pars.oblate./sqrt(pars.u0.^2+1);
    end
    pars.centers = [0 0 0; 5 0 0; 3.2 3.2 3.2];
    % pars.centers=[0 0 0;3 0 0;1.6 1.6 1.6];
    % pars.centers=[-0.5 -0.5 0;1.5 0 0;0.8 0.8 0.8];

    thetas = [0 pi/10 5*pi/3];
    phis = [0 0 pi/5];
    Ri=zeros(3,3,3);
    for i=1:3
        thetai=thetas(i);
        phii=phis(i);
        Riy=[cos(thetai) 0 sin(thetai); 0 1 0; -sin(thetai) 0 cos(thetai)];
        Riz=[cos(phii) -sin(phii) 0; sin(phii) cos(phii) 0; 0 0 1];
        Ri(:,:,i)=Riz*Riy;
    end
    pars.Rmat=Ri;
    pars.thetas=thetas;
    pars.phis=phis;
    
    pars.sigma=zeros(np,1,ns);
        
    % Point charges
    rng(13.5);
    nc = 2;
    scale = repmat(.5.*pars.a .*sqrt(pars.u0.^2-1),nc,1);
    scale = scale(:);
    c=pars.centers;
    Xptch = reshape(repmat(reshape(c',3,1,[]),1,nc),3,[],1)';
    % random point charges
    d = repmat(scale,1,3).*(rand(size(Xptch))-.5);
    Xptch = Xptch + d; %point charge locations, near centers of spheroids
    ptch = (2.*rand(ns*nc,1)-1); %charge value, between -1 and 1
    
    % if plt
    %     pars.plot(Xptch);
    %     title("spheroids and point charges")
    % end

    Y=pars.get_X;
    [NrY,~]=pars.get_Norm_rot(p);
    
    [~, gwt]=g_grid(p+1);
    wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
    wt = wt(:)';
    
    if ~neumann
        if useS
            CM = spheroidalMatVecKernel(pars,'SL',p);
        else
            CM=[];
            for i=1:ns
                if ~pars.oblate(i)
                    Yi = prolate_spheroid_shape(p,pars.u0(i),pars.a(i),'cart');
                else
                    Yi = oblate_spheroid_shape(p,pars.u0(i),pars.a(i),'cart');
                end
                ci=c(i,:);
                
                % Completion matrix operator: sum of 1/||x-c_i|| * integral of each sigma over
                % respective spheroid surface
                Sns = SurfaceSph(Yi(:));
                W = Sns.geoProp.W;
            
                RY = vecnorm(Y-repmat(ci,np*ns,1),2,2);
                CM = [CM (1./RY)*(W'.*wt)];
            end
        end
        % We will use this for surface boundary conditions
        truesolnSurf=PtChargePotential(ptch,Xptch,Y);

        % Construct DL on-surface matrices
        DM = spheroidalMatVecKernel(pars,'DL',p);
        
        % Final operator: 1/2 I + DM + CM
        K = .5*eye(ns*np) + DM +CM;
    else
        truesolnSurf=PtChargeFlux(ptch, Xptch,Y,NrY);

        % Construct SP on-surface matrices
        SPM= spheroidalMatVecKernel(pars,'SP',p);
        
        % Final operator: -1/2 I + SPM
        K = -.5*eye(ns*np) + SPM;
    end

    if useS
        fprintf("\n Condition number for p=%d of K=-0.5I+D+S is %f \n ", p,cond(K));
    else
        fprintf("\n Condition number for p=%d of K=-0.5I+D+C is %f \n", p,cond(K));
    end
    
    sigma_vec=K\truesolnSurf;
    sigma = reshape(sigma_vec,np,[],ns);
    pars.sigma = sigma;
    % 
    % addpath('./spheroidal/plot/utils/');
    % cmp = getPyPlot_cMap('rainbow', [], [], '"/opt/homebrew/bin/python3"');

    % if p==16
    %     figure;
    %     for sphind=1:3
    %         Yi=Y((sphind-1)*np+1:sphind*np,:);
    %         shades=pars.sigma(:,:,sphind);
    %         if norm(real(shades)-shades)>1e-8
    %             Error("density imaginary.");
    %         end
    %         plotb([Yi(:,1);Yi(:,2);Yi(:,3)],real(shades)); hold on;
    %     end
    %     colormap(cmp);
    %     colorbar; hold off;
    % end
    
    % Target points to test at
    [theta,phi]=gl_grid(p);
    v=cos(theta);
    Xcell=cell(1,ns);
    for i=1:ns
        u0=pars.u0(i);
        if ~pars.oblate(i)
            normal=1./sqrt((u0^2-1).*(u0^2-v.^2)).*[u0.*sqrt(u0.^2-1).*sqrt(1-v.^2).*cos(phi) u0.*sqrt(u0.^2-1).*sqrt(1-v.^2).*sin(phi) (u0.^2-1).*v];
            Yi = prolate_spheroid_shape(p,pars.u0(i),pars.a(i),'cart');
        else
            normal=1./sqrt((u0^2+1).*(u0^2+v.^2)).*[u0.*sqrt(u0.^2+1).*sqrt(1-v.^2).*cos(phi) u0.*sqrt(u0.^2+1).*sqrt(1-v.^2).*sin(phi) (u0.^2+1).*v];
            Yi = oblate_spheroid_shape(p,pars.u0(i),pars.a(i),'cart');
        end
        Xcell{i} = Yi + target_distances(i).*normal;
    end
    Xeval = pars.set_X_targets(Xcell);

    % if plt
    %     pars.plot(Xeval);
    %     title('spheroids and target points')
    % end
    
    if ~neumann
        if useS
            C = spheroidalMatVec(pars,'SL',Xeval);
        else
            C=zeros(size(Xeval,1),1);
            for i=1:ns
            
                if ~pars.oblate(i)
                    Yi = prolate_spheroid_shape(p,pars.u0(i),pars.a(i),'cart');
                else
                    Yi = oblate_spheroid_shape(p,pars.u0(i),pars.a(i),'cart');
                end
                ci=c(i,:);
                
                % Completion matrix operator: sum of 1/||x-c_i|| * integral of each sigma over
                % respective spheroid surface
                Sns = SurfaceSph(Yi(:));
                sigmaSurfInt = integrateOverS(Sns,sigma(:,:,i));
            
                RXeval = vecnorm(Xeval-repmat(ci,size(Xeval,1),1),2,2);
                
                C = C+ sigmaSurfInt./RXeval;
            end
        end

        soln=spheroidalMatVec(pars,'DL',Xeval) + C;

    else

        soln=spheroidalMatVec(pars,'SL',Xeval);

    end

    truesoln=PtChargePotential(ptch,Xptch,Xeval); 

    fprintf('p=%d: spectral method err=%.6e\n',p, max(abs(truesoln-soln)))  

    if test_quad
        % Quadrature resolution near singularity
        soln_quad=zeros(size(Xeval,1),1);
        for i=1:ns
            if ~pars.oblate(i)
                Yi = prolate_spheroid_shape(p,pars.u0(i),pars.a(i),'cart');
            else
                Yi = oblate_spheroid_shape(p,pars.u0(i),pars.a(i),'cart');
            end
            if neumann
                KEparams = Kernel_Eval_parameters('SL_L_3D',0,1,1,1,1e-8,2,400,1);
            else
                KEparams = Kernel_Eval_parameters('DL_L_3D',0,1,1,1,1e-8,2,400,1);
            end
            KEparams.dim = 3;
            % KEparams.X = Yi; 
            KEparams.X = Y((i-1)*np+1:i*np,:);
            Sns = SurfaceSph(Yi);
            Wns = Sns.geoProp.W; Wns= Wns.*wt';
            KEparams.W2 = Wns.';
            % Nrns = reshape(Sns.geoProp.nor.to_array,[],3);
            % KEparams.nor = Nrns;
            % KEval_mat=Kernel_Eval(Xeval,Yi,KEparams); 
            [Nrns,~] = pars.get_Norm_rot(p);
            KEparams.nor = Nrns((i-1)*np+1:i*np,:); 
            KEval_mat=Kernel_Eval(Xeval,Y((i-1)*np+1:i*np,:),KEparams); 
            soln_quad=soln_quad + KEval_mat*pars.sigma(:,:,i);
        end
        if ~neumann
            soln_quad=soln_quad+C;
        end
        fprintf('p=%d: quadrature err=%.6e\n',p, max(abs(truesoln-soln_quad))) 
    else
        soln_quad=zeros(size(Xeval,1),1);
    end

end

function test_three_plot(p,eta,u0,plt)
    

    % Set up spheroid system
    np=2*p*(p+1);
    ns=3;
    
    pars=SpheroidalParameters;
    pars.matvec_eta = eta;
    pars.isReal=0;
    pars.u0=u0;

    pars.a = 1./u0;
    pars.oblate = [0 0 0];
    
    pars.centers = [0 0 0; 5 0 0; 3.2 3.2 3.2];
    thetas = [0 pi/10 5*pi/3];
    phis = [0 0 pi/5];
    Ri=zeros(3,3,3);
    for i=1:3
        thetai=thetas(i);
        phii=phis(i);
        Riy=[cos(thetai) 0 sin(thetai); 0 1 0; -sin(thetai) 0 cos(thetai)];
        Riz=[cos(phii) -sin(phii) 0; sin(phii) cos(phii) 0; 0 0 1];
        Ri(:,:,i)=Riz*Riy;
    end
    pars.Rmat=Ri;
    pars.thetas=thetas;
    pars.phis=phis;
    
    pars.sigma=zeros(np,1,ns);
        
    % Point charges
    nc = 2;
    scale = repmat(.5.*pars.a .*sqrt(pars.u0.^2-1),nc,1);
    scale = scale(:);
    c=pars.centers;
    Xptch = reshape(repmat(reshape(c',3,1,[]),1,nc),3,[],1)';
    % random point charges
    d = repmat(scale,1,3).*(rand(size(Xptch))-.5);
    Xptch = Xptch + d; %point charge locations, near centers of spheroids
    ptch = (2.*rand(ns*nc,1)-1); %charge value, between -1 and 1

    if ~plt
        return;
    end

    Y=pars.get_X;

    % truesolnSurf=PtChargePotential(ptch,Xptch,Y);

    Y=pars.get_X;
    [NrY,~]=pars.get_Norm_rot(p);
    
    [~, gwt]=g_grid(p+1);
    wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
    wt = wt(:)';

    CM=[];
    for i=1:ns
        if ~pars.oblate(i)
            Yi = prolate_spheroid_shape(p,pars.u0(i),pars.a(i),'cart');
        else
            Yi = oblate_spheroid_shape(p,pars.u0(i),pars.a(i),'cart');
        end
        ci=c(i,:);
        
        % Completion matrix operator: sum of 1/||x-c_i|| * integral of each sigma over
        % respective spheroid surface
        Sns = SurfaceSph(Yi(:));
        W = Sns.geoProp.W;
    
        RY = vecnorm(Y-repmat(ci,np*ns,1),2,2);
        CM = [CM (1./RY)*(W'.*wt)];
    end

    % We will use this for surface boundary conditions
    truesolnSurf=PtChargePotential(ptch,Xptch,Y);

    % Construct DL on-surface matrices
    DM = spheroidalMatVecKernel(pars,'DL',p);
    
    % Final operator: 1/2 I + DM + CM
    K = .5*eye(ns*np) + DM +CM;
    
    sigma_vec=K\truesolnSurf;
    sigma = reshape(sigma_vec,np,[],ns);
    pars.sigma = sigma;



    addpath('./plot/utils/');
    cmp = getPyPlot_cMap('rainbow', [], [], '"/opt/homebrew/bin/python3"');

    for sphind=1:3
        Yi=Y((sphind-1)*np+1:sphind*np,:);
        % shades=truesolnSurf((sphind-1)*np+1:sphind*np,:);
        shades=pars.sigma(:,:,sphind);
        plotb([Yi(:,1);Yi(:,2);Yi(:,3)],real(shades)); hold on;
        colormap(cmp)
    end
    view(45,10);
    cb=colorbar; 
    set(cb,'Position',[0.3 0.25 0.01 0.3]);
    set(cb,'YTick',-0.1:0.1:0.1);
    hold off;
end

function [errs_spectral,errs_quad] = compare_quad(useS,if_neumann)
    % Compare exponential growth of error for smooth quadrature as distance
    % decrease, to maintained accuracy for spectral methods.
    % Set up: 3 prolate spheroids in space, with 2 random point charges
    % inside each. Use Dirichlet or Neumann problem formulation.

    parr = [4,8,12,16,24];  %p-values to test
    distance=1.*10.^(-(1:6)); %Distance from spheroid surface
    errs_spectral = cell(1,length(parr));
    errs_quad = cell(1,length(parr));
    for l=1:length(parr)
        p=parr(l);
        np=2*p*(p+1);
        this_p_errs_spec=zeros(np*3,length(distance));
        this_p_errs_quad=zeros(np*3,length(distance));
        for j=1:length(distance)
            [soln,soln_quad,truesoln,~] = test_three(p,10,[1.1 1.2 1.3],distance(j).*ones(1,3),useS,0,1,if_neumann,1);
            this_p_errs_spec(:,j) = log10(abs((soln-truesoln)./max(truesoln)));
            this_p_errs_quad(:,j) = log10(abs((soln_quad-truesoln)./max(truesoln)));
        end
        errs_spectral{l}=this_p_errs_spec;
        errs_quad{l}=this_p_errs_quad;
    end
    
    figure;
    for l=1:length(parr)
        subplot(1,length(parr),l)
        hold on;
        ylim([-16,10])
        boxchart(errs_spectral{l},BoxFaceColor='b');
        boxchart(errs_quad{l},BoxFaceColor='g');
        xlabel("-log10(distance)")
        if l==1
            ylabel("log10 relative error")
        end
        title("p="+string(parr(l)));
        hold off;
    end

end

% function [u0,a,oblate,centers,Rmats,ns_actual]=setup_large(ns,p)
%     % ns large: create random particles and make sure they don't
%     % overlap.
%       % NOTE: replaced by error-catching code for setup.
% 
%     % oblate=randi([0 1],1,ns);
%     % u0=rand(1,ns)+1;
% 
%     [sys_params,ns_actual] = call_randspheroids(ns,p);
%     centers=zeros(ns_actual,3);
%     Rmats=zeros(3,3,ns_actual);
%     oblate=zeros(1,ns_actual);
%     for ii2=1:ns_actual
%         centers(ii2,:)=sys_params(ii2).C;
%         Rmats(:,:,ii2)=sys_params(ii2).R;
%         oblate(ii2)=strcmp(sys_params(ii2).shape,'oblateZ');
%     end
%     u0=2/sqrt(3).*ones(size(oblate));
%     a=~oblate./u0 + oblate./sqrt(u0.^2+1);
% end





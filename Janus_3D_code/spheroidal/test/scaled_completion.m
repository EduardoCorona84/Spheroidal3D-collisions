%% minimizing cond(K)

% % Extreme epsilon -- slender body 
% eps=10.^(-1:-1:-4);

% Regular to more elongated spheroids
R = 2.^linspace(0.1,6.01,30);
% R=[1.1,1.5,2,4,8,20,50];
eps = 1./R;

p=8;

scale_list=zeros(size(eps));
cond_list = zeros(size(eps));

for eind=1:length(eps)
    
    ep=eps(eind);
    fprintf("\n epsilon=%d",ep)

    [x, condval] = fmincon(@(s) minfunc(ep,s,p,0),1,-1,-1e-10);

    scale_list(eind) = x;
    cond_list(eind) = condval;
end

% E = sqrt(1-eps.^2);
% etaC  = 1./( 2*pi.*eps.*(1+(1/2./ep./E).*asin(E)) );

% % SBT scaling of SL completion
% etaS = 1./2./eps./log(1./eps);
% etaS = 1./2./eps./log(1./2./eps);

% for moderate eps: attempts to fit minimizer result.
% etaS = 0.4.*eps.^(-0.7);
% etaS = 1./2./eps./log(2./eps);
etaS1 = 1./1.5./eps./log(3./eps);
etaS2 = 1./1.5./eps./log(4./eps);

figure(1);
loglog(eps,scale_list,"-*"); hold on;
loglog(eps,etaS1); 
loglog(eps,etaS2); hold off;
xlabel("$\epsilon$","interpreter","latex")
ylabel("$\eta$","interpreter","latex")
legend("minimizing cond(K)","$\frac{1}{1.5\epsilon\log(3/\epsilon)}$","$\frac{1}{1.5\epsilon\log(4/\epsilon)}$","interpreter","latex")

figure(2);
loglog(eps,cond_list,"-*"); hold on;
xlabel("$\epsilon$","interpreter","latex")
ylabel("cond(K)","interpreter","latex")

%% spectra of S, C, and I/2+D

% % Extreme epsilon -- slender body 
% eps=10.^(-1:-1:-3);

% Regular to more elongated spheroids
R=[1.1,1.5,2,4,8,20,50];
eps = 1./R;

p=16; np=2*p*(p+1);

conds = zeros(4,length(eps));

for eind=1:length(eps)

    ep=eps(eind);

    [CM,SM,DM,u0] = get_mats(ep,p); ns=1;
    % [CM,SM,DM,u0] = get_mats_two(ep,p); ns=2;

    % Analytic lambda_n^0 since Ynm is eigenfunction to D on surface.
    % Only works when one spheroid
    ind_n0 = zeros(1,p+1);
    for n=0:p
        if n==0
            ind_n0(n+1) = 1;
        else
            ind_n0(n+1) = ind_n0(n)+2*n;
        end
    end
    L=legendre_otc(p,u0,1,1,1);
    P=L{1}; Q=L{2}; dP=L{3}; dQ=L{4};
    lambda_surf=(u0^2-1)/2.*(P(ind_n0).*dQ(ind_n0)+dP(ind_n0).*Q(ind_n0));
    lambda_surf = lambda_surf + 0.5;

    k=0:p;

    % Get eigenvalues of matrices
    lambdaD = eig(.5*eye(np*ns) + DM);
    lambdaD = sort(lambdaD);

    lambdaS = eig(SM);
    lambdaS = sort(lambdaS, 'descend');

    lambdaC = eig(CM);
    lambdaC = sort(lambdaC, 'descend');

    % unscaled completion matrix
    KS2 = .5*eye(np*ns) + DM + SM;
    lambdaKS2 = eig(KS2);

    KC2 = .5*eye(np*ns) + DM + CM;
    lambdaKC2 = eig(KC2);
    
    % ------- using scalar from eigenvalue fitting -----------------------
    % % scaled completion matrix
    % etaS = 1/2/ep/log(1/2/ep);
    % KS = .5*eye(np*ns) + DM + etaS*SM;
    % lambdaKS = eig(KS);
    % 
    % E = sqrt(1-ep^2);
    % etaC  = 1/( 2*pi*(ep).*(1+(1/2/ep/E).*asin(E)) );
    % KC = .5*eye(np*ns) + DM + etaC*CM;
    % lambdaKC = eig(KC);
    % 
    % Sscale=ep.*log(1/ep./k(2:end));
    % Dscale=1/2*ep^2.*k(2:end).^2.*log(1/ep./k(2:end));
    % --------------------------------------------------------------------

    % -------- use scaling from minimizer -------------------------------
    [etaC,condKC] = fmincon(@(s) cond(.5*eye(np*ns)+DM+s*CM),1,-1,-1e-10);
    % [etaS,condKS] = fmincon(@(s) cond(.5*eye(np*ns)+DM+s*SM),1,-1,-1e-10);
    [etaS_uncon,condKS_uncon] = fminsearch(@(s) cond(.5*eye(np*ns)+DM+s*SM),1);

    KC = .5*eye(np*ns)+DM+etaC*CM;
    lambdaKC = eig(KC);
    lambdaKC = sort(lambdaKC);

    KS = .5*eye(np*ns)+DM+etaS_uncon*SM;
    lambdaKS = eig(KS);
    lambdaKS = sort(lambdaKS);
    % 
    % figure(7)
    % plot(ep,condKS,'*'); hold on;
    % plot(ep,condKS_uncon,"o");
    % 
    % figure(8)
    % plot(ep,etaS,'*'); hold on;
    % plot(ep,etaS_uncon,"o");
    % --------------------------------------------------------------------

    % store condition numbers
    conds(1,eind) = cond(KC);
    conds(2,eind) = cond(KC2);
    conds(3,eind) = cond(KS);
    conds(4,eind) = cond(KS2);

    % D eigenvalues and scaling
    figure(1)
    plot(k(:),real(lambdaD(1:p+1)),"-*"); hold on;
    % if eind==1
    %     plot(k(1:p/2+1),[0,Dscale(1:p/2)],"--");
    % else
    %     plot(k(:),[0,Dscale(1:p)],"--");
    % end
    % plot(k(:),lambda_surf(:),"-.");

    % S eigenvalues and scaling
    figure(2)
    plot(k(:),real(lambdaS(1:p+1)),"-*"); hold on;

    % KS_scaled vs unscaled eigenvalues 
    % (with different indexing: )
    %   (ind_n0) for assuming order is Ynm; 
    %   (1:p+1) for the first p+1 eigenvalues; 
    %   sort() then (1:p+1) for the first p+1 smallest.
    figure(3)
    plot(k(:),real(lambdaKS(1:p+1)),"-*"); hold on;
    % plot(k(:),real(lambdaKS2(1:p+1)),"-o"); hold on;

    % KC_scaled vs unscaled eigenvalues
    figure(4)
    plot(k(:),real(lambdaKC(1:p+1)),"-*"); hold on;
    % plot(k(:),real(lambdaKC2(1:p+1)),"-o"); hold on;

end
figure(1)
legend("$\epsilon$=0.91","","","","","","$\epsilon$=0.02","interpreter","latex");
% legend("$\epsilon$=0.1","$\epsilon\frac{1}{2}(k\epsilon)^2\log|k\epsilon|^{-1}$","analytical", ...
%     "$\epsilon$=0.01","$\epsilon\frac{1}{2}(k\epsilon)^2\log|k\epsilon|^{-1}$","analytical", ...
%     "$\epsilon$=0.001","$\epsilon\frac{1}{2}(k\epsilon)^2\log|k\epsilon|^{-1}$","analytical", ...
%     "interpreter","latex");
xlabel("n")
% title("I/2+D eigenvalues for Yn,m=0 for two spheroid system in mode 3 (- -)")
title("I/2+D eigenvalues single spheroid")
ylabel("$\lambda^D$","interpreter","latex")
hold off;

figure(2)
legend("$\epsilon$=0.91","","","","","","$\epsilon$=0.02","interpreter","latex");
% legend("$\epsilon$=0.1","$\epsilon\log(k \epsilon)^{-1}$", ...
%     "$\epsilon$=0.01","$\epsilon\log(k \epsilon)^{-1}$", ...
%     "$\epsilon$=0.001","$\epsilon\log(k \epsilon)^{-1}$", ...
%     "interpreter","latex");
% xlabel("n")
% title("largest 17 S eigenvalues for two spheroid system in mode 3 (- -)")
title("largest 17 S eigenvalues single spheroid")
ylabel("$\lambda^S$","interpreter","latex")
hold off;

figure(3)
% title("KS eigenvalues, first 17")
title("KS eigenvalues, smallest 17")
legend("$\epsilon$=0.91","","","","","","$\epsilon$=0.02","interpreter","latex");
% legend("$\epsilon$=0.91","$\epsilon$=0.91 unscaled", "","","","","","","","","","","$\epsilon$=0.02","",...
%     "interpreter","latex");
% legend("$\epsilon$=0.1","$\epsilon$=0.1, unscaled","$\epsilon$=0.01","$\epsilon$=0.01, unscaled","$\epsilon$=0.001","$\epsilon$=0.001, unscaled","interpreter","latex");
ylabel("$\lambda^{KS}$","interpreter","latex")
hold off;

figure(4)
% title("KC eigenvalues, first 17")
title("KC eigenvalues, smallest 17")
legend("$\epsilon$=0.91","","","","","","$\epsilon$=0.02","interpreter","latex");
% legend("$\epsilon$=0.91","$\epsilon$=0.91 unscaled", "","","","","","","","","","","$\epsilon$=0.02","",...
%     "interpreter","latex");
% legend("$\epsilon$=0.1","$\epsilon$=0.1, unscaled","$\epsilon$=0.01","$\epsilon$=0.01, unscaled","$\epsilon$=0.001","$\epsilon$=0.001, unscaled","interpreter","latex");
ylabel("$\lambda^{KC}$","interpreter","latex")
hold off;

% Condition number with minimized etaS, with or without constraints.
figure(5)
loglog(eps,conds(1,:)); hold on;
loglog(eps,conds(2,:)); hold off;
legend("KC scaled","KC unscaled")
xlabel("$\epsilon$","interpreter","latex")
ylabel("condition number")

figure(6)
loglog(eps,conds(3,:)); hold on;
loglog(eps,conds(4,:)); hold off;
legend("KS scaled","KS unscaled")
xlabel("$\epsilon$","interpreter","latex")
ylabel("condition number")

% figure(7)
% xlabel("$\epsilon$","interperter","latex")
% ylabel("condition number")
% legend("constrained min","unconstrained min")
% 
% figure(8)
% xlabel("$\epsilon$","interperter","latex")
% ylabel("$\eta$","interpreter","latex")
% legend("constrained min","unconstrained min")


%% Functions

% For system of two prolates in mode=3 (parallel minor axis) configuration,
% with SBT aspect ratio parameter 'eps', give condition number of final
% completed DL operator when completed using completion flow or SL.
function condK = minfunc(eps,eta,p,useS)
    np=2*p*(p+1);
    R = 1/eps;

    % two particles
    ns = 2;
    pars = set_up(3,p,1,R);

    Y = pars.get_X;
    [~, gwt]=g_grid(p+1);
    wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
    wt = wt(:)';

    if useS
        % Completion matrix operator using just SL
        CM = spheroidalMatVecKernel(pars,'SL',p);
    else
        % Completion matrix operator: sum of 1/||x-c_i|| * integral of each sigma over
        %       respective spheroid surface
        CM=[];
        for i=1:ns
            if ~pars.oblate(i)
                Yi = prolate_spheroid_shape(p,pars.u0(i),pars.a(i),'cart');
            else
                Yi = oblate_spheroid_shape(p,pars.u0(i),pars.a(i),'cart');
            end
            ci=pars.centers(i,:);
            Sns = SurfaceSph(Yi(:));
            W = Sns.geoProp.W;
            RY = vecnorm(Y-repmat(ci,np*ns,1),2,2);
            CM = [CM (1./RY)*(W'.*wt)];
        end
    end

    DM = spheroidalMatVecKernel(pars,'DL',p);
    K = .5*eye(np*ns) + DM + eta.*CM;
    condK = cond(K);

end

% Get BIO matrices for one prolate system of SBT parameter 'eps'.
function [CM,SM,DM,this_u0] = get_mats(eps,p)
    np=2*p*(p+1);
    R = 1/(eps);

    % Set up pars structure for one prolate of aspect ratio related to
    % 'eps'
    pars=SpheroidalParameters;
    this_u0=1/sqrt(1-1/R^2);
    pars.u0=this_u0;
    pars.a=1/this_u0;
    pars.matvec_eta = 5;
    pars.isReal=0;
    pars.oblate=0;
    pars.centers=[0 0 0];
    thetas=0;
    phis=0;
    Ri=zeros(3,3,1);
    for i=1:1
        thetai=thetas(i);
        phii=phis(i);
        Riy=[cos(thetai) 0 sin(thetai); 0 1 0; -sin(thetai) 0 cos(thetai)];
        Riz=[cos(phii) -sin(phii) 0; sin(phii) cos(phii) 0; 0 0 1];
        Ri(:,:,i)=Riz*Riy;
    end
    pars.Rmat=Ri;
    pars.thetas=thetas;
    pars.phis=phis;
    pars.sigma=zeros(np,1,1); pars.get_shc();
    
    SM = spheroidalMatVecKernel(pars,'SL',p);
    DM = spheroidalMatVecKernel(pars,'DL',p);
    % Completion matrix operator: sum of 1/||x-c_i|| * integral of each sigma over
    %       respective spheroid surface
    CM=[];
    Yi = prolate_spheroid_shape(p,pars.u0,pars.a,'cart');
    Sns = SurfaceSph(Yi(:));
    W = Sns.geoProp.W;
    [~, gwt]=g_grid(p+1);
    wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
    wt = wt(:)';
    RY = vecnorm(Yi,2,2);
    CM = [CM (1./RY)*(W'.*wt)];

end

function [CM,SM,DM,this_u0] = get_mats_two(eps,p)
    np=2*p*(p+1);
    ns = 2;

    R = 1/(eps);
    this_u0=1/sqrt(1-1/R^2);

    pars = set_up(3,p,1,R);
    Y = pars.get_X;
    [~, gwt]=g_grid(p+1);
    wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
    wt = wt(:)';

    SM = spheroidalMatVecKernel(pars,'SL',p);
    DM = spheroidalMatVecKernel(pars,'DL',p);

    % Completion matrix operator: sum of 1/||x-c_i|| * integral of each sigma over
    %       respective spheroid surface
    CM=[];
    for i=1:ns
        if ~pars.oblate(i)
            Yi = prolate_spheroid_shape(p,pars.u0(i),pars.a(i),'cart');
        else
            Yi = oblate_spheroid_shape(p,pars.u0(i),pars.a(i),'cart');
        end
        ci=pars.centers(i,:);
        
        Sns = SurfaceSph(Yi(:));
        W = Sns.geoProp.W;
    
        RY = vecnorm(Y-repmat(ci,np*ns,1),2,2);
        CM = [CM (1./RY)*(W'.*wt)];
    end
end

% Set up pars structure for two-prolate system with 'mode' configuration
%   mode: 1 (parallel major axis), 2 (T shape), 3 (parallel minor axis).
%   p: order of discretization on each.
%   betw_dist: separation from tip to tip
%   aR: aspect ratio of each prolate
function pars = set_up(mode,p,betw_dist,aR)
    np=2*p*(p+1);
    ns=2;
    
    pars=SpheroidalParameters;
    pars.matvec_eta = 5;
    pars.isReal=0;
    % pars.u0=[1.1,1.1];

    this_u0=1/sqrt(1-1/aR^2);
    pars.u0=[this_u0,this_u0];

    pars.a=1./pars.u0;
    pars.oblate=[0,0];

    % axis lengths for prolates
    Aaxis=sqrt(pars.u0.^2-1)./pars.u0;
    Caxis=[1,1];

    if mode==1 % ||
        pars.centers=[0 0 0;Aaxis(1)+Aaxis(2)+betw_dist 0 0];
        thetas = [0 0];
        phis = [0 0];
    elseif mode==2 % |-
        pars.centers=[0 0 0;Aaxis(1)+Caxis(2)+betw_dist 0 0];
        thetas = [0 pi/2];
        phis = [0 0];
    elseif mode==3 % - -
        pars.centers=[0 0 0;Caxis(1)+Caxis(2)+betw_dist 0 0];
        thetas = [pi/2 pi/2];
        phis = [0 0];
    else
        Error("mode not implemented. Choose 1 for ||, 2 for -|, 3 for - -")
    end

    Ri=zeros(3,3,ns);
    for i=1:ns
        thetai=thetas(i);
        phii=phis(i);
        Riy=[cos(thetai) 0 sin(thetai); 0 1 0; -sin(thetai) 0 cos(thetai)];
        Riz=[cos(phii) -sin(phii) 0; sin(phii) cos(phii) 0; 0 0 1];
        Ri(:,:,i)=Riz*Riy;
    end
    pars.Rmat=Ri;
    pars.thetas=thetas;
    pars.phis=phis;
    
    pars.sigma=zeros(np,1,ns); pars.get_shc();

end

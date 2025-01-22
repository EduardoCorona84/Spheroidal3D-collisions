% Exterior Dirichlet problem on two spheroids of u0=1.1, in three
% different configurations: 
%   long axis parallel (mode 1, ||); 
%   long axes perpendicular (mode 2, |-); 
%   short axes parallel (mode 3, --); 
% Dirichlet boundary conditions set as no-slip velocity on surface of 
% spheroids (u=-u_inf on surface), with u_inf=1. 

%% -----------------------------------------------------
% S and C solution comparison (when using higher p, solutions tend closer)

p=8; np=2*p*(p+1);
test_mode=1;
aRlist = [1.1,2,4];
naR = numel(aRlist);

sep_dist=0.05;

params=set_up(test_mode,p,sep_dist);

parsC=copy(params);
parsS=copy(params);

ns = size(params.centers,1);

[sigmaC,condC] = solve_noslip(p,parsC,0);
[sigmaS,condS] = solve_noslip(p,parsS,1);
% new sigma set in parsC/S

% Grid points to evaluate 
[XX,YY]=meshgrid(-0.6:0.02:1.6,-0.6:0.02:0.6);
Xeval=[XX(:),YY(:),zeros(size(XX(:),1),1)];

is_clear=ones(size(Xeval,1),1);
for i=1:2
    cntr=params.centers(i,:);
    u0=params.u0(i);
    a=params.a(i);
    ob=params.oblate(i);
    rot_mat=params.Rmat(:,:,i);
    is_clear = is_clear.*check_clear(cntr,u0,a,ob,rot_mat,Xeval);
end
Xeval=Xeval(is_clear==1,:);

solnC = eval_noslip(p,parsC,0,Xeval);
solnS = eval_noslip(p,parsS,1,Xeval);
solnC=real(solnC);
solnS=real(solnS);

max(abs(solnC-solnS)) % CHECKED: as p increases, C and S result in the same solution.
% ---------------------------------------------------------------

%% ---------------------------------------------------------------
% Self-convergence check for useS=true
plist=[4,8,12,16,24];
err_list=zeros(2,length(plist)-1);
test_mode=1; sep_dist=0.5;

for ii=1:length(plist)
    p=plist(length(plist)-ii+1);

    params=set_up(test_mode,p,sep_dist);

    parsC=copy(params);
    parsS=copy(params);

    [sigmaC,condC] = solve_noslip(p,parsC,0);
    [sigmaS,condS] = solve_noslip(p,parsS,1);

    solnC = eval_noslip(p,parsC,0,Xeval);
    solnS = eval_noslip(p,parsS,1,Xeval);
    solnC=real(solnC);
    solnS=real(solnS);

    if p==24
        truesolnC = solnC;
        truesolnS = solnS;
    else
        err_list(1,length(plist)-ii+1)=max(log10(abs((solnC-truesolnC)./max(truesolnC))));
        err_list(2,length(plist)-ii+1)=max(log10(abs((solnS-truesolnS)./max(truesolnS))));
    end

end
plot(plist(1:end-1),err_list(1,:),'-*'); hold on;
plot(plist(1:end-1),err_list(2,:),'-*'); hold off;
legend("D+C","D+S");
ylabel("log10(rel error)")
xlabel('p')
% ---------------------------------------------------------------

%% ----------------------------------------------------------------
% density change as sep_dist decreases
% clear;
test_mode=1;
ns = 2;

plist=[4,8,16];
aRlist=[1.1,2,4];
% sep_dist=10.^(0:-1:-2); % for condition number and density peak
% cond_list=zeros(2,length(sep_dist));
% den_list=zeros(2,length(sep_dist)); % only collect density at the closest point, and only for D+S

f=figure(1);
f.Position=[100 300 1400 400];
for jj=1:length(aRlist)

    aR = aRlist(jj);
    dist=1./aR;

    for ii=1:length(plist)
        p=plist(ii);

        np=2*p*(p+1);

        obl=1;

        params=set_up(aR,test_mode,p,dist,obl);

        if ~obl
            Aaxis=sqrt(params.u0.^2-1)./params.u0;
            Caxis=[1,1];
        else
            Caxis=params.u0./sqrt(params.u0.^2+1);
            Aaxis=[1,1];
        end
    
        % parsC=copy(params);
        parsS=copy(params);

        C = params.centers; ns = size(params.centers,1);
        m = 21; %number of charges per spheroid 

        Cch = repmat(params.centers',1,m)' + 0.8/aR.*(rand(m*ns,3)*2-1);
        
        Sch = 5*rand(ns*m,1);
        
        % [sigmaC,condC] = solve_noslip(p,parsC,0);
        [sigmaS,condS] = solve_noslip(p,parsS,1,Cch,Sch);

        % % plot density on surface
        % figure(1)
        % subplot(1,length(sep_dist),ii)
        % addpath('./spheroidal/plot/utils/');
        % cmp = getPyPlot_cMap('rainbow', [], [], '"/opt/homebrew/bin/python3"');
        % Y=params.get_X;
        % for sphind=1:ns
        %     Yi=Y((sphind-1)*np+1:sphind*np,:);
        %     shades=parsS.sigma(:,:,sphind);
        %     if norm(real(shades)-shades)>1e-8
        %         fprintf("\n ERROR: density imaginary.");
        %     end
        % 
        %     % TODO: adjust camera angle for far config.
        %     plotb([Yi(:,1);Yi(:,2);Yi(:,3)],real(shades)); hold on;
        % 
        % end
        % colormap(cmp);
        % cb=colorbar; hold off;
        % if ii==1
        %     set(cb,'Position',[0.365 0.25 0.01 0.3])
        % elseif ii==2
        %     set(cb,'Position',[0.65 0.25 0.01 0.3])
        % else
        %   set(cb,'Position',[0.93 0.25 0.01 0.3])
        % end
        % title("charge density")
        % 
        % % cond_list(1,ii)=condC;
        % cond_list(2,ii)=condS;
    % 
        % % density at the point closest to each other.
        % if test_mode==1
        %     den_list(1,ii)=real(sigmaS(ceil((p+1)/2)));
        %     den_list(2,ii)=real(sigmaS(ceil((p+1)/2)+(p+1)*p+np));
        % elseif test_mode==2
        %     den_list(1,ii)=real(sigmaS(ceil((p+1)/2)));
        %     den_list(2,ii)=real(sigmaS(1+np));
        % elseif test_mode==3
        %     den_list(1,ii)=real(sigmaS(1+p));
        %     den_list(2,ii)=real(sigmaS(1+np));
        % end
    
        % plot u(x) to check
        % figure(2)
        % subplot(1,length(sep_dist),ii)
        [XX,YY] = meshgrid(-1:0.01:4.5,-1:0.01:1);
        Xeval = [XX(:),zeros(size(XX(:),1),1),YY(:)];
        is_clear1 = check_clear(parsS.centers(1,:),parsS.u0(1),parsS.a(1),parsS.oblate(1),parsS.Rmat(:,:,1),Xeval);
        is_clear2 = check_clear(parsS.centers(2,:),parsS.u0(2),parsS.a(2),parsS.oblate(1),parsS.Rmat(:,:,2),Xeval);
        is_clear = is_clear1 .* is_clear2;
        Xeval = Xeval(is_clear==1,:);
    
        usol = eval_noslip(p,parsS,1,Xeval);
        UU2 = zeros(size(XX(:),1),1);
        UU2(is_clear==1) = usol(:);
        UU2 = reshape(UU2,size(XX,1),size(XX,2));
        UU2 = real(UU2);
    
        truesolnS = -PtChargePotential(Sch,Cch,Xeval);
    
        ZZ = zeros(size(XX));
        z = ZZ(:); z(is_clear==1) = log10(abs(usol-truesolnS)./max(abs(truesolnS)));
        z(z==0)=NaN;
        ZZ2 = reshape(z,size(XX,1),size(XX,2));
        % errz{ii}=z; 
        % 
        figure(1);
        subplot(3,3,(jj-1)*3+ii);
        contourf(XX,YY,ZZ2,-10:1:0); axis equal; clim([-10,0]);  colorbar;
        hold on; 
        plot(Cch(:,1),Cch(:,3),'xk');
        hold off; 
        % keyboard;

    end

end

% figure;
% % semilogx(sep_dist,cond_list(1,:),'-*'); hold on;
% semilogx(sep_dist,cond_list(2,:),'-o'); hold off;
% ylabel("condition number");
% xlabel("log10(separation of particles)");
% % legend("D+C","D+S");
% % 
% figure;
% semilogx(sep_dist,den_list(1,:),'-*'); hold on;
% semilogx(sep_dist,den_list(2,:),'-o'); hold off;
% ylabel("sigma");
% xlabel("log10(separation of particles)");
% legend("sph1","sph2");
% title("density at point closest to each other")

% ---------------------------------------------------------------


function pars = set_up(aR,mode,p,betw_dist,obl)

    if nargin==4
        obl=0;
    end

    np=2*p*(p+1);
    ns=2;
    
    pars=SpheroidalParameters;
    pars.matvec_eta = 5;
    pars.isReal=0;

    if obl
        this_u0=1./sqrt(aR.^2-1);
    else
        this_u0=1/sqrt(1-1/aR^2);
    end
   
    pars.u0=[this_u0,this_u0];

    if ~obl
        pars.a=1./pars.u0;
        pars.oblate=[0,0];
        % axis lengths for prolates
        Aaxis=sqrt(pars.u0.^2-1)./pars.u0;
        Caxis=[1,1];
    else
        pars.a = 1./sqrt(pars.u0.^2+1);
        pars.oblate=[1,1];
        Caxis=pars.u0./sqrt(pars.u0.^2+1);
        Aaxis=[1,1];
    end

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
    
    pars.sigma=zeros(np,1,ns);

end


function [sigma_vec,condK]=solve_noslip(p,pars,useS,Xptch,ptch)
    % Solve for density on surface of particles in pars, using completion
    % term C or S; report sigma as a vector and condition number of K
    % solved.

    c = pars.centers;
    ns = size(c,1);
    np = 2*p*(p+1);
    Y = pars.get_X;
    [~, gwt]=g_grid(p+1);
    wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
    wt = wt(:)';
    
    % axis lengths for prolates
    if pars.oblate(1)==0
        Aaxis=sqrt(pars.u0.^2-1)./pars.u0;
        Caxis=[1,1];
        aR = Caxis(1)/Aaxis(1); 
    else
        Caxis=pars.u0./sqrt(pars.u0.^2+1);
        Aaxis=[1,1];
        aR = Aaxis(1)/Caxis(1); 
    end

    % E=sqrt(1-ep^2); 
    % SA = 2*pi*(ep)*(1+(aR/E)*asin(E));
    % etaC = 1/SA; 

    etaS = 1;
    etaC = 1;

    if useS
        CM = spheroidalMatVecKernel(pars,'SL',p);
        CM = etaS * CM;
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
        CM = etaC * CM;
    end

    % Impose uniform flow at far field to the right and no-slip condition on surface: 
    % (1/2I+K)*sigma=-u_inf
    % truesolnSurf=1*ones(np*ns,1);
    % %%%%%%%%% If one ptcl is moving toward the other, u2 = -u_inf + u_rel
    % %%%%%%%%% Note: hardcoded ns=2.
    % truesolnSurf(np+1:end) = 0.7*ones(np,1);
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % % Point charge in between spheroids
    % Xptch = sum(pars.centers,1)/2;
    % ptch = 5;
    truesolnSurf = -PtChargePotential(ptch,Xptch,Y);

    % Set of point charges randomly distributed around spheroids, creates
    % background potential. Impose Dirichlet BC s.t. total charge on
    % surface is 0: 
    % 0 = u_from_surface + u_from_ptch
    % rng(1);
    % nch = 20;
    % nelig = 0;
    % Xch = []; pch=[];
    % while nelig<nch
    %     Xptch = rand(nch,3)*10-3;
    %     ptch = rand(nch,1)*5;
    %     is_clear1 = check_clear(pars.centers(1,:),pars.u0(1),pars.a(1),0,pars.Rmat(:,:,1),Xptch);
    %     is_clear2 = check_clear(pars.centers(2,:),pars.u0(2),pars.a(2),0,pars.Rmat(:,:,2),Xptch);
    %     is_clear = is_clear1 .* is_clear2;
    %     nelig = nelig + sum(is_clear);
    %     Xch = [Xch;Xptch(is_clear==1,:)];
    %     pch = [pch;ptch(is_clear==1)];
    % end
    % truesolnSurf = -PtChargePotential(pch,Xch,Y);

    % Construct DL on-surface matrices
    DM = spheroidalMatVecKernel(pars,'DL',p);
    
    % Final operator: 1/2 I + DM + CM
    K = .5*eye(ns*np) + DM + CM;
    
    condK=cond(K);
    
    sigma_vec=K\truesolnSurf;
    sigma = reshape(sigma_vec,np,[],ns);
    pars.sigma = sigma;
end

function soln = eval_noslip(p,pars,useS,Xeval)
    
    c = pars.centers;
    ns = size(c,1);
    np = 2*p*(p+1);
    sigma=pars.sigma;

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

    % soln=spheroidalMatVec(pars,'DL',Xeval) + pars.etaC.*C;
    soln=spheroidalMatVec(pars,'DL',Xeval) + C;

end

function is_clear = check_clear(cntr,u0,a,ob,rot_mat,xtrg)
    % check if points xtrg are inside a given spheroid

    % input:
    %       cntr: 1 x 3 position of spheroid.
    %       u0, a, ob: parameters of the spheroid i.
    %       xtrg: n x 3 points to determine.

    % output: 
    %       is_clear: n x 1 boolean of whether each point is "clear", i.e.
    %       outside the given spheroid.

    xtrg_rot=(rot_mat'*(xtrg-cntr)')'+cntr;

    xdis=xtrg_rot(:,1)-cntr(1);
    ydis=xtrg_rot(:,2)-cntr(2);
    zdis=xtrg_rot(:,3)-cntr(3);
    if ob
        formula=(xdis.^2+ydis.^2)./a^2./(u0^2+1)+zdis.^2./a^2./u0^2;
    else
        formula=(xdis.^2+ydis.^2)./a^2./(u0^2-1)+zdis.^2./a^2./u0^2;
    end
    is_clear=formula>1;
end


function pcp = PtChargePotential(ptch,Xptch,Y)
M=length(ptch);
% np=length(Y);
npch=size(Y,1);
Rptch=zeros(npch,M);
for i=1:M
    Rptch(:,i)=sqrt(sum((Xptch(i,:)-Y).^2,2));
end
pcp=1./(4*pi*Rptch)*ptch; 
end



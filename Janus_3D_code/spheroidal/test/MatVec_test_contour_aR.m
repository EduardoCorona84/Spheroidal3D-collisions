% Exterior Dirichlet problem on two prolate spheroids of u0=1.1, in three
% different configurations: 
%   long axis parallel (mode 1, ||); 
%   long axes perpendicular (mode 2, |-); 
%   short axes parallel (mode 3, --); 
%   four spheroids (mode 4)
% Dirichlet boundary conditions set as no-slip velocity on surface of 
% spheroids (u=-u_inf on surface), with u_inf=1. 

% Report error contour plot as distance decreases and aR increases. Report
% also GMRES iterations for matrix solve.

rng(13.5);
% plist=[4,8,16,32];
plist=[4,8,16];
aRlist=[1.1,2,4,8,16,32,64];
naR = length(aRlist); 
test_mode=4; ns=4;

err_list=zeros(2,length(plist)-1);
errz = cell(length(plist),1); 
GMiter_table = zeros(naR*3,4); % columns: operators; rows: different R.

f=figure(1);
f.Position=[100 300 1200,700];
for jj=1:naR
    aR = aRlist(jj);
    sep_dist=2/aR;

    for ii=1:length(plist)
        p=plist(ii);
        params=set_up(test_mode,p,sep_dist,aR,ns);

        if ii==1

            Cpot = [0.1 1 4]; Spot=0; 
            phi_inf = @(X) evalfpot(X,Cpot,Spot); 
            np=2*p*(p+1); 
            
            Aaxis=sqrt(params.u0.^2-1)./params.u0; % prolate only
            Caxis=ones(1,ns);
            
            aspectR = Caxis./Aaxis; 
            ep = 1./aspectR;
            
            if aspectR(1) > 3
                etaS=1./(2*ep.*log(1./ep)); 
            else
                etaS=0.5; 
            end
            
            A = ep; E=sqrt(1-ep.^2); 
            SA = 2*pi*(ep).*(1+(aspectR./E).*asin(E));
            etaC = 1./SA; 
            
            C = params.centers; ns = size(params.centers,1);
            m = 21; %number of charges per spheroid 
            P = 2*rand(ns*m,3)-1; 
            Pr = min(Aaxis)*repmat(2*rand(ns*m,1)-1,1,3).*(P./repmat(sqrt(sum(P.^2,2)),1,3)); 
            
            if test_mode==1 || test_mode==4
                Pr(1,:)=0.8*min(Aaxis)*[1,0,0];
                Pr(2,:)=-0.7*min(Aaxis)*[1,0,0];
                Pr(3,:)=0.75*min(Caxis)*[0,0.01,1];
                Pr(4,:)=-0.72*min(Caxis)*[0.01,0,1];
                Pr(5,:)=-0.70*min(Caxis)*[0,0.01,1];
                Pr(6,:)=0.75*min(Caxis)*[0.02,0,1];
            elseif test_mode==3
                Pr(1,:)=0.8*min(Caxis)*[1,0,0];
                Pr(2,:)=-0.7*min(Caxis)*[1,0,0];
                Pr(3,:)=0.75*min(Aaxis)*[0,0.01,1];
                Pr(4,:)=-0.72*min(Aaxis)*[0,0.01,1];
                Pr(5,:)=-0.70*min(Aaxis)*[0,0.01,1];
                Pr(6,:)=0.75*min(Aaxis)*[0,0.02,1];
            end
            
            Cch = repmat(params.centers',1,m)'+Pr; 
            Sch = 5*rand(ns*m,1);

        end

        if jj==1
            % Grid points to evaluate
            N = 100000; n=round(sqrt(N)); 

            switch test_mode
                case 3
                    gx = linspace(params.centers(1,1)-1.5,params.centers(ns,1)+1.5,n) ; 
                    my = min(sqrt(8/aR),sqrt(2));
                    gy = linspace(-my,my,n); gy = sign(gy).*(gy.^2); 
                    [XX,YY]=meshgrid(gx,gy);
                    Xeval=[XX(:),YY(:),zeros(size(XX(:),1),1)];
                case 1
                    %mx = sqrt(3/aR + sep_dist/2); %min(sqrt(8/aR),sqrt(2)); 
                    %cx=1/aR + sep_dist/2;
                    %gx = linspace(-mx+cx,cx+mx,np); gx = sign(gx-cx).*((gx-cx).^2) + cx;  
                    gx = linspace(params.centers(1,1)-2/aR,params.centers(ns,1)+2/aR,n);
                    gy = linspace(-1.5,1.5,n); 
                    [XX,YY]=meshgrid(gx,gy);
                    Xeval=[XX(:),zeros(size(XX(:),1),1),YY(:)];
                case 2
                    gx = linspace(-1.6,3.6+sep_dist,n) ; 
                    my = min(sqrt(8/aR),1);
                    gy = linspace(-my,my,n); gy = sign(gy).*(gy.^2);
                    [XX,YY]=meshgrid(gx,gy);
                    Xeval=[XX(:),YY(:),zeros(size(XX(:),1),1)];
                case 4
                    gx = linspace(min(params.centers(:,1))-2/aR,max(params.centers(:,1))+2/aR,n);
                    gy = linspace(min(params.centers(:,3))-1.5,max(params.centers(:,3))+1.5,n); 
                    [XX,YY]=meshgrid(gx,gy);
                    Xeval=[XX(:),zeros(size(XX(:),1),1),YY(:)];
            end

        end
    
        if ii==3 % only store GMRES values for p=16
            % without scaling
            parsC1= copy(params);
            parsS1= copy(params);
            [sigmaS1,condS1,GMS1] = solve_noslip(p,parsS1,1,Cch,Sch,1,phi_inf);
            [sigmaC1,condC1,GMC1] = solve_noslip(p,parsC1,0,Cch,Sch,1,phi_inf);
            
            if GMC1(1)~=1 || GMS1(1)~=1
                error("\n number of inner iter not 1.")
            end
        end
        
        % with scaling
        parsC = copy(params);
        parsS = copy(params);
        [sigmaC,condC,GMC] = solve_noslip(p,parsC,0,Cch,Sch,etaC,phi_inf);
        [sigmaS,condS,GMS] = solve_noslip(p,parsS,1,Cch,Sch,etaS,phi_inf);

        if GMC(1)~=1 || GMS(1)~=1
            error("\n number of inner iter not 1.")
        end

        if ii==3
            fprintf("\n cond(C)=%d, cond(etaC)=%d, cond(S)=%d, cond(etaS)=%d \n",condC1,condC,condS1,condS)

            GMiter_table(jj,1) = GMC1(2);
            GMiter_table(jj,2) = GMC(2);
            GMiter_table(jj,3) = GMS1(2);
            GMiter_table(jj,4) = GMS(2);
        end

        if jj<4 % only plot contour if R=1.1,2,4.
            is_clear=ones(size(Xeval,1),1);
            for i=1:ns
                cntr=params.centers(i,:);
                u0=params.u0(i);
                a=params.a(i);
                ob=params.oblate(i);
                rot_mat=params.Rmat(:,:,i);
                is_clear = is_clear.*check_clear(cntr,u0,a,ob,rot_mat,Xeval);
            end
            Xeval_loc=Xeval(is_clear==1,:);

            solnC = eval_noslip(p,parsC,0,Xeval_loc,etaC,phi_inf);
            solnS = eval_noslip(p,parsS,1,Xeval_loc,etaS,phi_inf);
            solnC=real(solnC);
            solnS=real(solnS);

            truesolnC = evalfpot(Xeval_loc,Cch,Sch) + phi_inf(Xeval_loc);
            truesolnS = truesolnC; 

            err_list(1,ii)=mean(log10(abs((solnC-truesolnC)./max(truesolnC))));
            err_list(2,ii)=mean(log10(abs((solnS-truesolnS)./max(truesolnS))));

            ZZ = zeros(size(XX));
            z = ZZ(:); z(is_clear==1) = log10(abs(solnS-truesolnS)./max(abs(truesolnS)));
            z(z==0)=NaN;
            ZZ2 = reshape(z,size(XX,1),size(XX,2));
            errz{ii}=z; 

            figure(1);
            % subplot(3,4,(jj-1)*4+ii);
            subplot(3,3,(jj-1)*3+ii);
            contourf(XX,YY,ZZ2,-10:1:0); axis equal; clim([-10,0]); 
            if jj==1 && ii==numel(plist)
                cb=colorbar;
                set(cb,'Position',[0.92 0.1 0.02 0.8])
            end
            hold on; 
            plot(Cch(:,1),Cch(:,3),'xk');
            hold off; 
            %surf(XX,YY,ZZ2); 
            %clim([-10 0]); shading interp;
            pause(0.00001);
        end
    end

    % figure; 
    % plot(plist,err_list(1,:),'-*'); hold on;
    % plot(plist,err_list(2,:),'-*'); hold off;
    % legend("D+C","D+S");
    % ylabel("log10(rel error)")
    % xlabel('p')
    % title("aR="+string(aR))

end
% keyboard;

% ---------------------------------------------------------------



function Kx = evalfpot(Xtrg,Cch,Sch)

    %p = pars.p; C = pars.centers; u0=pars.u0; a=pars.a; 
    %ns=size(C,1); np=2*p*(p+1); 
    %X = zeros(np*ns,3); Sns=cell(ns,1); W2 = zeros(np*ns,1); 

    pot='SL_L_3D'; 
    KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-14,2,400,1);
    KEparams.dim = 3;
    KEparams.W2=ones(size(Cch,1),1).'; 
    %{
    [~, gwt]=g_grid(p+1);
    wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
    wt = wt(:);

    for i=1:ns
        Xsph = prolate_spheroid_shape(p,u0(i),a(i));
        X((1:np)+ns*(i-1),:) = Xsph + repmat(C(i,:),np,1);
        Sns{i}=SurfaceSph(X((1:np)+ns*(i-1),:));
        Wns = Sns{i}.geoProp.W; Wns= Wns.*wt;
        W2((1:np)+ns*(i-1)) = Wns; 
    end

    KEparams.X = X;
    KEparams.W2 = Wns.';
    %}

    Kx = Kernel_Eval(Xtrg,Cch,KEparams)*Sch;

end

function pars = onesph_set_up(p,aR,oblate)
    np=2*p*(p+1);
    ns=1;
    
    pars=SpheroidalParameters;
    pars.matvec_eta = 5;
    pars.isReal=0;
    this_u0=1/sqrt(1-1/aR^2);
    pars.u0=this_u0;
    pars.a=1./pars.u0;
    pars.oblate=oblate;

    % axis lengths for prolates
    Aaxis=sqrt(pars.u0.^2-1)./pars.u0;
    Caxis=[1,1];

    aR = Caxis(1)/Aaxis(1); ep=1/aR; 
    pars.etaS = 1/(2*ep*log(1/ep));

    E=sqrt(1-ep^2); 
    SA = 2*pi*(ep)*(1+(aR/E)*asin(E));
    pars.etaC = 1/SA; 

    pars.centers=[0 0 0];
    thetas=0; phis=0; 

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

function pars = set_up(mode,p,betw_dist,aR,ns)
    np=2*p*(p+1);
    if nargin<5
        ns=2;
    end
    
    pars=SpheroidalParameters;
    pars.matvec_eta = 5;
    pars.isReal=0;
    this_u0=1/sqrt(1-1/aR^2);
    pars.u0=this_u0*ones(1,ns);
    pars.a=1./pars.u0;
    pars.oblate=zeros(1,ns);

    % axis lengths for prolates
    Aaxis=sqrt(pars.u0.^2-1)./pars.u0;
    Caxis=ones(1,ns);

    aR = Caxis(1)/Aaxis(1); ep=1/aR; 
    pars.etaS = 1/(2*ep*log(1/ep));

    E=sqrt(1-ep^2); 
    SA = 2*pi*(ep)*(1+(aR/E)*asin(E));
    pars.etaC = 1/SA; 

    if mode==1 % ||
        pars.centers=[0 0 0;Aaxis(1)+Aaxis(2)+betw_dist 0 0];
        thetas = [0 0];
        phis = [0 0];
        if ns>2
            for k=3:ns
                pars.centers=[pars.centers; ...
                    (pars.centers(k-1,1)+Aaxis(k-1)+betw_dist+Aaxis(k)) 0 0];
                thetas(k)=0; phis(k)=0;
            end
        end

    elseif mode==2 % |-
        pars.centers=[0 0 0;Aaxis(1)+Caxis(2)+betw_dist 0 0];
        thetas = [0 pi/2];
        phis = [0 0];
        if ns>2
            for k=3:ns
                pars.centers=[pars.centers; ...
                    (pars.centers(k-1,1)+Axis(k-1)+betw_dist+Aaxis(k)) 0 0];
                thetas(k)=0; phis(k)=0; 
            end
        end
    elseif mode==3 % - -
        pars.centers=[0 0 0;Caxis(1)+Caxis(2)+betw_dist 0 0];
        thetas = [pi/2 pi/2];
        phis = [0 0];
        if ns>2
            for k=3:ns
                pars.centers=[pars.centers; ...
                    (pars.centers(k-1,1)+Caxis(k-1)+betw_dist+Caxis(k)) 0 0];
                thetas(k)=pi/2; phis(k)=0; 
            end            
        end
    elseif mode==4
        % lattice 
        ns2 = ceil(sqrt(ns))^2; %nearest square 
        ms = sqrt(ns); 
        latc = zeros(ns2,3); 
        Cz = 0; 

        for i=1:ms
            Cx = 0;   
            for j=1:ms
                k = (i-1)*ms+j; 
                latc(k,:) = [Cx 0 Cz];
                if j<ms
                    Cx = Cx + Aaxis(k)+betw_dist+Aaxis(k+1);
                end
            end
            if i<ms
                Cz = Cz + Caxis(k)+betw_dist+Caxis(k+1);
            end
        end

        tt = linspace(0,pi/2,ns2)';
        latc(:,1) = latc(:,1) + 0.5*(betw_dist/2)*cos(tt);
        latc(:,3) = latc(:,3) + 0.5*(betw_dist/2)*sin(tt);

        pars.centers = latc(1:ns,:); 
        thetas = zeros(1,ns);
        phis = zeros(1,ns); 
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


function [sigma_vec,condK,iter]=solve_noslip(p,pars,useS,Cch,Sch,eta,phi_inf)
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
    
    if useS
        pars.etaS = mean(eta); 
        CM = pars.etaS*spheroidalMatVecKernel(pars,'SL',p);
    else
        pars.etaC = mean(eta); 
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
        CM = pars.etaC*CM; 
    end

    % Impose uniform flow at far field to the right and no-slip condition on surface: 
    % (K)*sigma=-u_inf
    truesolnSurf = evalfpot(Y,Cch,Sch) + phi_inf(Y);
    % truesolnSurf=4*(2*rand(np*ns,1)-1);

    % Construct DL on-surface matrices
    DM = spheroidalMatVecKernel(pars,'DL',p);
    
    % Final operator: 1/2 I + DM + CM
    K = .5*eye(ns*np) + DM +CM;
    
    condK=cond(K);
    
    % keyboard;

    % sigma_vec_prev=K\truesolnSurf;

    warning('off')

    [sigma_vec,flag,relres,iter] = gmres(K,truesolnSurf,[],1e-10,1e6);
    % fprintf("\n iters = %d", iter)
    if flag
        error("\n GMRES did not converge.")
    % else
    %     fprintf("\n GMRES residue: %.3g", relres);
    end

    % fprintf("\n norm between backslash and gmres: %.3g\n", norm(sigma_vec_prev-sigma_vec));
    % keyboard;

    sigma = reshape(sigma_vec,np,[],ns);
    pars.sigma = sigma;
end

function soln = eval_noslip(p,pars,useS,Xeval,eta,phi_inf)

    c = pars.centers;
    ns = size(c,1);
    np = 2*p*(p+1);
    sigma=pars.sigma;

    if useS
        pars.etaS = mean(eta); 
        C = pars.etaS*spheroidalMatVec(pars,'SL',Xeval);
    else
        pars.etaC=mean(eta); 
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
        C = pars.etaC*C; 
    end

    soln=spheroidalMatVec(pars,'DL',Xeval) + C;
    %soln = soln + phi_inf(Xeval); 

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
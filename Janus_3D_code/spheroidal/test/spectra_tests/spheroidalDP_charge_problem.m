%{
    Helper code to test spheroidalDP.
    This code tests the code on a simple charge problem.
%}
function [soln, fluxsoln, truesoln, trueflux, sigma_vec,condK] = spheroidalDP_charge_problem(p, eta, ns, u0, target_distances, plt, neumann, interior)
    %{
        Inputs:
            - p -> order
            - eta -> 
            - ns -> number of spheroids
            - u0 -> eccentricities of spheroids
            - target_distances ->
            - plt (true/false) -> whether to plot the points; a holdover 
            from old code.
            - neumann (true/false) -> whether to do a neumann problem or 
            not; note that we get an integral equation of the first kind, 
            which is known to have bad conditioning
            - interior (true/false) -> whether to do an exterior problem or 
            an interior problem

        Outputs:
            soln ->
            truesoln -> analytical solution
            truefluxSurf -> 
            sigma_vec -> 
            condK -> condition number of BIE matrix
        
        The idea is:
        (1) set up system of 3 spheroids, 2 close and 1 far
        (2) put a few point charges around the center of each spheroid.
              - Calculate the potential explicitly.
              - Determine potential on the boundary of each spheroid for BCs, f
        (3) Use MatVec self-evaluation (no X) to construct matrix operator D.
        (4) Solve (1/2*I + D + Completion)*sigma = f
              * Completion term will need some care, since it's multiple particles.
              For each particle it will be something like sum_i[ 1/|| x- c_i||
              integral(sigma_i) ]
        (5) Compare sigma with true point charge potential at targets that
            are 'target_distances' away from surface of each spheroid.
    %}

    %%% Backwards compatibility
    S_scale = 1;
    Xeval_p = p;
    mix_obl = false;
    useS = false;

    %%% Set up spheroid system
    np=2*p*(p+1);
    
    if length(u0)~=ns
        fprintf("\n input u0 length does not match number of spheroids indicated.");
        ns=length(u0);
    end
    
    pars=SpheroidalParameters;
    pars.matvec_eta = eta;
    pars.isReal = true;
    pars.u0=u0;

    if ns==3
        if mix_obl
            obl = [1,0,1];
        else
            obl = [0,0,0];
        end
        pars.oblate = obl;
        a = 1./u0;
        a(obl==1) = 1./sqrt(u0(obl==1).^2+1);
        pars.a = a;
        
        pars.centers = [0 0 0; 5 0 0; 3.2 3.2 3.2];
        % pars.centers = [0 0 0; 2.5 0 0; 1.6 1.6 1.6];
        thetas = [0 pi/10 5*pi/3];
        phis = [0 0 pi/5];
    else
        if ns~=1
            fprintf("\n Number of spheroids given not implemented here. Setting ns=1.")
            ns=1;
        end

        if mix_obl
            pars.a = 1./sqrt(u0^2+1);
            pars.oblate=1;
        else
            pars.a=1./u0;
            pars.oblate=0;
        end
        
        pars.centers=[0 0 0];
        thetas=0;
        phis=0;
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
        
    %%% Point charges
    % 2 point charges per spheroid
    nc = 2;
    c = pars.centers;
    Xptch = reshape(repmat(reshape(c',3,1,[]),1,nc),3,[],1)';

    %%% Placement of point charges
    if ns==3
        % 1/2 of the minor radius of the prolate spheroids
        scale = repmat(.5.*pars.a .*sqrt(pars.u0.^2-1),nc,1);
        scale = scale(:);
        d = repmat(scale,1,3) .* (rand(size(Xptch))-.5);
    else
        d_z=0.5.*(rand(nc,1)-.5);
        d_xy=0.01.*(rand(nc,2)-0.5);
        d = [d_xy d_z];
    end

    Xptch = Xptch + d; % Point charge locations, near centers of spheroids
    ptch = (2.*rand(ns*nc,1)-1); % Charge value, between -1 and 1

    ptch = ptch - sum(ptch)/(ns*nc); % Ensures that the BC is 0 for compatibility condition
    
    if plt
        pars.plot(Xptch);
        title("spheroids and point charges")
    end

    Y = pars.get_X;
    [NrY,~]=pars.get_Norm_rot(p); % Necessary to account of the spheroid's rotation
    
    [~, gwt]=g_grid(p+1);
    wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
    wt = wt(:)';
    
    if ~neumann
        if useS
            CM = S_scale.*spheroidalMatVecKernel(pars,'SL',p);
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
    else % Neumann problem
        error("not implemented.");
        if FIRST_KIND_FLAG
            truesolnSurf=PtChargeFlux(ptch, Xptch, Y, NrY);
            K = spheroidalMatVecKernel(pars,'DP',p);
        else
            truesolnSurf=PtChargeFlux(ptch, Xptch,Y,NrY);

            % Construct SP on-surface matrices
            SPM= spheroidalMatVecKernel(pars,'SP',p);
            
            % Final operator: -1/2 I + SPM
            K = -.5*eye(ns*np) + SPM;
        end
    end

    condK = cond(K);
    
    sigma_vec=K\truesolnSurf;
    sigma = reshape(sigma_vec,np,[],ns);
    pars.sigma = sigma;
    
    %%% Target points to test at
    [theta,phi]=gl_grid(Xeval_p);
    v=cos(theta);
    Xcell=cell(1,ns);
    for i=1:ns
        u0=pars.u0(i);
        if ~pars.oblate(i)
            normal=1./sqrt((u0^2-1).*(u0^2-v.^2)).*[u0.*sqrt(u0.^2-1).*sqrt(1-v.^2).*cos(phi) u0.*sqrt(u0.^2-1).*sqrt(1-v.^2).*sin(phi) (u0.^2-1).*v];
            Yi = prolate_spheroid_shape(Xeval_p,pars.u0(i),pars.a(i),'cart');
        else
            normal=1./sqrt((u0^2+1).*(u0^2+v.^2)).*[u0.*sqrt(u0.^2+1).*sqrt(1-v.^2).*cos(phi) u0.*sqrt(u0.^2+1).*sqrt(1-v.^2).*sin(phi) (u0.^2+1).*v];
            Yi = oblate_spheroid_shape(Xeval_p,pars.u0(i),pars.a(i),'cart');
        end
        Xcell{i} = Yi + target_distances(i).*normal;
    end
    Xeval = pars.set_X_targets(Xcell);

    if plt
        pars.plot(Xeval);
        title('spheroids and target points')
    end
    
    if ~neumann
        if useS
            C = S_scale.*spheroidalMatVec(pars,'SL',Xeval);
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
        if FIRST_KIND_FLAG
            soln = spheroidalMatVec(pars,'DL',Xeval);
        else
            soln = spheroidalMatVec(pars,'SL',Xeval);
        end
    end

    truesoln=PtChargePotential(ptch,Xptch,Xeval);
    fprintf('p=%d: DL comparison = %.6e\n',p, norm(truesoln - soln) ./ norm(soln)); 

    %%% Compare flux solution evaluated away from the surface
    trueflux=PtChargeFlux(ptch, Xptch, Xeval, NrY);
    fluxsoln = spheroidalMatVec(pars, 'DP', Xeval, NrY);

    fprintf('p=%d: matvec comparison -- relative error = %.6e\n',p, norm(trueflux-fluxsoln) / norm(trueflux));
end

function pcp = PtChargePotential(ptch,Xptch,Y)
    %{
        Calculates the electric potential due to a set of point charges.
    %}
    M=length(ptch);
    np=size(Y,1);
    Rptch=zeros(np,M);

    % \|x - y\|
    for i=1:M
        Rptch(:,i)=sqrt(sum((Xptch(i,:)-Y).^2,2));
    end

    % \sum_j q_j/\|x - y\|
    pcp=1./(4*pi*Rptch)*ptch; 
end



function flux=PtChargeFlux(ptch,Xptch,Y,NrY)
    M=length(ptch);
    np=length(Y);
    Rptch=zeros(np,M);
    EdotN=zeros(np,M);
    for i=1:M
        % \|x - y\|
        r=Xptch(i,:)-Y;
        Rptch(:,i)=sqrt(sum(r.^2,2));
        
        EdotN(:,i)=dot(NrY,r./Rptch(:,i),2);
    end
    flux=EdotN./(4*pi*Rptch.^2)*ptch;
end

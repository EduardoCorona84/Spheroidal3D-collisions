function [Stk_x,Stk_y,Stk_z]=L2StkMatVec(pars, pot, sigma_x, sigma_y, sigma_z, X, nu_trg)
    %{
        X -> global frame
        nu_trg -> global frame
        densities are in local frame
    %}
    if isempty(pars.u0)
        error("No surface parameter u_0 given")
    end
    if isempty(pars.a)
        error("No scale (a) given")
    end
    if isempty(pars.matvec_eta)
        error("No matvec_eta provided.");
    end

    DEVELOPMENT_FLAG = true;
    DENSE_EVALUATION_FLAG = true; % Determine if FMM is to be used
    u0=pars.u0;
    a=pars.a;
    ns=size(pars.centers,1);
    if_oblate=pars.oblate;
    p=pars.p;
    np=2*p*(p+1);

    if length(u0)==1
        u0=u0*ones(1,ns); 
    end
    
    if length(a)==1
        a=a*ones(1,ns); 
    end

    assert(size(sigma_x, 3) == ns, "Density should be defined on every spheroid surface.");

    if ~DEVELOPMENT_FLAG
    if nargin==4
        error("Not implemented.");
    else
        if (strcmp(pot, 'TLP') && nargin < 7)
            error("Input for normal vector is necessary for the traction of the Stokes SLP.");
        end

        nt=size(X,1);

        [sep,Xt]=pars.separate_targets(X);
        X_spectral = cell(1,ns);
        if strcmp(pot, 'TLP')
            nu_spectral = cell(1,ns);
            Nu_t = pars.get_nu_targets(nu_trg);
        end

        % Separate target points and normal vectors into nearby and far
        for i=1:ns
            is_near = (sep(:,i) == 0);
            X_spectral{i} = Xt(is_near,:,i);
            if strcmp(pot, 'TLP')
                nu_spectral{i} = Nu_t(is_near,:,i);
            end
        end
        
        %%%
        %%% Do "nearby" evaluation spectrally
        %%%
        if strcmp(pot, 'SLP')
            [L2Stk_spectral_cell_x, L2Stk_spectral_cell_y, L2Stk_spectral_cell_z]=L2Stk(X_spectral, pars, sigma_x, sigma_y, sigma_z, ns);
        elseif strcmp(pot, 'DLP')
            [L2Stk_spectral_cell_x, L2Stk_spectral_cell_y, L2Stk_spectral_cell_z]=L2StkDLP(X_spectral, pars, sigma_x, sigma_y, sigma_z, ns);
        elseif strcmp(pot, 'TLP')
            [L2Stk_spectral_cell_x, L2Stk_spectral_cell_y, L2Stk_spectral_cell_z]=L2StkTLP(X_spectral, nu_spectral, pars, sigma_x, sigma_y, sigma_z, ns, false);
        else
            error("Invalid potential given.");
        end

        %%%
        %%% Do "far" evaluation with smooth quadrature
        %%%
        Stk_x=zeros(nt,1); Stk_y=Stk_x; Stk_z=Stk_x;

        % Calculate effect from each particle on each target point
        for i=1:ns
            L2Stk_local_i = [L2Stk_spectral_cell_x{i}, L2Stk_spectral_cell_y{i}, L2Stk_spectral_cell_z{i}];
            L2Stk_global_i = L2Stk_local_i * pars.Rmat(:,:,i)';

            % Add contributions from near interactions
            Stk_x(sep(:,i)==0,:) = Stk_x(sep(:,i)==0,:) + L2Stk_global_i(:,1); 
            Stk_y(sep(:,i)==0,:) = Stk_y(sep(:,i)==0,:) + L2Stk_global_i(:,2); 
            Stk_z(sep(:,i)==0,:) = Stk_z(sep(:,i)==0,:) + L2Stk_global_i(:,3); 
        
            % Get targets for smooth quadrature
            X_smooth = Xt(sep(:,i)==1,:,i);
            nt_smooth=size(X_smooth,1);

            if nt_smooth>0
                fprintf("Kernel_Eval being used in L2StkMatVec...\n");
                if ~if_oblate(i)
                    Xself=prolate_spheroid_shape(p,u0(i),a(i));
                else
                    Xself=oblate_spheroid_shape(p,u0(i),a(i));
                end
                nu_self = get_norm_vecs(p, u0(i), if_oblate(i));
                Sns=SurfaceSph(Xself);

                if strcmp(pot, 'SLP')
                    KEpot = 'SL_Stk_3D';
                elseif strcmp(pot, 'DLP')
                    KEpot = 'DL_Stk_3D';
                elseif strcmp(pot, 'TLP')
                    KEpot = 'TSL_Stk_3D';
                else
                    error("Invalid potential input.");
                end

                KEparams = Kernel_Eval_parameters(KEpot,0,1,1,1,1e-8,2,400,1);
                KEparams.dim = 3; KEparams.mu=1;
                [~, gwt]=g_grid(p+1);
                wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
                wt = wt(:);
                Wns = Sns.geoProp.W; Wns= Wns.*wt;
                Wv = repmat(Wns,1,3)'; Wv=Wv(:);
                Xv = reshape(repmat(Xself,1,3)',3,[])';
                Xtrg_ii = reshape(repmat(X_smooth,1,3)',3,[])';
                KEparams.X = Xv;
                KEparams.W2 = Wv.';
                KEparams.cj=repmat((1:3)',np,1);
                KEparams.ci=repmat((1:3)',nt_smooth,1);

                %%% Normal vector handling
                if strcmp(pot, 'TLP')
                    Nu_smooth = Nu_t(sep(:,i)==1,:,i);
                    KEparams.nor = reshape(repmat(Nu_smooth,1,3)',3,[])';
                elseif strcmp(pot, 'DLP')
                    KEparams.nor = reshape(repmat(nu_self,1,3)',3,[])';
                end

                LP_Kernel = Kernel_Eval(Xtrg_ii,Xv,KEparams);

                sig = reshape([sigma_x(:,:,i),sigma_y(:,:,i),sigma_z(:,:,i)].',[],1);
                LP = LP_Kernel*sig;
        
                LP=reshape(LP,3,[]).';
                Stk_x(sep(:,i)==1,:) = Stk_x(sep(:,i)==1,:) + LP(:,1);
                Stk_y(sep(:,i)==1,:) = Stk_y(sep(:,i)==1,:) + LP(:,2);
                Stk_z(sep(:,i)==1,:) = Stk_z(sep(:,i)==1,:) + LP(:,3);
            end
        end
    end
    end % END DEVELOPMENT FLAG

    if DEVELOPMENT_FLAG
    nt = size(X, 1);
    Stk_total_global = zeros(nt, 3);

    % For every source spheroid...
    for i = 1:ns
        % Get the rotation matrix for the current source
        Ri = pars.Rmat(:,:,i);

        % Transform target points and normals into the local frame of source i
        X_trg_local_i = (X - pars.centers(i,:)) * Ri;
        if strcmp(pot, 'TLP')
            nu_trg_local_i = nu_trg * Ri;
        end

        params_i = pars.copy();
        params_i.u0 = pars.u0(i);
        params_i.a = pars.a(i);
        params_i.oblate = pars.oblate(i);
        params_i.centers = [0 0 0];
        params_i.Rmat = eye(3);

        % Evaluate the potential from source i onto all targets
        if strcmp(pot, 'SLP')
            [vx, vy, vz] = L2Stk({X_trg_local_i}, params_i, sigma_x(:,:,i), sigma_y(:,:,i), sigma_z(:,:,i), 1);
        elseif strcmp(pot, 'DLP')
            [vx, vy, vz] = L2StkDLP({X_trg_local_i}, params_i, sigma_x(:,:,i), sigma_y(:,:,i), sigma_z(:,:,i), 1);
        elseif strcmp(pot, 'TLP')
            % NOTE: The L2StkTLP function needs to accept normals as an argument
            [vx, vy, vz] = L2StkTLP({X_trg_local_i}, {nu_trg_local_i}, params_i, sigma_x(:,:,i), sigma_y(:,:,i), sigma_z(:,:,i), 1, false);
        else
            error("Invalid potential given.");
        end

        Stk_contribution_local_i = [vx{1}, vy{1}, vz{1}];

        % Rotate from local frame i to global frame
        Stk_contribution_global_i = Stk_contribution_local_i * Ri';
        
        Stk_total_global = Stk_total_global + Stk_contribution_global_i;
    end

    Stk_x = Stk_total_global(:,1);
    Stk_y = Stk_total_global(:,2);
    Stk_z = Stk_total_global(:,3);
    end % END DEVELOPMENT FLAG

    if pars.isReal
        Stk_x = real(Stk_x);
        Stk_y = real(Stk_y);
        Stk_z = real(Stk_z);
    end
end
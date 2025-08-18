function [Stk_x, Stk_y, Stk_z] = L2StkDLP(Xeval, pars, sigma_x, sigma_y, sigma_z, ns)
    %{
        An implementation of the Laplace to Stokes double layer potential.

        Inputs
            Xeval       -   target points
            pars        -   parameters needed for to calculate spheroidal Laplace LPs
            sigma_x     -  
            sigma_y     -
            sigma_z     -
            ns          -   number of spheroidal bodies
    %}
    %%% First, let's generate the necessary vectors for the normal derivatives
    % of the Laplace layer potentials.
    if isempty(Xeval) % Self-evaluation
        % For each body, 
        nu_x_spectral=repmat([1,0,0],size(sigma_x,1),size(sigma_x,2),ns);
        nu_y_spectral=repmat([0,1,0],size(sigma_x,1),size(sigma_x,2),ns);
        nu_z_spectral=repmat([0,0,1],size(sigma_x,1),size(sigma_x,2),ns);
    else % All to targets
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
    end

    %% Calculate associated layer potentials
    % i = 1
    pars.sigma = sigma_x; pars.get_shc();
    [DPsigXdX,DPsigXdY,DPsigXdZ] = spheroidalDP(pars,Xeval,nu_x_spectral,nu_y_spectral,nu_z_spectral);

    % i = 2
    pars.sigma = sigma_y; pars.get_shc();
    [DPsigYdX,DPsigYdY,DPsigYdZ] = spheroidalDP(pars,Xeval,nu_x_spectral,nu_y_spectral,nu_z_spectral);

    % i = 3
    pars.sigma = sigma_z; pars.get_shc();
    [DPsigZdX,DPsigZdY,DPsigZdZ] = spheroidalDP(pars,Xeval,nu_x_spectral,nu_y_spectral,nu_z_spectral);

    %%% Next, let us calculate y_k \sigma_k \boldsymbol{n} + y_k n_k \boldsymbol{\sigma}
    %%% for k = 1, 2, 3.
    y_dot_sig=zeros(size(sigma_x));
    for i=1:ns
        if ~pars.oblate(i)
            Xloc=prolate_spheroid_shape(pars.p,pars.u0(i),pars.a(i));
        else
            Xloc=oblate_spheroid_shape(pars.p,pars.u0(i),pars.a(i));
        end
        y_dot_sig(:,:,i)=sigma_x(:,:,i).*Xloc(:,1)+sigma_y(:,:,i).*Xloc(:,2)+sigma_z(:,:,i).*Xloc(:,3);
    end

    pars.sigma = y_dot_sig; pars.get_shc;
    [ydotsig_dx,ydotsig_dy,ydotsig_dz] = spheroidalDP(pars,Xeval,nu_x_spectral,nu_y_spectral,nu_z_spectral);

    %%% Get normal vectors (at source points) necessary for densities
    % TODO: change this for multiple spheroids
    norm_vecs = get_norm_vecs(pars.p, pars.u0, pars.oblate);
    nx_src = norm_vecs(:,1);
    ny_src = norm_vecs(:,2);
    nz_src = norm_vecs(:,3);

    %%% Now, we must calculate the third quantity, which involves NINE terms.
    % For x-component: \nabla \cdot (S_L[n_x \sigma])
    pars.sigma = nx_src .* sigma_x; pars.get_shc;
    [SPxx, ~, ~] = spheroidalSP(pars, Xeval, nu_x_spectral, nu_y_spectral, nu_z_spectral);
    pars.sigma = nx_src .* sigma_y; pars.get_shc;
    [~, SPxy, ~] = spheroidalSP(pars, Xeval, nu_x_spectral, nu_y_spectral, nu_z_spectral);
    pars.sigma = nx_src .* sigma_z; pars.get_shc;
    [~, ~, SPxz] = spheroidalSP(pars, Xeval, nu_x_spectral, nu_y_spectral, nu_z_spectral);

    % For y-component: \nabla \cdot (S_L[n_y \sigma])
    pars.sigma = ny_src .* sigma_x; pars.get_shc;
    [SPyx, ~, ~] = spheroidalSP(pars, Xeval, nu_x_spectral, nu_y_spectral, nu_z_spectral);
    pars.sigma = ny_src .* sigma_y; pars.get_shc;
    [~, SPyy, ~] = spheroidalSP(pars, Xeval, nu_x_spectral, nu_y_spectral, nu_z_spectral);
    pars.sigma = ny_src .* sigma_z; pars.get_shc;
    [~, ~, SPyz] = spheroidalSP(pars, Xeval, nu_x_spectral, nu_y_spectral, nu_z_spectral);

    % For z-component: \nabla \cdot (S_L[n_z \sigma])
    pars.sigma = nz_src .* sigma_x; pars.get_shc;
    [SPzx, ~, ~] = spheroidalSP(pars, Xeval, nu_x_spectral, nu_y_spectral, nu_z_spectral);
    pars.sigma = nz_src .* sigma_y; pars.get_shc;
    [~, SPzy, ~] = spheroidalSP(pars, Xeval, nu_x_spectral, nu_y_spectral, nu_z_spectral);
    pars.sigma = nz_src .* sigma_z; pars.get_shc;
    [~, ~, SPzz] = spheroidalSP(pars, Xeval, nu_x_spectral, nu_y_spectral, nu_z_spectral);

    %%% Finally, we calculate the Stokes double layer potentials in term of
    %%% the Laplace potentials.
    if isempty(Xeval)
        % Verify that the total number of surface discretization points is split
        % up evenly among the spheroidal bodies.
        [Xloc,~]=pars.get_X();
        np_div=size(Xloc,1)/ns;
        np=round(np_div);
        if abs(np_div-np)>1e-8
            error("\n size of target on surface is not multiple of ns\n")
        end

        second_sum_term_x = SPxx + SPxy + SPxz;
        second_sum_term_y = SPyx + SPyy + SPyz;
        second_sum_term_z = SPzx + SPzy + SPzz;

        if isa(sigma_x,"cell")
            Stk_x=cell(1,ns); Stk_y=Stk_x; Stk_z=Stk_x;

            for i=1:ns % Loop over each body
                Xloc_i=Xloc((i-1)*np+1:i*np,:);

                %%% X term
                % x_k times first coordinate of DP[sigma_k] for k=1,2,3
                first_sum_term = Xloc_i(:,1).*DPsigXdX{i} + Xloc_i(:,2).*DPsigYdX{i} + Xloc_i(:,3).*DPsigZdX{i};
                Stk_x{i}=-first_sum_term + ydotsig_dx{i} + -second_sum_term_x{i};
                
                %%% Y term
                first_sum_term = Xloc_i(:,1).*DPsigXdY{i} + Xloc_i(:,2).*DPsigYdY{i} + Xloc_i(:,3).*DPsigZdY{i};
                Stk_y{i}=-first_sum_term + ydotsig_dy{i} - second_sum_term_y{i};
                
                %%% Z term
                first_sum_term = Xloc_i(:,1).*DPsigXdZ{i} + Xloc_i(:,2).*DPsigYdZ{i} + Xloc_i(:,3).*DPsigZdZ{i};
                Stk_z{i}=-first_sum_term + ydotsig_dz{i} - second_sum_term_z{i};
            end
        else
            Stk_x=zeros(size(sigma_x,1),size(sigma_x,2),ns); Stk_y=Stk_x; Stk_z=Stk_x;
            for i=1:ns
                Xloc_i=Xloc((i-1)*np+1:i*np,:);
                
                %%% X term
                % x_k times first coordinate of DP[sigma_k] for k=1,2,3
                first_sum_term = Xloc_i(:,1).*DPsigXdX(:,:,i) + Xloc_i(:,2).*DPsigYdX(:,:,i) + Xloc_i(:,3).*DPsigZdX(:,:,i);
                Stk_x(:,:,i)=-first_sum_term + ydotsig_dx(:,:,i) - second_sum_term_x(:,:,i);
                
                %%% Y term
                first_sum_term = Xloc_i(:,1).*DPsigXdY(:,:,i) + Xloc_i(:,2).*DPsigYdY(:,:,i) + Xloc_i(:,3).*DPsigZdY(:,:,i);
                Stk_y(:,:,i)=-first_sum_term + ydotsig_dy(:,:,i) - second_sum_term_y(:,:,i);
                
                %%% Z term
                first_sum_term = Xloc_i(:,1).*DPsigXdZ(:,:,i) + Xloc_i(:,2).*DPsigYdZ(:,:,i) + Xloc_i(:,3).*DPsigZdZ(:,:,i);
                Stk_z(:,:,i)=-first_sum_term + ydotsig_dz(:,:,i) - second_sum_term_z(:,:,i);
            end
        end
    else
        % return cell/array of contribution from each spheroid separately.
        if isa(Xeval,"cell")
            Stk_x=cell(1,ns); Stk_y=Stk_x; Stk_z=Stk_x;
            for i=1:ns
                second_sum_term_x = SPxx{i} + SPxy{i} + SPxz{i};
                second_sum_term_y = SPyx{i} + SPyy{i} + SPyz{i};
                second_sum_term_z = SPzx{i} + SPzy{i} + SPzz{i};

                %%% X term
                first_sum_term = Xeval{i}(:,1).*DPsigXdX{i} + Xeval{i}(:,2).*DPsigYdX{i} + Xeval{i}(:,3).*DPsigZdX{i};
                Stk_x{i}=-first_sum_term + ydotsig_dx{i} - second_sum_term_x;
                
                %%% Y term
                first_sum_term = Xeval{i}(:,1).*DPsigXdY{i} + Xeval{i}(:,2).*DPsigYdY{i} + Xeval{i}(:,3).*DPsigZdY{i};
                Stk_y{i}=-first_sum_term + ydotsig_dy{i} - second_sum_term_y;
                
                %%% Z term
                first_sum_term = Xeval{i}(:,1).*DPsigXdZ{i} + Xeval{i}(:,2).*DPsigYdZ{i} + Xeval{i}(:,3).*DPsigZdZ{i};
                Stk_z{i}=-first_sum_term + ydotsig_dz{i} - second_sum_term_z;
            end
        else
            second_sum_term_x = SPxx + SPxy + SPxz;
            second_sum_term_y = SPyx + SPyy + SPyz;
            second_sum_term_z = SPzx + SPzy + SPzz;
            Stk_x=zeros(size(Xeval,1),size(sigma_x,2),ns); Stk_y=Stk_x; Stk_z=Stk_x;
            for i=1:ns
                %%% X term
                first_sum_term = Xeval{i}(:,1).*DPsigXdX(:,:,i) + Xeval{i}(:,2).*DPsigYdX(:,:,i) + Xeval{i}(:,3).*DPsigZdX(:,:,i);
                Stk_x{i}=-first_sum_term + ydotsig_dx(:,:,i) - second_sum_term_x(:,:,i);
                
                %%% Y term
                first_sum_term = Xeval{i}(:,1).*DPsigXdY(:,:,i) + Xeval{i}(:,2).*DPsigYdY(:,:,i) + Xeval{i}(:,3).*DPsigZdY(:,:,i);
                Stk_y{i}=-first_sum_term + ydotsig_dy(:,:,i) - second_sum_term_y(:,:,i);
                
                %%% Z term
                first_sum_term = Xeval{i}(:,1).*DPsigXdZ(:,:,i) + Xeval{i}(:,2).*DPsigYdZ(:,:,i) + Xeval{i}(:,3).*DPsigZdZ(:,:,i);
                Stk_z{i}=-first_sum_term + ydotsig_dz(:,:,i) - second_sum_term_z(:,:,i);
            end
        end
    end
end
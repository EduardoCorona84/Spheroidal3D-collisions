function [Stk_x, Stk_y, Stk_z] = L2StkTLP(X_eval, nu_eval, pars, sigma_x, sigma_y, sigma_z, ns)
    %{
        An implementation of the Laplace to Stokes traction layer potential.

        Inputs
            X_eval      -   target points
            nu_eval     -   normal vector at each target point
            pars        -   parameters needed for to calculate spheroidal Laplace LPs
            sigma_x     -  
            sigma_y     -
            sigma_z     -
            ns          -   number of spheroidal bodies
    %}

    % Flag to enable alternative calculation and avoid the use of the
    % second derivative.
    ALTERNATE_TO_SPP_FLAG = false;

    %%% First, let's generate the necessary vectors for the normal derivatives
    % of the Laplace layer potentials.
    if isempty(X_eval) % Self-evaluation
        % For each body, 
        nu_x_spectral=repmat([1,0,0],size(sigma_x,1),size(sigma_x,2),ns);
        nu_y_spectral=repmat([0,1,0],size(sigma_x,1),size(sigma_x,2),ns);
        nu_z_spectral=repmat([0,0,1],size(sigma_x,1),size(sigma_x,2),ns);
    else % All to targets
        error_msg = "The evaluation points need to have type " + ...
            "'cell' of size 1 x ns. Each cell should contain the target points " + ...
            "for the associated body, and each cell should be of size nt x 3, where " + ...
            "nt is the number of target points on the body.\n" + ...
            "If you have a nt x 3 x ns matrix, then use squeeze(num2cell(Xeval, [1, 2])).' " + ...
            "to reduce back to a 1 x ns cell.";
        assert(isa(X_eval, "cell"), error_msg);
        assert(strcmp(class(X_eval), class(nu_eval)), "Evaluation points and normal vectors at those points should be of the same type.");
        assert(size(X_eval,1) == 1 && size(X_eval, 2) == ns, error_msg)
        nu_x_spectral=cell(1,ns); nu_y_spectral=nu_x_spectral; nu_z_spectral=nu_x_spectral; 
        for i=1:ns
            nu_x_spectral{i} = repmat([1,0,0],size(X_eval{i},1),1);
            nu_y_spectral{i} = repmat([0,1,0],size(X_eval{i},1),1);
            nu_z_spectral{i} = repmat([0,0,1],size(X_eval{i},1),1);
        end
    end

    %%% Extract components of target normal vectors for later
    nx_trg = cell(size(nu_eval)); ny_trg = nx_trg; nz_trg = nx_trg;
    for i=1:numel(nu_eval)
        nu = nu_eval{i};
        nx_trg{i} = nu(:,1);
        ny_trg{i} = nu(:,2);
        nz_trg{i} = nu(:,3);
    end

    %%%%
    %%%% Calculate layer potentials
    %%%%
    %%% First, handle double-derivative (gradient of divergence of vector SLP)
    [graddivSL_sig_U, graddivSL_sig_V, graddivSL_sig_PHI] = spheroidalgraddivSL(pars, sigma_x, sigma_y, sigma_z, X_eval);
 
    sigma = struct();
    sigma.x = sigma_x; sigma.y = sigma_y; sigma.z = sigma_z;

    % Handle y_j * \vec{\sigma} for j = 1,2,3.
    xyz_fields = {'x', 'y', 'z'};
    y_times_sig = struct();
    for j = 1:3 % Index for Xloc
        for k = 1:3 % Index for sigma
            field_name = sprintf('%s%d', xyz_fields{k}, j); % 'x1', 'y2', etc.
            y_times_sig.(field_name) = zeros(size(sigma.x));
        end
    end

    for i=1:ns
        if ~pars.oblate(i)
            Xloc = prolate_spheroid_shape(pars.p, pars.u0(i), pars.a(i));
        else
            Xloc = oblate_spheroid_shape(pars.p, pars.u0(i), pars.a(i));
        end
    
        for j = 1:3
            for k = 1:3
                sig_field = xyz_fields{k};
                output_field = sprintf('%s%d', sig_field, j); % 'x1', 'y2', etc.
                % y_times_sig.z3 means y3*sig_z
                y_times_sig.(output_field)(:,:,i) = sigma.(sig_field)(:,:,i) .* Xloc(:,j);
            end
        end
    end
    
    [graddivSL_y1timessig_U, graddivSL_y1timessig_V, graddivSL_y1timessig_PHI] = spheroidalgraddivSL(pars, y_times_sig.x1, y_times_sig.y1, y_times_sig.z1, X_eval);
    [graddivSL_y2timessig_U, graddivSL_y2timessig_V, graddivSL_y2timessig_PHI] = spheroidalgraddivSL(pars, y_times_sig.x2, y_times_sig.y2, y_times_sig.z2, X_eval);
    [graddivSL_y3timessig_U, graddivSL_y3timessig_V, graddivSL_y3timessig_PHI] = spheroidalgraddivSL(pars, y_times_sig.x3, y_times_sig.y3, y_times_sig.z3, X_eval);

    graddivSL_sph = {
        {graddivSL_y1timessig_U, graddivSL_y1timessig_V, graddivSL_y1timessig_PHI},
        {graddivSL_y2timessig_U, graddivSL_y2timessig_V, graddivSL_y2timessig_PHI},
        {graddivSL_y3timessig_U, graddivSL_y3timessig_V, graddivSL_y3timessig_PHI},
    };

    % Now, need to handle the normal derivatives (this is necessary since the code for S''
    % is for the gradient, and not for the normal derivative).
    graddivSL = struct();
    for j=1:3
        for k=1:3
            density_field = sprintf('y%dsig', j);
            graddivSL.(density_field) = struct();
        end
    end

    for k=1:ns
        S_eval_k = cart2spheroidal(X_eval{k}, pars.a(k), pars.oblate(k));

        [nu_x_sph, ~] = cartNu2spheroidal(nu_x_spectral{k}, S_eval_k, pars.a(k), pars.oblate(k));
        [nu_y_sph, ~] = cartNu2spheroidal(nu_y_spectral{k}, S_eval_k, pars.a(k), pars.oblate(k));
        [nu_z_sph, ~] = cartNu2spheroidal(nu_z_spectral{k}, S_eval_k, pars.a(k), pars.oblate(k));
        nu_sph_vecs = { nu_x_sph, nu_y_sph, nu_z_sph };

        for j = 1:3
            for m = 1:3 % X,Y,Z index
                nu_sph = nu_sph_vecs{m};
                density_field = sprintf('y%dsig', j);

                % etc. graddivSL.y#sig.x
                % Converts (U, V, PHI) to (X, Y, Z)
                % The notation graddivSL_sph{j}{1}{k} means the j-th
                % spheroidal coordinate of graddivSL_y1timessigma for the 
                % k-th particle.
                graddivSL.(density_field).(xyz_fields{m}) = nu_sph(:,1).*graddivSL_sph{j}{1}{k} + nu_sph(:,2).*graddivSL_sph{j}{2}{k} + nu_sph(:,3).*graddivSL_sph{j}{3}{k};
            end
        end

        for m=1:3
            nu_sph = nu_sph_vecs{m}; 
            graddivSL.sig.(xyz_fields{m}) = nu_sph(:,1).*graddivSL_sig_U{k} + nu_sph(:,2).*graddivSL_sig_V{k} + nu_sph(:,3).*graddivSL_sig_PHI{k};
        end
    end

    %%% Now, handle first derivatives.
    %%% In particular, we just want e_k \cdot S^L[\sigma_j] for j=1,2,3 and k = 1,2,3.
    pars.sigma = sigma_x; pars.get_shc;
    [SP_sigmax_X, SP_sigmax_Y, SP_sigmax_Z] = spheroidalSP(pars, X_eval, nu_x_spectral, nu_y_spectral, nu_z_spectral);

    pars.sigma = sigma_y; pars.get_shc;
    [SP_sigmay_X, SP_sigmay_Y, SP_sigmay_Z] = spheroidalSP(pars, X_eval, nu_x_spectral, nu_y_spectral, nu_z_spectral);

    pars.sigma = sigma_z; pars.get_shc;
    [SP_sigmaz_X, SP_sigmaz_Y, SP_sigmaz_Z] = spheroidalSP(pars, X_eval, nu_x_spectral, nu_y_spectral, nu_z_spectral);

    %%% Finally, we calculate the Stokes traction layer potentials in term of
    %%% the Laplace potentials.
    
    if isempty(X_eval)
        error("not implemented.");
        % Verify that the total number of surface discretization points is split
        % up evenly among the spheroidal bodies.
        [Xloc,~]=pars.get_X();
        np_div=size(Xloc,1)/ns;
        np=round(np_div);
        if abs(np_div-np)>1e-8
            error("\n size of target on surface is not multiple of ns\n")
        end

        Stk_x=cell(1,ns); Stk_y=Stk_x; Stk_z=Stk_x;
        for i=1:ns % Loop over each body
            Xloc_i=Xloc((i-1)*np+1:i*np,:);
            n_dot_x = dot(Xloc{i}, nu_eval{i}, 2);

            SP_sum_term_x = nx_trg{i}.*SP_sigmax_X{i} + ny_trg{i}.*SP_sigmax_Y{i} + nz_trg{i}.*SP_sigmax_Z{i};
            SP_sum_term_y = nx_trg{i}.*SP_sigmay_X{i} + ny_trg{i}.*SP_sigmay_Y{i} + nz_trg{i}.*SP_sigmay_Z{i};
            SP_sum_term_z = nx_trg{i}.*SP_sigmaz_X{i} + ny_trg{i}.*SP_sigmaz_Y{i} + nz_trg{i}.*SP_sigmaz_Z{i};

            %%% X term
            xnx_term_X = n_dot_x.*graddivSL.sig.x;
            graddivSL_sum_term_X = nx_trg{i}.*graddivSL.y1sig.x + ny_trg{i}.*graddivSL.y2sig.x + nz_trg{i}.*graddivSL.y3sig.x;
            Stk_x{i} = SP_sum_term_x + graddivSL_sum_term_X - xnx_term_X;
            
            %%% Y term
            xnx_term_Y = n_dot_x.*graddivSL.sig.y;
            graddivSL_sum_term_Y = nx_trg{i}.*graddivSL.y1sig.y + ny_trg{i}.*graddivSL.y2sig.y + nz_trg{i}.*graddivSL.y3sig.y;
            Stk_y{i} = SP_sum_term_y + graddivSL_sum_term_Y - xnx_term_Y;
            
            %%% Z term
            xnx_term_Z = n_dot_x.*graddivSL.sig.z;
            graddivSL_sum_term_Z = nx_trg{i}.*graddivSL.y1sig.z + ny_trg{i}.*graddivSL.y2sig.z + nz_trg{i}.*graddivSL.y3sig.z;
            Stk_z{i} = SP_sum_term_z + graddivSL_sum_term_Z - xnx_term_Z;
        end
    else
        % Returns the cell of contribution from each spheroid separately.
        Stk_x=cell(1,ns); Stk_y=Stk_x; Stk_z=Stk_x;
        for i=1:ns
            n_dot_x = dot(X_eval{i}, nu_eval{i}, 2);

            SP_sum_term_x = nx_trg{i}.*SP_sigmax_X{i} + ny_trg{i}.*SP_sigmax_Y{i} + nz_trg{i}.*SP_sigmax_Z{i};
            SP_sum_term_y = nx_trg{i}.*SP_sigmay_X{i} + ny_trg{i}.*SP_sigmay_Y{i} + nz_trg{i}.*SP_sigmay_Z{i};
            SP_sum_term_z = nx_trg{i}.*SP_sigmaz_X{i} + ny_trg{i}.*SP_sigmaz_Y{i} + nz_trg{i}.*SP_sigmaz_Z{i};

            %%% X term
            xnx_term_X = n_dot_x.*graddivSL.sig.x;
            graddivSL_sum_term_X = nx_trg{i}.*graddivSL.y1sig.x + ny_trg{i}.*graddivSL.y2sig.x + nz_trg{i}.*graddivSL.y3sig.x;
            Stk_x{i} = SP_sum_term_x + graddivSL_sum_term_X - xnx_term_X;
            
            %%% Y term
            xnx_term_Y = n_dot_x.*graddivSL.sig.y;
            graddivSL_sum_term_Y = nx_trg{i}.*graddivSL.y1sig.y + ny_trg{i}.*graddivSL.y2sig.y + nz_trg{i}.*graddivSL.y3sig.y;
            Stk_y{i} = SP_sum_term_y + graddivSL_sum_term_Y - xnx_term_Y;
            
            %%% Z term
            xnx_term_Z = n_dot_x.*graddivSL.sig.z;
            graddivSL_sum_term_Z = nx_trg{i}.*graddivSL.y1sig.z + ny_trg{i}.*graddivSL.y2sig.z + nz_trg{i}.*graddivSL.y3sig.z;
            Stk_z{i} = SP_sum_term_z + graddivSL_sum_term_Z - xnx_term_Z;
        end
    end
end
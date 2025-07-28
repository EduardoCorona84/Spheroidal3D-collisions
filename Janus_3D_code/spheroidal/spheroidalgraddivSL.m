function [graddivSL_x, graddivSL_y, graddivSL_z,varargout]=spheroidalgraddivSL(params, sigma_x, sigma_y, sigma_z, X)
    %--------------------------------------------------------------------%
    % spheroidalgraddivSL computes the gradient of the divergence of the Laplace vector SLP.
    % Note that the input density should be a vector, and the output will be a vector.
    %    (o) params:
    %       (o) p = order
    %       (o) u0 = 1/eccentricity of the spheroid surface.
    %       (o) a = scale of spheroid
    %       (o) isReal = Can be set if real output is expected. then imaginary 
    %           parts are removed 
    %       (o) sigma_coefficients: spherical harmonic coefficients of sigma
    %    (o) (sigma_x, sigma_y, sigma_z) = each sigma_j vector denotes the vector
    %        density on the surface, all as a function of (theta, phi).
    %    (o) X = (x,y,z) coordinates of target points. Is a cell of length
    %        ns or an nt x 3 x ns matrix.
    %    (o) nu = (nu_x,nu_y,nu_z) coordinates of normals at target points. Is
    %    a cell fo length ns or an nt x 3 x ns matrix.
    %    (o) nu_2,... = each input same format as nu, and prompts an additional
    %    entry to the output. If nargin==3, then nargout==1.
    %    
    %   Returns
    %--------------------------------------------------------------------%
    
    DEBUG_FLAG = true;
    
    % Pre-processing
    %--------------------------------------------------------------------%
    if isempty(params.u0)
        error("No surface parameter u_0 given")
    end
    if isempty(params.a)
        error("No surface parameter 'a' given")
    end
    
    %%% Calculate the sigma coefficients for the vector density
    %%% A little hacky, but gets the job done...
    params.sigma = sigma_x; params.get_shc; shc_x = params.sigma_coefficients;
    params.sigma = sigma_y; params.get_shc; shc_y = params.sigma_coefficients;
    params.sigma = sigma_z; params.get_shc; shc_z = params.sigma_coefficients;
    [sp,nf,ns]=size(params.sigma_coefficients);

    p=params.p;
    u0=params.u0;
    a=params.a;
    isReal=params.isReal;
    oblate = params.oblate;

    %%% Now, we need to do a basis transformation since we use
    %%% the spectral representation for the Laplace scalar SLP.
    [Gshc_x, Gshc_y, Gshc_z] = SL_basis_transform(p, sp, u0, shc_x, shc_y, shc_z, oblate);
    
    if length(u0)==1
        u0=u0*ones(1,ns); 
    end
    if length(a)==1
        a=a*ones(1,ns); 
    end
    
    % No target points given; do self evaluation
    % -----------------------------------------------------------------------
    if nargin==1
        error("not implemented.");
        %if X not given, calculate on-surface points; nu taken to be surface
        %normal
    
        [theta_x,phi_x]=gl_grid(p);
        v_x=cos(theta_x);
        nt=length(theta_x);
        
        Y=zeros(nt,sp);
        for n=0:p  %loop over terms in spheroidal harmonic expansion
            Yn = Ynm(n,[],acos(v_x)',phi_x);
            Y(:,n^2+1:(n+1)^2)=Yn;
        end
    
        % Calculate spectra on surface with outward normal
        [spectra_surf,~,~]=DPspectrum(p,u0,params.a,oblate);
    
        %reshape so that each column of DP_coefs corresponds to each function 
        % on each spheroid
        spectra_surf = reshape(spectra_surf,sp,1,ns);
        spectra_matrix = repmat(spectra_surf,1,nf,1);
        DP_coefs = spectra_matrix .* shc; 
        DP_coefs = reshape(DP_coefs, sp, [], 1);
    
        %multiply coefficients with spheroidal harmonics
        DP=Y*DP_coefs;
    
        % reshape to match original shape of sigma
        DP = reshape(DP,[],nf,ns);
    
        %%%%%%%%%%%%
        % assuming nf=1
        coef_mat=zeros(np,1,ns);
        for i=1:ns
            coef_mat(:,:,i)=oblate(i).*1./sqrt(u0(i)^2+v_x.^2)+~oblate(i).*1./sqrt(u0(i)^2-v_x.^2);
        end
        DP = coef_mat.*DP;    
        %%%%%%%%%%%%%%%%
    
        if isReal
            DP=real(DP);
        end
    %  Self Eval, arbitrary normal.
    % ------------------------------------------------------------------------
    elseif isempty(X)
        error("not implemented.");
        % When we want to get surface dS/dnu with arbitrary normal,
        % to avoid error in converting coordinates, enter X=[] with nu vectors
        % to perform on-surface calculations.
    
        if isa(nu,"cell")
            nu_t=nu;
        else
            nu_t=mat2cell(nu,size(nu,1),size(nu,2),ones(1,size(nu,3)));
            nu_t=reshape(nu_t,1,length(nu_t));
        end
    
        DP=cell(1,ns);
    
        [theta2,phi_k]=gl_grid(p);
        v_k=cos(theta2);
        nt_r=length(v_k);
        Yr=zeros(nt_r*2,sp);
    
        for n=0:p  %loop over terms in spheroidal harmonic expansion
            Yn=Ynm(n,[],real(acos(v_k))',phi_k); % v_x_r exceeds [-1,1] by 1e-8, but acos() returns imaginary values. Impose real values for Ynm.
            Yr(1:nt_r,n^2+1:(n+1)^2)=Yn;

            Yn1=Ynm(n+1,-n:n,real(acos(v_k))',phi_k); 
            yn1_scale = sqrt((2*n+1)/(2*n+3).*(n+(-n:n)+1)./(n-(-n:n)+1));

            Yr(nt_r+1:end,n^2+1:(n+1)^2) = yn1_scale.*Yn1;
        end
        %%%%%%%%%%%%%%%%%%%%%

        for k=1:ns  %loop over each spheroid surface we want to evaluate
            [~,Xself_k]=params.get_X(k);
            S=cart2spheroidal(Xself_k,a(k),oblate(k));
    
            u_k=u0(k).*ones(size(v_k));
    
            nu_list=cell(1,nvarin+1); nu_sph_list=nu_list;
            nu_list{1}=nu_t{k};
            nu_sph_list{1}=cartNu2spheroidal(nu_t{k},S,a(k),oblate(k));
            for nu_ind=1:nvarin
                nu_list{nu_ind+1}=nu_extra_cells{nu_ind}{k};
                [nu_sph_temp,~]=cartNu2spheroidal(nu_list{nu_ind+1},S,a(k),oblate(k));
                nu_sph_list{nu_ind+1}=nu_sph_temp;
            end
    
            for nu_ind=1:nvarin+1
                [spectra_nm_prime,spectra_nm,spectra_n1m] = graddivSL_away(p,u0(k),a,u_k,v_k,nu_sph_list{nu_ind},oblate(k));
                
                FYr = (spectra_nm_prime + spectra_nm).*Yr(1:nt_r,:) + spectra_n1m.*Yr(nt_r+1:end,:);
                DP_k=FYr*shc(:,:,k);
    
                if isReal
                    DP_k=real(DP_k);
                end
                
                if nu_ind==1
                    DP{k}=DP_k;
                else
                    varargout{nu_ind-1}{k}=DP_k;
                end
            end
    
        end
    
        if ~isa(nu, "cell")
            DP = cell2mat(reshape(DP,1,1,ns));
        end
    elseif nargin > 1 % Target points given
        %%% Input processing and input validation
        if isa(X, "cell")
            Xt=X;
        else
            Xt = mat2cell(X,size(X,1),size(X,2),ones(1,size(X,3)));
            Xt = reshape(Xt,1,length(Xt));
        end
        
        if length(Xt) ~= ns
            error("dimensions of target input and number of surfaces do not match. X should be a 1 x N cell or nt x 3 x N matrix.")
        end
    
        %%% Actual computation
        % Calculate spectra for interior, exterior, and surface
        for k=1:ns % Loop over each spheroid surface
            
            Xt_k = Xt{k};
            [nt_k,d] = size(Xt_k);
            DP_k_x = zeros(nt_k, nf); DP_k_y = zeros(nt_k, nf); DP_k_z = zeros(nt_k, nf);
    
            % Grab the sigma coefficients associated with the spheroid
            Gshc_x_k = Gshc_x(:,:,k);
            Gshc_y_k = Gshc_y(:,:,k);
            Gshc_z_k = Gshc_z(:,:,k);

            if nt_k > 0
                if d~=3
                    error("Dimensions of target point array should be N x 3.")
                end

                % Convert targets to spheroidal coords
                S=cart2spheroidal(Xt_k,a(k),oblate(k));
                u_x=S(:,1);
                
                % Split up interior/surface/exterior
                S_int=S(u_x < u0(k),:);
                S_surf=S(u_x == u0(k),:);
                S_ext=S(u_x > u0(k),:);

                % If there's no points in a particular region, don't do any
                % computations.
                do_int=~isempty(S_int);
                do_surf=~isempty(S_surf); 
                do_ext=~isempty(S_ext);
            
                % Split coordinates into 3 regions.
                regions=[do_int,do_surf,do_ext];
                Sregions={S_int,S_surf,S_ext};

                % For each of the three regions, we need to take care of
                % the x, y, and z components.
                DPregions=cell(3,3,1);
                
                % Loop over each region
                for r=1:3
                    if regions(r)
                        Sr=Sregions{r};
                        u_x_r=Sr(:,1);
                        v_x_r=Sr(:,2);
                        phi_x_r=Sr(:,3);
                        v_x_r_real=real(v_x_r);
                        if abs(v_x_r_real-v_x_r)>1e-10
                            fprintf("v_x_r imaginary\n ")
                        end
                        v_x_r=v_x_r_real;

                        [Ucomponent, Vcomponent, PHIcomponent] = graddivSL_away(p,u0(k),a,u_x_r,v_x_r,phi_x_r,Gshc_x_k,Gshc_y_k,Gshc_z_k,oblate(k));
                        
                        % Now, we are done, so store the information for the associated region.
                        DPregions{r, 1} = Ucomponent;
                        DPregions{r, 2} = Vcomponent;
                        DPregions{r, 3} = PHIcomponent;
                    end
                end
            
                % Recombine all DPs from interior/exterior/surface
                DP_k_x(u_x < u0(k),:) = DPregions{1,1};
                DP_k_x(u_x == u0(k),:) = DPregions{2,1};
                DP_k_x(u_x > u0(k),:) = DPregions{3,1};

                DP_k_y(u_x < u0(k),:) = DPregions{1,2};
                DP_k_y(u_x == u0(k),:) = DPregions{2,2};
                DP_k_y(u_x > u0(k),:) = DPregions{3,2};

                DP_k_z(u_x < u0(k),:) = DPregions{1,3};
                DP_k_z(u_x == u0(k),:) = DPregions{2,3};
                DP_k_z(u_x > u0(k),:) = DPregions{3,3};
        
                if isReal
                    DP_k_x = real(DP_k_x);
                    DP_k_y = real(DP_k_y);
                    DP_k_z = real(DP_k_z);
                end
            end
            
            %%% Evaluation for particle k
            graddivSL_x{k} = DP_k_x;
            graddivSL_y{k} = DP_k_y;
            graddivSL_z{k} = DP_k_z;
        end

        if ~isa(X, "cell")
            graddivSL_x = cell2mat(reshape(graddivSL_x,1,1,ns));
            graddivSL_y = cell2mat(reshape(graddivSL_y,1,1,ns));
            graddivSL_z = cell2mat(reshape(graddivSL_z,1,1,ns));
        end
    end
end
    
    
function [lambda_nm_prime, lambda_nm, lambda_n1m] = DPspectrum(p, u0, a, oblate)
    %{
        Calculates the coefficients for D'^-, D'^+, and D' on the surface of the spheroid. 
    %}
    sp = (p+1)^2;
    ii = (1:sp)';
    nn = floor(sqrt(ii-1));
    mm = ii - nn.^2 - nn - 1;
    anm_base = factorial(nn-mm)./factorial(nn+mm) .* ((-1) .^ (mm));

    L = [];
    if ~oblate
        anm = anm_base .* -1 .* (u0.^2 + 1); % cnm
        L = legendre_otc(p,1j.*u0,1,1,1);
    else
        error("Not implemented.");
        anm = anm_base .* (u0.^2 - 1); % bnm
        L = legendre_otc(p,u0,1,1,1);
    end

    P = L{1}; Q = L{2}; dP = L{3}; dQ = L{4};

    % Coefficient of D' of surface
    lambda_nm_prime = anm.' .* sqrt(u0.^2 - 1) .* dP.' .* dQ.' ./ a;
    lambda_nm = anm .* ((dQ .* P + dP .* Q) ./ 2) ./ a; % Is this right?
    lambda_n1m = lambda_nm;
end
     
function [Ucomponent, Vcomponent, PHIcomponent] = graddivSL_away(p,u0,a,u,v,phi,Gshc_x,Gshc_y,Gshc_z,oblate)
    %{
        Obtains the coefficients associated with Y_k^m for k = n, n+1, n+2
        Note that unlike all of the similar functions (e.g. spheroidalSP, etc.), this
        handles ALL coefficients associated with the spheroidal harmonics.

        Inputs
            - p
            - u0
            - a
            - u_x
            - v_x
            - Gshc_x, Gshc_y, Gshc_z
            - oblate
    %}
    sp=(p+1)^2;
    ii = (1:sp)'; n = floor(sqrt(ii-1)); m=ii-n.^2-n-1;
    nt_r=length(u);

    % Necessary for looping later
    gshc_types = {'gshcx', 'gshcy', 'gshcz'};
    Gshc_data = {Gshc_x, Gshc_y, Gshc_z}; 

    % Note that the factor of 3 here allows us to store
    % Y_n, Y_{n+1}, and Y_{n+2} in the same vector.
    Yr=zeros(nt_r*3,sp);

    % Note that unlike the other spheroidal code, the coefficients associated
    % with the Ynm's are handled in the spectral coefficient code.
    for n=0:p  % Loop over terms in spheroidal harmonic expansion
        % v_x_r exceeds [-1,1] by 1e-8, but acos() returns imaginary values. Impose real values for Ynm.
        %%% Store Yn harmonics
        Yn=Ynm(n, [], real(acos(v))', phi);
        Yr(1:nt_r, n^2+1:(n+1)^2)=Yn;

        %%% Store Yn+1 harmonics
        Yn1=Ynm(n+1, -n:n, real(acos(v))', phi);
        Yr(nt_r+1:2*nt_r, n^2+1:(n+1)^2) = Yn1;

        %%% Store Yn+2 harmonics
        Yn2=Ynm(n+2, -n:n, real(acos(v))', phi); 
        Yr(2*nt_r+1:end, n^2+1:(n+1)^2) = Yn2;
    end

    %%% Calculate fnm and fnm' and fnm''
    if oblate
        error("not implemented.");
    else
        bnm = a*factorial(n-m)./factorial(n+m) .* ((-1) .^ (m)) .* sqrt(u0.^2-1);
        L = legendre_otc(p,u0,1,1,1);
        if abs(u)-u0 < 1e-14 % interior
            gnm = L{4}; % Q'(u_0)
        elseif abs(u)-u0 > 1e-14 % exterior
            gnm = L{3}; % P'(u_0)
        end
        [Fr, Fp, Fpp] = solid_harmonic_prime(p, u0, u, oblate);
        common_coeffs = a.^(-2).*bnm.*gnm;
    end

    %%% Define helper function
    function gshc_coeff = calculate_gshc_coeff(coeffs_for_gshc)
        %{
            Calculates the coefficient associated with a GSHC given the
            coefficients for the Associated Legendre function (i.e. f, f', 
            and f'').
        %}
        get_Y_coeff = @(f_coeffs) f_coeffs{1}.*Fr + f_coeffs{2}.*Fp + f_coeffs{3}.*Fpp;
    
        % Grab f_coeffs for each Y type
        f_nm = coeffs_for_gshc.Ynm;
        f_n1m = coeffs_for_gshc.Yn1m;
        f_n2m = coeffs_for_gshc.Yn2m;
    
        % Calculate the coefficient for each Ynm harmonic
        Ynm_coeff  = get_Y_coeff(f_nm);
        Yn1m_coeff = get_Y_coeff(f_n1m);
        Yn2m_coeff = get_Y_coeff(f_n2m);
    
        gshc_coeff = Ynm_coeff.*Yr(1:nt_r,:) + Yn1m_coeff.*Yr(nt_r+1:2*nt_r,:) ...
                     + Yn2m_coeff.*Yr(2*nt_r+1:end,:);
    end

    if ~oblate
        if norm(abs(u)-u0)<1e-14 % on surface with arbitrary nu
            error("on surface not implemented.");
        else % off-surface.
            n = n'; m = m';
            %%% Note that below is very granually split so that it's easier
            %%% to debug/fix size issues (hence, the ugliness).
            % Working backwards, the procedure is this (so step 3 is first 
            % and step 1 is last):
            % 1. Handle each sigma coefficient
            % 2. Within each sigma coefficient, handle the Ynm coefficients
            % 3. Within each Ynm coefficient, handle the f coefficients
            coeffs = spheroidalgraddivSLcoefficients(u, v, n, m, oblate);

            %%% --- U COMPONENT ---
            Ucomponent = 0;
            for i = 1:3
                type = gshc_types{i};
                Gshc_coeff = calculate_gshc_coeff(coeffs.U.(type));
                Ucomponent = Ucomponent + (common_coeffs' .* Gshc_coeff) * Gshc_data{i};
            end

            %%% --- V COMPONENT ---
            Vcomponent = 0;
            for i = 1:3
                type = gshc_types{i};
                Gshc_coeff = calculate_gshc_coeff(coeffs.V.(type));
                Vcomponent = Vcomponent + (common_coeffs' .* Gshc_coeff) * Gshc_data{i};
            end

            %%% --- PHI COMPONENT ---
            PHIcomponent = 0;
            for i = 1:3
                type = gshc_types{i};
                Gshc_coeff = calculate_gshc_coeff(coeffs.PHI.(type));
                PHIcomponent = PHIcomponent + (common_coeffs' .* Gshc_coeff) * Gshc_data{i};
            end
        end
    else
        error("oblate not implemented.");
    end
end

function [Fr, Fp, Fpp]=solid_harmonic_prime(p, u0, u_x, oblate)
    %{
        Solid spheroidal harmonics can be written as f_n^m(u)Y_n^m(v, phi).
        This function returns what Fr = f_n^m is (and its derivative as Fp), depending 
        on whether we are in the exterior or interior.
    %}
    Fr = ones(size(u_x,1),(p+1)^2);
    Fp = ones(size(u_x,1),(p+1)^2);
    Fpp = ones(size(u_x,1),(p+1)^2);

    if oblate
        u_x = 1j.*u_x;
    end

    if abs(u_x)-u0 < -1e-14 % Interior
        PQ=legendre_otc(p,u_x,1,2,2);
        P=PQ{1}; dP=PQ{3}; ddP = PQ{5};
        Fr=P.'; Fp=dP.'; Fpp=ddP.';
    elseif abs(u_x)-u0 > 1e-14 % Exterior
        PQ=legendre_otc(p,u_x,1,2,2);
        Q=PQ{2}; dQ=PQ{4}; ddQ=PQ{6};
        Fr=Q.'; Fp=dQ.'; Fpp=ddQ.';
    end
end

function [Gshc_x, Gshc_y, Gshc_z] = SL_basis_transform(p, sp, u0, shc_x, shc_y, shc_z, oblate)
    %{
        Do basis transformations for all possible u0 values we might have.
        For reference, see spheroidalSL.m and spheroidalSP.m.
    %}
    all_u0s = unique(u0);
    Gshc_x = zeros(size(shc_x)); Gshc_y = zeros(size(shc_x)); Gshc_z = zeros(size(shc_x));

    % Change to basis that diagonalizes the single layer operator
    for i=1:length(all_u0s)
        this_u0=all_u0s(i);

        %Calculate basis transformation matrix G: Ynm -> Ynm / sqrt(u0^2-v^2)
        %or /sqrt(u0^2+v^2) for oblate
        this_G_pro = Gmatrix(p,this_u0,0,0);
        this_G_obl = Gmatrix(p,this_u0,0,1);

        % Keep track of indices and shc of particles with the same u0
        this_index = find(u0==this_u0); 
        obl_index = nonzeros(this_index.*oblate(this_index));
        pro_index = nonzeros(this_index.*~oblate(this_index));

        %%% X
        obl_shc_x = reshape(shc_x(:,:,obl_index),sp,[],1);
        pro_shc_x = reshape(shc_x(:,:,pro_index),sp,[],1);

        % Apply basis transformation to spheroidal harmonic coefficients
        obl_Gshc_x = this_G_obl\obl_shc_x;
        pro_Gshc_x = this_G_pro\pro_shc_x;
        
        Gshc_x(:,:,obl_index) = reshape(obl_Gshc_x,sp,[],length(obl_index));
        Gshc_x(:,:,pro_index) = reshape(pro_Gshc_x,sp,[],length(pro_index));

        %%% Y
        obl_shc_y = reshape(shc_y(:,:,obl_index),sp,[],1);
        pro_shc_y = reshape(shc_y(:,:,pro_index),sp,[],1);

        % Apply basis transformation to spheroidal harmonic coefficients
        obl_Gshc_y = this_G_obl\obl_shc_y;
        pro_Gshc_y = this_G_pro\pro_shc_y;
        
        Gshc_y(:,:,obl_index) = reshape(obl_Gshc_y,sp,[],length(obl_index));
        Gshc_y(:,:,pro_index) = reshape(pro_Gshc_y,sp,[],length(pro_index));

        %%% Z
        obl_shc_z = reshape(shc_z(:,:,obl_index),sp,[],1);
        pro_shc_z = reshape(shc_z(:,:,pro_index),sp,[],1);

        % Apply basis transformation to spheroidal harmonic coefficients
        obl_Gshc_z = this_G_obl\obl_shc_z;
        pro_Gshc_z = this_G_pro\pro_shc_z;
        
        Gshc_z(:,:,obl_index) = reshape(obl_Gshc_z,sp,[],length(obl_index));
        Gshc_z(:,:,pro_index) = reshape(pro_Gshc_z,sp,[],length(pro_index));
    end
end
function [DP,varargout]=spheroidalDP(params,X,nu,varargin)
    %--------------------------------------------------------------------%
    % spheroidalDP computes the normal derivative of the Laplace double layer 
    % potential of a density 'sigma' at a target point.
    %    (o) params:
    %       (o) sigma = density on the spheroid surface, as a function of (theta,phi). 
    %           See gl_grid. cos(theta) are gauss-Legendre nodes and phi are 
    %           equispaced. sigma is size [np,nf,ns] and each column of sigma 
    %           corresponds to a different spheroid surface.
    %       (o) p = order
    %       (o) u0 = 1/eccentricity of the spheroid surface.
    %       (o) a = scale of spheroid
    %       (o) isReal = Can be set if real output is expected. then imaginary 
    %           parts are removed 
    %       (o) sigma_coefficients: spherical harmonic coefficients of sigma
    %    (o) X = (x,y,z) coordinates of target points. Is a cell of length
    %        ns or an nt x 3 x ns matrix.
    %    (o) nu = (nu_x,nu_y,nu_z) coordinates of normals at target points. Is
    %    a cell fo length ns or an nt x 3 x ns matrix.
    %    (o) nu_2,... = each input same format as nu, and prompts an additional
    %    entry to the output. If nargin==3, then nargout==1.
    %    
    % Returns DP as a matrix if X is a matrix or not provided, or as a cell if
    % X is a cell.
    % 
    %--------------------------------------------------------------------%
    
    
    % Pre-processing
    %--------------------------------------------------------------------%
    if isempty(params.sigma) && isempty(params.sigma_coefficients)
        error("No surface density and no spheroidal harmonic coefficients given. Must provide at least one.")
    end
    if isempty(params.u0)
        error("No surface parameter u_0 given")
    end
    if isempty(params.a)
        error("No surface parameter 'a' given")
    end
    
    p=params.p;
    u0=params.u0;
    a=params.a;
    isReal=params.isReal;
    sigma=params.sigma;
    oblate = params.oblate;
    
    % sp - Number of spheroidal harmonic coefficients
    % nf - 
    % ns - Number of spheroidal bodies
    [sp,nf,ns]=size(params.sigma_coefficients);
    [np,nf2,ns2]=size(sigma);
    if nargin <3
        nu = repmat([1,0,0],1,1,ns);
    end
    
    if nargin>3
        if (nargout-1)~=(nargin-3)
            error("Extra nu given and expecting extra output have mismatch in size.")
        end
    end
    
    nvarin=nargin-3;
    more_nu=nvarin>0;
    if more_nu
        % if extra nu inputted, store them in cell nu_extra to be used later. 
        nu_extra=cell(1,nvarin);
        for nind=1:nvarin
            nu_extra{nind}=varargin{nind};
        end
    end
    
    if sp ~= (p+1)^2 || (nf~=nf2 && ~isempty(sigma)) || (ns2~=ns && ~isempty(sigma))
        fprintf("Updating spheroidal harmonic coefficients.\n")
        params.get_shc;
    end
    
    shc = params.sigma_coefficients;
    [sp,nf,ns]=size(params.sigma_coefficients);
    
    if length(u0)==1
        u0=u0*ones(1,ns); 
    end
    if length(a)==1
        a=a*ones(1,ns); 
    end
    
    % No target points given; do self evaluation
    % -----------------------------------------------------------------------
    if nargin==1
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
        % Note that this is not the basis transformation that was done in
        % spheroidalSP. Instead, this just accounts for the missing scale
        % factor in the spectrum.
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
        % When we want to get surface dS/dnu with arbitrary normal,
        % to avoid error in converting coordinates, enter X=[] with nu vectors
        % to perform on-surface calculations.
    
        if isa(nu,"cell")
            nu_t=nu;
        else
            nu_t=mat2cell(nu,size(nu,1),size(nu,2),ones(1,size(nu,3)));
            nu_t=reshape(nu_t,1,length(nu_t));
        end
    
        if more_nu
            nu_extra_cells=cell(1,nvarin);
            for nu_ind=1:nvarin
                nu_temp=nu_extra{nu_ind};
                if isa(nu_temp,"cell")
                    nu_extra_cells{nu_ind}=nu_temp;
                else
                    nu_temp_cell=mat2cell(nu_temp,size(nu_temp,1),size(nu_temp,2),ones(1,size(nu_temp,3)));
                    nu_extra_cells{nu_ind}=reshape(nu_temp_cell,1,length(nu_temp_cell));
                end
            end
        end
    
        DP=cell(1,ns);
    
        if more_nu
            varargout=cell(1,nvarin);
            for nu_ind=1:nvarin
                varargout{nu_ind}=cell(1,ns);
            end
        end
    

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
                [spectra_nm_prime,spectra_nm,spectra_n1m] = DPspectrum_away(p,u0(k),a(k),u_k,v_k,nu_sph_list{nu_ind},oblate(k));
                
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
    
            if more_nu
                for nu_ind=1:nvarin
                    varargout{nu_ind}=cell2mat(reshape(varargout{nu_ind},1,1,ns));
                end
            end
        end
    elseif nargin > 1 % Target points given
        if isa(X, "cell")
            Xt=X;
        else
            Xt = mat2cell(X,size(X,1),size(X,2),ones(1,size(X,3)));
            Xt = reshape(Xt,1,length(Xt));
        end
    
        if isa(nu,"cell")
            nu_t=nu;
        else
            nu_t=mat2cell(nu,size(nu,1),size(nu,2),ones(1,size(nu,3)));
            nu_t=reshape(nu_t,1,length(nu_t));
        end
    
        if more_nu
            nu_extra_cells=cell(1,nvarin);
            for nu_ind=1:nvarin
                nu_temp=nu_extra{nu_ind};
                if isa(nu_temp,"cell")
                    nu_extra_cells{nu_ind}=nu_temp;
                else
                    nu_temp_cell=mat2cell(nu_temp,size(nu_temp,1),size(nu_temp,2),ones(1,size(nu_temp,3)));
                    nu_extra_cells{nu_ind}=reshape(nu_temp_cell,1,length(nu_temp_cell));
                end
            end
        end
        
        if length(Xt) ~= ns
            error("dimensions of target input and number of surfaces do not match. X should be a 1 x N cell or nt x 3 x N matrix.")
        end
    
        DP=cell(1,ns);
    
        if more_nu
            varargout=cell(1,nvarin);
            for nu_ind=1:nvarin
                varargout{nu_ind}=cell(1,ns);
            end
        end
    
        % % Starting the actual computation
        %--------------------------------------------------------------------%
    
        % Calculate spectra for interior, exterior, and surface
        for k=1:ns  %loop over each spheroid surface we want to evaluate
            
            Xtk = Xt{k};
            nu_k=nu_t{k};
            [ntk,d]=size(Xtk);
            DPk = zeros(ntk,nf);
            if more_nu
                DP_extra_at_k=zeros(ntk,nf,nvarin);
                DP_extra_at_k=mat2cell(DP_extra_at_k,ntk,nf,ones(1,nvarin));
            end
    
            if ntk > 0
                if d~=3
                    error("Dimensions of target point array should be N x 3.")
                end

                % Convert targets to spheroidal coords
                S=cart2spheroidal(Xtk,a(k),oblate(k));

                u_x=S(:,1);
                if min(u_x)-1<1e-12 && oblate(k)==0
                    fprintf("\n warning: u_x~1. check spheroid #%d\n ",k);
                end

                %%% split up interior/surface/exterior
                indices_interior = (u_x < u0(k) - 1e-14);
                indices_surface = (abs(u_x-u0(k)) <= 1e-14);
                indices_exterior = (u_x > u0(k) + 1e-14);

                S_int = S(indices_interior,:);
                S_surf = S(indices_surface,:);
                S_ext = S(indices_exterior,:);

                Nu_int = nu_k(indices_interior,:);
                Nu_surf = nu_k(indices_surface,:);
                Nu_ext = nu_k(indices_exterior,:);
    
                if more_nu
                    Nu_int_extra=cell(1,nvarin);
                    Nu_surf_extra=cell(1,nvarin);
                    Nu_ext_extra=cell(1,nvarin);
                    for nu_ind=1:nvarin
                        Nu_int_extra{nu_ind}=nu_extra_cells{nu_ind}{k}(u_x<u0(k),:);
                        Nu_surf_extra{nu_ind}=nu_extra_cells{nu_ind}{k}(u_x==u0(k),:);
                        Nu_ext_extra{nu_ind}=nu_extra_cells{nu_ind}{k}(u_x>u0(k),:);
                    end
                end
                
                % If there's no points in a particular region, don't do any
                % computations
                do_int=~isempty(S_int);
                do_surf=~isempty(S_surf); 
                do_ext=~isempty(S_ext);
            
                % Split coordinates into 3 regions: We have a DP for each region
                regions=[do_int,do_surf,do_ext];
                Sregions={S_int,S_surf,S_ext};
                Nuregions={Nu_int,Nu_surf,Nu_ext};
                DPregions=cell(3,1);
                if more_nu
                    Nuregions_extra={Nu_int_extra,Nu_surf_extra,Nu_ext_extra};
                    DPregions_extra={cell(1,nvarin),cell(1,nvarin),cell(1,nvarin)}; %DPregions_extra{1}={DP_int from nu_1},{DP_int from nu_2},...
                end
                
                % Loop over each region
                for r=1:3
                    if regions(r)
                        Sr=Sregions{r};
                        u_x_r=Sr(:,1);
                        v_x_r=Sr(:,2);
                        phi_x_r=Sr(:,3);
                        nu_r_cart=Nuregions{r};
                        [nu_r_sph,~] = cartNu2spheroidal(nu_r_cart,Sr,a(k),oblate(k));
                        v_x_r_real=real(v_x_r);
                        if abs(v_x_r_real-v_x_r)>1e-10
                            fprintf("v_x_r imaginary\n ")
                        end
                        v_x_r=v_x_r_real;
    
                        [spectra_nm_prime,spectra_nm,spectra_n1m] = DPspectrum_away(p,u0(k),a(k),u_x_r,v_x_r,nu_r_sph,oblate(k));
                
                        [Fr,Fp]=solid_harmonic_prime(p,u0(k),u_x_r,oblate(k));
                        nt_r=length(u_x_r);
                        Yr=zeros(nt_r*2,sp);
                
                        for n=0:p  %loop over terms in spheroidal harmonic expansion
                            % Yn=Ynm(n,[],acos(v_x_r)',phi_x_r);
                            Yn=Ynm(n,[],real(acos(v_x_r))',phi_x_r); % v_x_r exceeds [-1,1] by 1e-8, but acos() returns imaginary values. Impose real values for Ynm.
                            Yr(1:nt_r,n^2+1:(n+1)^2)=Yn;
                            % Yn1 = Ynm(n+1,-n:n,acos(v_x_r)',phi_x_r);
                            Yn1=Ynm(n+1,-n:n,real(acos(v_x_r))',phi_x_r); 
                            yn1_scale=sqrt((2*n+1)/(2*n+3).*(n+(-n:n)+1)./(n-(-n:n)+1));
                            Yr(nt_r+1:end,n^2+1:(n+1)^2) = yn1_scale.*Yn1;
                        end
    
                        % Solid spheroidal harmonics
                        FYr = (spectra_nm_prime.*Fp + spectra_nm.*Fr).*Yr(1:nt_r,:)...
                            +spectra_n1m.*Fr.*Yr(nt_r+1:end,:);
    
                        % Spheroidal harmonic coefficients of S'[fac*sigma]
                        DPregions{r}=FYr*shc(:,:,k);
    
                        % if there are more nu vectors to calculate DP with,
                        % use already computed Fr,Yr at target points, and
                        % evaluate spectrum using new normal
                        if more_nu
                            nu_r_cell=Nuregions_extra{r};
                            DP_r_cell=cell(1,nvarin);
                            for nu_ind=1:nvarin
                                nu_r_cart=nu_r_cell{nu_ind};
                                [nu_r_sph,~] = cartNu2spheroidal(nu_r_cart,Sr,a(k),oblate(k));
                                [spectra_nm_prime,spectra_nm,spectra_n1m] = DPspectrum_away(p,u0(k),a(k),u_x_r,v_x_r,nu_r_sph,oblate(k));
                                FYr = (spectra_nm_prime.*Fp + spectra_nm.*Fr).*Yr(1:nt_r,:)+spectra_n1m.*Fr.*Yr(nt_r+1:end,:);
                                DP_r_cell{nu_ind}=FYr*shc(:,:,k);
                            end
                            DPregions_extra{r}=DP_r_cell;
                        end
    
                    end
                end
            
                % Recombine all DPs from interior/exterior/surface
                DPk(indices_interior,:) = DPregions{1};
                DPk(indices_surface,:) = DPregions{2};
                DPk(indices_exterior,:) = DPregions{3};
    
                if more_nu
                    for nu_ind=1:nvarin
                        DPk_temp=DP_extra_at_k{nu_ind};
                        DPk_temp(indices_interior,:) = DPregions_extra{1}{nu_ind};
                        DPk_temp(indices_surface,:) = DPregions_extra{2}{nu_ind};
                        DPk_temp(indices_exterior,:) = DPregions_extra{3}{nu_ind};
                        DP_extra_at_k{nu_ind}=DPk_temp;
                        if isReal
                            DP_extra_at_k{nu_ind}=real(DPk_temp);
                        end
                    end
                end
        
                if isReal
                    DPk=real(DPk);
                end
            end
            
            %evluation for particle k
            DP{k} = DPk;
    
            if more_nu
                for nu_ind=1:nvarin
                    varargout{nu_ind}{k}=DP_extra_at_k{nu_ind};
                end
            end
                    
    
        end
    
        if ~isa(X, "cell")
            DP = cell2mat(reshape(DP,1,1,ns));
    
            if more_nu
                for nu_ind=1:nvarin
                    varargout{nu_ind}=cell2mat(reshape(varargout{nu_ind},1,1,ns));
                end
            end
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
    if oblate
        anm = anm_base .* -1 .* (u0.^2 + 1); % cnm
        L = legendre_otc(p,1j.*u0,1,1,1);
    else
        anm = anm_base .* (u0.^2 - 1); % bnm
        L = legendre_otc(p,u0,1,1,1);
    end

    P = L{1}; Q = L{2}; dP = L{3}; dQ = L{4};

    % Coefficient of D' of surface
    lambda_nm_prime = anm.' .* sqrt(u0.^2 - 1) .* dP.' .* dQ.' ./ a;
    lambda_nm = anm .* ((dQ .* P + dP .* Q) ./ 2) ./ a;
    lambda_n1m = lambda_nm;
end
     
function [lambda_nm_prime, lambda_nm, lambda_n1m] = DPspectrum_away(p,u0,a,u_x,v_x,nu,oblate)
    %{
        Calculates the normal derivative of the Laplace DLP with arbitrary normal vector.
        Resulting vectors are of size N x sp, where N is the number of target points.

        nu -> [nu_u,nu_v,nu_phi] Nx x 3 normal vector in spheroidal basis

        Note that terms relating to f_n^mY_n^m are ignored and are handled outside of this function.
    %}
    sp=(p+1)^2;
    ii = (1:sp)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;
    anm_base = factorial(nn-mm)./factorial(nn+mm) .* ((-1) .^ (mm));
    nu_u = nu(:,1); nu_v = nu(:,2); nu_phi = nu(:,3);

    if ~oblate % PROLATE CASE
        anm = anm_base .* (u0.^2-1); % bnm
        L = legendre_otc(p,u0,1,1,1);
        if abs(u_x)-u0 < -1e-14 % interior
            gnm = L{4}; % Q'(u_0)
        elseif abs(u_x)-u0 > 1e-14 % exterior
            gnm = L{3}; % P'(u_0)
        end

        if norm(abs(u_x)-u0)<1e-14 % on surface with arbitrary nu
            [lambda_nm_prime,lambda_nm,lambda_n1m] = DPspectrum(p,u0,a,oblate);
            lambda_nm_prime = lambda_nm_prime ./ sqrt(u_x.^2-v_x.^2) .* nu_u;
            
            lambda_nm_term1 = ((nn'+1) .* v_x) ./ sqrt((u_x.^2 - v_x.^2).*(1 - v_x.^2)) .* nu_v;
            lambda_nm_term2 = 1j.*mm'.*nu_phi./sqrt((u_x.^2-1).*(1-v_x.^2));
            lambda_nm = (lambda_nm.') .* (lambda_nm_term1 + lambda_nm_term2);

            lambda_n1m = -(lambda_n1m.') .* (nn'-mm'+1).*nu_v./sqrt((u_x.^2-v_x.^2).*(1-v_x.^2));
        else % off-surface.
            lambda_nm_prime = anm.' .* gnm.' .* sqrt((u_x.^2-1)./(u_x.^2-v_x.^2)) .* nu_u ./ a;

            lambda_nm_term1 = ((nn'+1) .* v_x) ./ sqrt((u_x.^2 - v_x.^2).*(1 - v_x.^2)) .* nu_v;
            lambda_nm_term2 = 1j .* mm' ./ sqrt((u_x.^2-1).*(1-v_x.^2)) .* nu_phi;
            lambda_nm = anm.' .* gnm.' .*(lambda_nm_term1 + lambda_nm_term2) ./ a;
            
            lambda_n1m = -anm.' .* gnm.' .* (nn'-mm'+1)./sqrt((u_x.^2-v_x.^2).*(1-v_x.^2)) .*nu_v ./ a;
        end
    else % OBLATE CASE
        anm = -1 .* anm_base .* (u0.^2+1); % cnm
        L = legendre_otc(p,1j.*u0,1,1,1);
        if abs(u_x)-u0 <- 1e-14 % interior
            gnm=L{4}; % Q'(iu_0)
        elseif abs(u_x)-u0>1e-14 % exterior
            gnm=L{3}; % P'(iu_0)
        end
        if norm(abs(u_x)-u0)<1e-14
            [lambda_nm_prime,lambda_nm,lambda_n1m]=DPspectrum(p,u0,a,oblate);
            lambda_nm_prime = lambda_nm_prime ./ sqrt(u_x.^2+v_x.^2) .* nu_u;
            lambda_nm = lambda_nm.'.*((nn'+1).*v_x.*nu_v./sqrt((u_x.^2+v_x.^2).*(1-v_x.^2))+1j.*mm'.*nu_phi./sqrt((u_x.^2+1).*(1-v_x.^2)));
            lambda_n1m = -lambda_n1m.'.*(nn'-mm'+1).*nu_v./sqrt((u_x.^2+v_x.^2).*(1-v_x.^2));
        else
            lambda_nm_prime = 1j.*anm.'.*gnm.'.*sqrt((u_x.^2+1)./(u_x.^2+v_x.^2)).*nu_u ./ a;
            lambda_nm = anm.'.*gnm.'.*((nn'+1).*v_x.*nu_v./sqrt((u_x.^2+v_x.^2).*(1-v_x.^2))+1j.*mm'.*nu_phi./sqrt((u_x.^2+1).*(1-v_x.^2))) ./ a;
            lambda_n1m = -anm.'.*gnm.'.*(nn'-mm'+1).*nu_v./sqrt((u_x.^2+v_x.^2).*(1-v_x.^2)) ./ a;
        end
    end
end

function [Fr, Fp]=solid_harmonic_prime(p, u0, u_x, oblate)
    %{
        Solid spheroidal harmonics can be written as f_n^m(u)Y_n^m(v, phi).
        This function returns what Fr = f_n^m is (and its derivative as Fp), depending 
        on whether we are in the exterior or interior.
    %}
    Fr=ones(size(u_x,1),(p+1)^2); Fp=ones(size(u_x,1),(p+1)^2);

    if oblate
        u_x = 1j.*u_x;
    end

    if abs(u_x)-u0 < -1e-14 % Interior
        PQ=legendre_otc(p,u_x,1,1);
        P=PQ{1}; dP=PQ{3};
        Fr=P.'; Fp=dP.';
    elseif abs(u_x)-u0 > 1e-14 % Exterior
        PQ=legendre_otc(p,u_x,1,1,1);
        Q=PQ{2}; dQ=PQ{4};
        Fr=Q.'; Fp=dQ.';
    end
end
    
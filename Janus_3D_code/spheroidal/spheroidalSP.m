function [SP,varargout]=spheroidalSP(params,X,nu,varargin)
%--------------------------------------------------------------------%
% spheroidalSP computes the normal derivative of the laplace single Layer 
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
%    a cell of length ns or an nt x 3 x ns matrix.
%    (o) nu_2,... = each input same format as nu, and prompts an additional
%    entry to the output. If nargin==3, then nargout==1.
%    
% Returns SP as a matrix if X is a matrix or not provided, or as a cell if
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

% Do basis transformations for all possible u0 values we might have.
% -----------------------------------------------------------------------
all_u0s=unique(u0);
Gshc = zeros(size(shc));
Gmats= zeros(sp,sp,ns);

% Change to basis that diagonalizes the single layer operator
for i=1:length(all_u0s)
    this_u0=all_u0s(i);

    %Calculate basis transformation matrix G: Ynm -> Ynm / sqrt(u0^2-v^2)
    %or /sqrt(u0^2+v^2) for oblate
    this_G_pro = Gmatrix(p,this_u0,0,0);
    this_G_obl = Gmatrix(p,this_u0,0,1);

    % keep track of indices and shc of particles with the same u0
    this_index = find(u0==this_u0); 
    obl_index = nonzeros(this_index.*oblate(this_index));
    pro_index = nonzeros(this_index.*~oblate(this_index));
    obl_shc = reshape(shc(:,:,obl_index),sp,[],1);
    pro_shc = reshape(shc(:,:,pro_index),sp,[],1);

    %Apply basis transformation to spheroidal harmonic coefficients
    obl_Gshc = this_G_obl\obl_shc;
    pro_Gshc = this_G_pro\pro_shc;
    
    Gshc(:,:,obl_index) = reshape(obl_Gshc,sp,[],length(obl_index));
    Gshc(:,:,pro_index) = reshape(pro_Gshc,sp,[],length(pro_index));
    Gmats(:,:,obl_index) = repmat(this_G_obl,1,1,length(obl_index));
    Gmats(:,:,pro_index) = repmat(this_G_pro,1,1,length(pro_index));
end

% Gshc = shc ./ sqrt(u_0(k).^2-v.^2);

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
        Yn=Ynm(n,[],acos(v_x)',phi_x);
        Y(:,n^2+1:(n+1)^2)=Yn;
    end

    % Calculate spectra on surface with outward normal
    [spectra_surf,~,~]=SPspectrum(p,u0,oblate);

    %reshape so that each column of SP_coefs corresponds to each function 
    % on each spheroid
    spectra_surf = reshape(spectra_surf,sp,1,ns);
    spectra_matrix = repmat(spectra_surf,1,nf,1);
    SP_coefs = spectra_matrix.*Gshc; 
    %%%%%%%%%
    % SP_coefs = pagemtimes(Gmats,SP_coefs);
    %%%%%%%%%%
    SP_coefs = reshape(SP_coefs,sp,[],1);

    %multiply coefficients with spheroidal harmonics
    SP=Y*SP_coefs;

    % reshape to match original shape of sigma
    SP = reshape(SP,[],nf,ns);

    %%%%%%%%%%%%
    % assuming nf=1
    coef_mat=zeros(np,1,ns);
    for i=1:ns
        coef_mat(:,:,i)=oblate(i).*1./sqrt(u0(i)^2+v_x.^2)+~oblate(i).*1./sqrt(u0(i)^2-v_x.^2);
    end
    % pro_coef=1./sqrt(u0^2-v_x.^2);
    % obl_coef=1./sqrt(u0^2+v_x.^2);
    % if size(oblate,1)~=1
    %     oblate=reshape(oblate,1,[]);
    % end
    % coef_mat=obl_coef*oblate+pro_coef*~oblate;
    % coef_mat=reshape(coef_mat,size(coef_mat,1),1,size(coef_mat,2));
    SP = coef_mat.*SP;

    %%%%%%%%%%%%%%%%

    
    if isReal
        SP=real(SP);
    end

%  Self Eval, arbitrary normal.
% ------------------------------------------------------------------------
elseif isempty(X)
    % When want to get surface dS/dnu with arbitrary normal,
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

    SP=cell(1,ns);

    if more_nu
        varargout=cell(1,nvarin);
        for nu_ind=1:nvarin
            varargout{nu_ind}=cell(1,ns);
        end
    end

    % TODO: why does the phase shift of phi not matter?
    %%%%%%%%%%%
    [theta2,phi_k]=gl_grid(p);
    v_k=cos(theta2);

    nt_r=length(v_k);
    Yr=zeros(nt_r*2,sp);

    for n=0:p  %loop over terms in spheroidal harmonic expansion
        % Yn=Ynm(n,[],acos(v_x_r)',phi_x_r);
        Yn=Ynm(n,[],real(acos(v_k))',phi_k); % v_x_r exceeds [-1,1] by 1e-8, but acos() returns imaginary values. Impose real values for Ynm.
        Yr(1:nt_r,n^2+1:(n+1)^2)=Yn;
        % Yn1 = Ynm(n+1,-n:n,acos(v_x_r)',phi_x_r);
        Yn1=Ynm(n+1,-n:n,real(acos(v_k))',phi_k); 
        yn1_scale=sqrt((2*n+1)/(2*n+3).*(n+(-n:n)+1)./(n-(-n:n)+1));
        Yr(nt_r+1:end,n^2+1:(n+1)^2) = yn1_scale.*Yn1;
    end
    %%%%%%%%%%%%%%%%%%%%%

    for k=1:ns  %loop over each spheroid surface we want to evaluate
        [~,Xself_k]=params.get_X(k);
        S=cart2spheroidal(Xself_k,a(k),oblate(k));
        % v_k=S(:,2);
        % if abs(v_k-real(v_k))>1e-10
        %     error("v_k imaginary\n")
        % end
        % v_k=real(v_k);
        % phi_k=S(:,3);

        u_k=u0(k).*ones(size(v_k));

        nu_list=cell(1,nvarin+1); nu_sph_list=nu_list;
        nu_list{1}=nu_t{k};
        nu_sph_list{1}=cartNu2spheroidal(nu_t{k},S,a(k),oblate(k));
        for nu_ind=1:nvarin
            nu_list{nu_ind+1}=nu_extra_cells{nu_ind}{k};
            [nu_sph_temp,~]=cartNu2spheroidal(nu_list{nu_ind+1},S,a(k),oblate(k));
            nu_sph_list{nu_ind+1}=nu_sph_temp;
        end
      
        %%%%%%%%%%%%% where Yr was %%%%%%%%%%%

        for nu_ind=1:nvarin+1
            [spectra_nm_prime,spectra_nm,spectra_n1m] = SPspectrum_away(p,u0(k),u_k,v_k,nu_sph_list{nu_ind},oblate(k));
            
            FYr = (spectra_nm_prime + spectra_nm).*Yr(1:nt_r,:)+spectra_n1m.*Yr(nt_r+1:end,:);
            SP_k=FYr*Gshc(:,:,k);

            if isReal
                SP_k=real(SP_k);
            end
            
            if nu_ind==1
                SP{k}=SP_k;
            else
                varargout{nu_ind-1}{k}=SP_k;
            end
        end

    end

    if ~isa(nu, "cell")
        SP = cell2mat(reshape(SP,1,1,ns));

        if more_nu
            for nu_ind=1:nvarin
                varargout{nu_ind}=cell2mat(reshape(varargout{nu_ind},1,1,ns));
            end
        end
    end


% Target points given
% -----------------------------------------------------------------------
elseif nargin > 1

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

    SP=cell(1,ns);

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
        SPk = zeros(ntk,nf);
        if more_nu
            SP_extra_at_k=zeros(ntk,nf,nvarin);
            SP_extra_at_k=mat2cell(SP_extra_at_k,ntk,nf,ones(1,nvarin));
        end

        if ntk > 0
            if d~=3
                error("Dimensions of target point array should be N x 3.")
            end
            
            % Convert targets to spheroidal coords
            S=cart2spheroidal(Xtk,a(k),oblate(k));
            u_x=S(:,1);

            indices_interior = (u_x < u0(k) - 9e-12);
            indices_surface = (abs(u_x-u0(k)) <= 9e-12);
            indices_exterior = (u_x > u0(k) + 9e-12);
            
            % split up interior/surface/exterior
            S_int=S(indices_interior,:);
            S_surf=S(indices_surface,:);
            S_ext=S(indices_exterior,:);
            Nu_int=nu_k(indices_interior,:);
            Nu_surf=nu_k(indices_surface,:);
            Nu_ext=nu_k(indices_exterior,:);

            if more_nu
                Nu_int_extra=cell(1,nvarin);
                Nu_surf_extra=cell(1,nvarin);
                Nu_ext_extra=cell(1,nvarin);
                for nu_ind=1:nvarin
                    Nu_int_extra{nu_ind}=nu_extra_cells{nu_ind}{k}(indices_interior,:);
                    Nu_surf_extra{nu_ind}=nu_extra_cells{nu_ind}{k}(indices_surface,:);
                    Nu_ext_extra{nu_ind}=nu_extra_cells{nu_ind}{k}(indices_exterior,:);
                end
            end
            
            % If there's no points in a particular region, don't do any
            % computations
            do_int=~isempty(S_int);
            do_surf=~isempty(S_surf); 
            do_ext=~isempty(S_ext);
        
            % Split coordinates into 3 regions: We have a SP for each region
            regions=[do_int,do_surf,do_ext];
            Sregions={S_int,S_surf,S_ext};
            Nuregions={Nu_int,Nu_surf,Nu_ext};
            SPregions=cell(3,1);
            if more_nu
                Nuregions_extra={Nu_int_extra,Nu_surf_extra,Nu_ext_extra};
                SPregions_extra={cell(1,nvarin),cell(1,nvarin),cell(1,nvarin)}; %SPregions_extra{1}={SP_int from nu_1},{SP_int from nu_2},...
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

                    % if more_nu
                    %     nu_r_cell=Nuregions_extra{r};
                    %     nu_r_extra=cell(1,nvarin);
                    %     for nu_ind=1:nvarin
                    %         nu_r_cart_extra=nu_r_cell{nu_ind};
                    %         [nu_r_sph_extra,~] = cartNu2spheroidal(nu_r_cart_extra,Sr,a(k),oblate(k));
                    %         nu_r_extra{nu_ind}=nu_r_sph_extra;
                    %     end
                    % end

                    % %% testing
                    % Xcurr=spheroidal2cart(Sr,a(k),oblate);
                    % quiver3(Xcurr(:,1),Xcurr(:,2),Xcurr(:,3),nu_r_cart(:,1),nu_r_cart(:,2),nu_r_cart(:,3))
                    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    [spectra_nm_prime,spectra_nm,spectra_n1m] = SPspectrum_away(p,u0(k),u_x_r,v_x_r,nu_r_sph,oblate(k));
            
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

                    % solid spheroidal harmonics
                    % FYr=Fr.*Yr;
                    % FYr_matrix = repmat(FYr,1,nf);
                    FYr = (spectra_nm_prime.*Fp + spectra_nm.*Fr).*Yr(1:nt_r,:)...
                        +spectra_n1m.*Fr.*Yr(nt_r+1:end,:);

                    % Spheroidal harmonic coefficients of S'[fac*sigma]
                    SPregions{r}=FYr*Gshc(:,:,k);

                    % if there are more nu vectors to calculate SP with,
                    % use already computed Fr,Yr at target points, and
                    % evaluate spectrum using new normal
                    if more_nu
                        nu_r_cell=Nuregions_extra{r};
                        SP_r_cell=cell(1,nvarin);
                        for nu_ind=1:nvarin
                            nu_r_cart=nu_r_cell{nu_ind};
                            [nu_r_sph,~] = cartNu2spheroidal(nu_r_cart,Sr,a(k),oblate(k));
                            [spectra_nm_prime,spectra_nm,spectra_n1m] = SPspectrum_away(p,u0(k),u_x_r,v_x_r,nu_r_sph,oblate(k));
                            FYr = (spectra_nm_prime.*Fp + spectra_nm.*Fr).*Yr(1:nt_r,:)+spectra_n1m.*Fr.*Yr(nt_r+1:end,:);
                            SP_r_cell{nu_ind}=FYr*Gshc(:,:,k);
                        end
                        SPregions_extra{r}=SP_r_cell;
                    end

                end
            end
        
            % Recombine all SPs from interior/exterior/surface
            SPk(indices_interior,:) = SPregions{1};
            SPk(indices_surface,:) = SPregions{2};
            SPk(indices_exterior,:) = SPregions{3};

            if more_nu
                for nu_ind=1:nvarin
                    SPk_temp=SP_extra_at_k{nu_ind};
                    SPk_temp(indices_interior,:) = SPregions_extra{1}{nu_ind};
                    SPk_temp(indices_surface,:) = SPregions_extra{2}{nu_ind};
                    SPk_temp(indices_exterior,:) = SPregions_extra{3}{nu_ind};
                    SP_extra_at_k{nu_ind}=SPk_temp;
                    if isReal
                        SP_extra_at_k{nu_ind}=real(SPk_temp);
                    end
                end
            end
    
            if isReal
                SPk=real(SPk);
            end
        end
        
        %evluation for particle k
        SP{k} = SPk;

        if more_nu
            for nu_ind=1:nvarin
                varargout{nu_ind}{k}=SP_extra_at_k{nu_ind};
            end
        end
                

    end

    if ~isa(X, "cell")
        SP = cell2mat(reshape(SP,1,1,ns));

        if more_nu
            for nu_ind=1:nvarin
                varargout{nu_ind}=cell2mat(reshape(varargout{nu_ind},1,1,ns));
            end
        end
    end

end
%--------------------------------------------------------------------%


end


function [lambda_nm_prime,lambda_nm,lambda_n1m]=SPspectrum(p,u0,oblate)
% Note May 20 2024: Changed output format because off surface spectra not
% diagonal. Outputs are coefficients for fnm'Ynm, fnmYnm, fnmY(n+1)m as the 
% average of exterior limit and interior limit.
    sp=(p+1)^2;
    ii = (1:sp)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;
    anm=factorial(nn-mm)./factorial(nn+mm).*(-1).^(mm).*(oblate.*1j.*sqrt(u0.^2+1)+~oblate.*sqrt(u0.^2-1));
    % anm=factorial(nn-mm)./factorial(nn+mm).*(-1).^(mm).*(oblate.*(-1).*(u0.^2+1)+~oblate.*(u0.^2-1));
    L_obl=legendre_otc(p,1j.*u0,1,1,1);
    u0_pro=u0+oblate; % if i-th u0 for oblate, will add 1 to make sure it is >1. actual value won't matter because only taking indices corresponding to prolates from this PQ.
    L_pro=legendre_otc(p,u0_pro,1,1,1);
    % L_pro=legendre_otc(p,u0,1,1,1);
    Pp=L_pro{1}; Qp=L_pro{2}; dPp=L_pro{3}; dQp=L_pro{4};
    Po=L_obl{1}; Qo=L_obl{2}; dPo=L_obl{3}; dQo=L_obl{4};
    % lambda_int=anm.*(Qo*diag(oblate)+Qp*diag(~oblate));
    lambda_nm_prime=anm.*(oblate.*1j.*sqrt(u0.^2+1)+~oblate.*sqrt(u0.^2-1))./2.*((Po.*dQo+dPo.*Qo)*diag(oblate)+(Pp.*dQp+dPp.*Qp)*diag(~oblate));
    % lambda_ext=anm.*(Po*diag(oblate)+Pp*diag(~oblate)); 
    lambda_nm=anm.*(Po.*Qo*diag(oblate)+Pp.*Qp*diag(~oblate));
    lambda_n1m=anm.*(Po.*Qo*diag(oblate)+Pp.*Qp*diag(~oblate));
end
 
function [lambda_nm_prime, lambda_nm,lambda_n1m] = SPspectrum_away(p,u0,u_x,v_x,nu,oblate)
    % nu = [nu_u,nu_v,nu_phi] Nx x 3 normal vector in spheroidal basis
    sp=(p+1)^2;
    ii = (1:sp)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;
    anm_base=factorial(nn-mm)./factorial(nn+mm).*(-1).^mm;
    nu_u = nu(:,1); nu_v = nu(:,2); nu_phi = nu(:,3);
    if ~oblate
        
        if norm(abs(u_x)-u0)<9e-12 % on surface with arbitrary nu
            [lambda_nm_prime,lambda_nm,lambda_n1m]=SPspectrum(p,u0,oblate);
            lambda_nm_prime = lambda_nm_prime.'./sqrt(u_x.^2-v_x.^2).*nu_u;
            lambda_nm = lambda_nm.'.*((nn'+1).*v_x.*nu_v./sqrt((u_x.^2-v_x.^2).*(1-v_x.^2))+1j.*mm'.*nu_phi./sqrt((u_x.^2-1).*(1-v_x.^2)));
            lambda_n1m = -lambda_n1m.'.*(nn'-mm'+1).*nu_v./sqrt((u_x.^2-v_x.^2).*(1-v_x.^2));
        else
            anm=anm_base.*sqrt(u0.^2-1);
            PQ0 = legendre_otc(p,u0,1,1,1);
            if abs(u_x)-u0<1e-14 % interior
                gnm=PQ0{2};
            elseif abs(u_x)-u0>1e-14 % exterior
                gnm=PQ0{1};
            end
            lambda_nm_prime = anm.'.*gnm.'.*sqrt((u_x.^2-1)./(u_x.^2-v_x.^2)).*nu_u;
            lambda_nm = anm.'.*gnm.'.*((nn'+1).*v_x.*nu_v./sqrt((u_x.^2-v_x.^2).*(1-v_x.^2))+1j.*mm'.*nu_phi./sqrt((u_x.^2-1).*(1-v_x.^2)));
            lambda_n1m = -anm.'.*gnm.'.*(nn'-mm'+1).*nu_v./sqrt((u_x.^2-v_x.^2).*(1-v_x.^2));
        end

    else
        
        if norm(abs(u_x)-u0)<9e-12
            [lambda_nm_prime,lambda_nm,lambda_n1m]=SPspectrum(p,u0,oblate);
            lambda_nm_prime = lambda_nm_prime.'./sqrt(u_x.^2+v_x.^2).*nu_u;
            lambda_nm = lambda_nm.'.*((nn'+1).*v_x.*nu_v./sqrt((u_x.^2+v_x.^2).*(1-v_x.^2))+1j.*mm'.*nu_phi./sqrt((u_x.^2+1).*(1-v_x.^2)));
            lambda_n1m = -lambda_n1m.'.*(nn'-mm'+1).*nu_v./sqrt((u_x.^2+v_x.^2).*(1-v_x.^2));
        else
            anm=1j.*anm_base.*sqrt(u0.^2+1);
            PQ0 = legendre_otc(p,1j.*u0,1);
            if abs(u_x)-u0<1e-14 % interior
                gnm=PQ0{2};
            elseif abs(u_x)-u0>1e-14 % exterior
                gnm=PQ0{1};
            end
            lambda_nm_prime = 1j.*anm.'.*gnm.'.*sqrt((u_x.^2+1)./(u_x.^2+v_x.^2)).*nu_u;
            lambda_nm = anm.'.*gnm.'.*((nn'+1).*v_x.*nu_v./sqrt((u_x.^2+v_x.^2).*(1-v_x.^2))+1j.*mm'.*nu_phi./sqrt((u_x.^2+1).*(1-v_x.^2)));
            lambda_n1m = -anm.'.*gnm.'.*(nn'-mm'+1).*nu_v./sqrt((u_x.^2+v_x.^2).*(1-v_x.^2));
        end
    end

    
end

% Storage before oblate changes.
%{
function [lambda_nm_prime, lambda_nm,lambda_n1m] = SPspectrum_away(p,u0,u_x,v_x,nu,oblate)
    % nu = [nu_u,nu_v,nu_phi] Nx x 3 normal vector in spheroidal basis
    sp=(p+1)^2;
    ii = (1:sp)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;
    anm_base=factorial(nn-mm)./factorial(nn+mm).*(-1).^mm;
    nu_u = nu(:,1); nu_v = nu(:,2); nu_phi = nu(:,3);
    if ~oblate
        
        if norm(abs(u_x)-u0)<1e-14 % on surface with arbitrary nu
            [lambda_nm_prime,lambda_nm,lambda_n1m]=SPspectrum(p,u0,oblate);
            lambda_nm_prime = lambda_nm_prime.'./sqrt(u_x.^2-v_x.^2).*nu_u;
            lambda_nm = lambda_nm.'.*((nn'+1).*v_x.*nu_v./sqrt((u_x.^2-v_x.^2).*(1-v_x.^2))+1j.*mm'.*nu_phi./sqrt((u_x.^2-1).*(1-v_x.^2)));
            lambda_n1m = -lambda_n1m.'.*(nn'-mm'+1).*nu_v./sqrt((u_x.^2-v_x.^2).*(1-v_x.^2));
        else
            anm=anm_base.*sqrt(u0.^2-1);
            PQ0 = legendre_otc(p,u0,1,1,1);
            if abs(u_x)-u0<1e-14 % interior
                gnm=PQ0{2};
            elseif abs(u_x)-u0>1e-14 % exterior
                gnm=PQ0{1};
            end
            lambda_nm_prime = anm.'.*gnm.'.*sqrt((u_x.^2-1)./(u_x.^2-v_x.^2)).*nu_u;
            lambda_nm = anm.'.*gnm.'.*((nn'+1).*v_x.*nu_v./sqrt((u_x.^2-v_x.^2).*(1-v_x.^2))+1j.*mm'.*nu_phi./sqrt((u_x.^2-1).*(1-v_x.^2)));
            lambda_n1m = -anm.'.*gnm.'.*(nn'-mm'+1).*nu_v./sqrt((u_x.^2-v_x.^2).*(1-v_x.^2));
        end

    else
        
        if norm(abs(u_x)-u0)<1e-14
            [lambda_nm_prime,lambda_nm,lambda_n1m]=SPspectrum(p,u0,oblate);
            lambda_nm_prime = lambda_nm_prime./sqrt(u_x.^2+v_x.^2).*nu_u;
            lambda_nm = lambda_nm.*((nn'+1).*v_x.*nu_v./sqrt((u_x.^2+v_x.^2).*(1-v_x.^2))+1j.*mm'.*nu_phi./sqrt((u_x.^2+1).*(1-v_x.^2)));
            lambda_n1m = -lambda_n1m.*(nn'-mm'+1).*nu_v./sqrt((u_x.^2+v_x.^2).*(1-v_x.^2));
        else
            anm=1j.*anm_base.*sqrt(u0.^2+1);
            PQ0 = legendre_otc(p,1j.*u0,1);
            if abs(u_x)-u0<1e-14 % interior
                gnm=PQ0{2};
            elseif abs(u_x)-u0>1e-14 % exterior
                gnm=PQ0{1};
            end
            lambda_nm_prime = 1j.*anm.'.*gnm.'.*sqrt((u_x.^2+1)./(u_x.^2+v_x.^2)).*nu_u;
            lambda_nm = anm.'.*gnm.'.*((nn'+1).*v_x.*nu_v./sqrt((u_x.^2+v_x.^2).*(1-v_x.^2))+1j.*mm'.*nu_phi./sqrt((u_x.^2+1).*(1-v_x.^2)));
            lambda_n1m = -anm.'.*gnm.'.*(nn'-mm'+1).*nu_v./sqrt((u_x.^2+v_x.^2).*(1-v_x.^2));
        end
    end

    
end
%}

function [F,Fp]=solid_harmonic_prime(p,u0,u_x,oblate)
    F=ones(size(u_x,1),(p+1)^2);
    Fp=ones(size(u_x,1),(p+1)^2);
    if nargin < 4
        oblate = false;
    end
    if oblate
        u_x = 1j.*u_x;
    end
    % if abs(u_x)-u0 < 1e-14 
    if abs(u_x)-u0 < -1e-14 
        PQ=legendre_otc(p,u_x,1,1);
        P=PQ{1}; dP=PQ{3};
        F=P.'; Fp=dP.';
    elseif abs(u_x)-u0 > 1e-14
        PQ=legendre_otc(p,u_x,1,1,1);
        Q=PQ{2}; dQ=PQ{4};
        F=Q.'; Fp=dQ.';
    end
end

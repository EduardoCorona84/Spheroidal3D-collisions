function DL=spheroidalDL(params,X)
%--------------------------------------------------------------------%
% spheroidalDL computes the laplace Double Layer potential of a density 
% 'sigma' on the surface of a spheroid.
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
%    
% Returns DL as a matrix if X is a matrix or not provided, or as a cell if
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
oblate=params.oblate;

[sp,nf,ns]=size(params.sigma_coefficients);
[~,nf2,ns2]=size(sigma);

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
    
    %if X not given, calculate on-surface points
    [theta_x,phi_x]=gl_grid(p);
    v_x=cos(theta_x);
    nt=length(theta_x);
    
    Y=zeros(nt,sp);
    for n=0:p  %loop over terms in spheroidal harmonic expansion
        Yn=Ynm(n,[],acos(v_x)',phi_x);
        Y(:,n^2+1:(n+1)^2)=Yn;
    end

    % Calculate spectra for interior, exterior, and surface
    [~,spectra_surf,~]=DLspectrum(p,u0,oblate);

    %reshape so that each column of DL_coefs corresponds to each function 
    % on each spheroid
    spectra_surf = reshape(spectra_surf,sp,1,ns);
    spectra_matrix = repmat(spectra_surf,1,nf,1);
    DL_coefs = spectra_matrix.*shc;

    DL_coefs = reshape(DL_coefs,sp,[],1);

    %multiply coefficients with spheroidal harmonics
    DL=Y*DL_coefs;

    % reshape to match original shape of sigma
    DL = reshape(DL,[],nf,ns);

    if isReal
        DL=real(DL);
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
    
    if length(Xt) ~= ns
        error("dimensions of target input and number of surfaces do not match. X should be a 1 x N cell or nt x 3 x N matrix.")
    end

    DL=cell(1,ns);

    % % Starting the actual computation
    %--------------------------------------------------------------------%

    % Calculate spectra for interior, exterior, and surface
    [spectra_int,spectra_surf,spectra_ext] = DLspectrum(p,u0,oblate);
    
    for k=1:ns  %loop over each spheroid surface we want to evaluate
        
        Xtk = Xt{k};
        [ntk,d]=size(Xtk);
        DLk = zeros(ntk,nf);

        if ntk > 0
            if d~=3
                error("Dimensions of target point array should be N x 3.")
            end
        
            spectra_regions={spectra_int(:,k),spectra_surf(:,k),spectra_ext(:,k)};

            % Convert targets to spheroidal coords
            S=cart2spheroidal(Xtk,a(k),oblate(k));
    
            u_x=S(:,1);
            
            % split up interior/surface/exterior
            % S_int=S(u_x < u0(k),:);
            % S_surf=S(u_x == u0(k),:);
            % S_ext=S(u_x > u0(k),:);
            S_int=S(u_x < u0(k)-1e-14,:);
            S_surf=S(abs(u_x-u0(k))<=1e-14,:);
            S_ext=S(u_x > u0(k)+1e-14,:);
            
            % If there's no points in a particular region, don't do any
            % computations
            do_int=~isempty(S_int);
            do_surf=~isempty(S_surf); 
            do_ext=~isempty(S_ext);
        
            % Split coordinates into 3 regions: We have a DL for each region
            regions=[do_int,do_surf,do_ext];
            Sregions={S_int,S_surf,S_ext};
            DLregions=cell(3,1);
            
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
            
                    Fr=solid_harmonic(p,u0(k),u_x_r,oblate(k));
                    nt_r=length(u_x_r);
                    Yr=zeros(nt_r,sp);
            
                    for n=0:p  %loop over terms in spheroidal harmonic expansion
                        % Yn=Ynm(n,[],acos(v_x_r)',phi_x_r);
                        Yn=Ynm(n,[],real(acos(v_x_r))',phi_x_r); % v_x_r exceeds [-1,1] by 1e-8, but acos() returns imaginary values. Impose real values for Ynm.
                        Yr(:,n^2+1:(n+1)^2)=Yn;
                    end

                    % solid spheroidal harmonics
                    FYr=Fr.*Yr;

                    % Spheroidal harmonic coefficients of D[sigma]
                    spectra_matrix = repmat(spectra_regions{r},1,nf);
                    DLcoefs_r = spectra_matrix.*shc(:,:,k);
                    
                    % multiply spheroidal harmonic coefficients and
                    % harmonic functions together
                    DLregions{r}=FYr*DLcoefs_r;
                end
            end
        
            % Recombine all DLs from interior/exterior/surface
            % DLk(u_x < u0(k),:) = DLregions{1};
            % DLk(u_x == u0(k),:) = DLregions{2};
            % DLk(u_x > u0(k),:) = DLregions{3};
            DLk(u_x < u0(k)-1e-14,:) = DLregions{1};
            DLk(abs(u_x-u0(k))<=1e-14,:) = DLregions{2};
            DLk(u_x > u0(k)+1e-14,:) = DLregions{3};
    
            if isReal
                DLk=real(DLk);
            end
        end
        
        %evluation for particle k
        DL{k} = DLk;

    end

    if ~isa(X, "cell")
        DL = cell2mat(reshape(DL,1,1,ns));
    end

end
%--------------------------------------------------------------------%


% % Starting the actual computation
% %--------------------------------------------------------------------%
% 
% for k=1:ns  %loop over each spheroid surface we want to evaluate
% 
%     Xtk = Xt(sep(:,k)==0,:,k);
% 
%     % Calculate spectra for interior, exterior, and surface
%     [spectra_int,spectra_surf,spectra_ext]=DLspectrum(p,u0(k));
%     spectra_regions={spectra_int,spectra_surf,spectra_ext};
%     
%     % sort targets by increasing u so we can split up interior/surface/exterior
%     % Convert targets to spheroidal coords
%     if strcmp(target_coords, 'cart')
%         S=cart2spheroidal(Xtk,1/u0(k));
%     elseif strcmp(target_coords, 'spheroidal')
%         S=Xtk;
%     else 
%         error("Invalid coordinates given. Use 'cart' for cartesian or 'spheroidal' for spheroidal.")
%     end
%     
%     u_x=S(:,1);
%     [u_x_sort,j_inc_u]=sort(u_x);
% 
%     %we will use this later to resort targets back into original order
%     j_unsort(j_inc_u)=1:length(u_x); 
%     
%     % Sorted targets in spheroidal coordinates
%     S_sort=S(j_inc_u,:);
%     S_int=S_sort(u_x<u0(k),:);
%     S_surf=S_sort(u_x==u0(k),:);
%     S_ext=S_sort(u_x>u0(k),:);
% 
%     % If there's no points in a particular region, don't do any
%     % computations
%     do_int=~isempty(S_int); 
%     do_surf=~isempty(S_surf); 
%     do_ext=~isempty(S_ext);
% 
%     % Split coordinates into 3 regions: We have a DL for each region
%     regions=[do_int,do_surf,do_ext];
%     Sregions={S_int,S_surf,S_ext};
%     DLregions=cell(3,1);
%     
%     % Loop over each region
%     for r=1:3
%         if regions(r)
%             Sr=Sregions{r};
%             u_x_r=Sr(:,1);
%             v_x_r=Sr(:,2);
%             phi_x_r=Sr(:,3);
%     
%             Fr=solid_harmonic(p,u0(k),u_x_r);
%             nt_r=length(u_x_r);
%             Yr=zeros(nt_r,sp);
%     
%             for n=0:p  %loop over terms in spheroidal harmonic expansion
%                 Yn=Ynm(n,[],acos(v_x_r)',phi_x_r);
%                 Yr(:,n^2+1:(n+1)^2)=Yn;
%             end
%     
%             FYr=Fr.*Yr;
%             DLcoefs_r=spectra_regions{r}.*shc(:,k);
%             DLregions{r}=FYr*DLcoefs_r;
%         end
%     end
%     
%     % For debugging:
% %----------------------------------
% %         sprintf('size Y: ')
% %         disp(size(Y))
% %         sprintf('size spectra: ')
% %         disp(size(spectra))
% %         sprintf('size shc: ')
% %         disp(size(shc))
% %         sprintf('size FY: ')
% %         disp(size(FY))
% %---------------------------------
% 
%     % Recombine all DLs from interior/exterior/surface and resort to
%     % original order of targets
%     DLk_unsorted=cat(1,DLregions{1},DLregions{2},DLregions{3});
%     DL(sep(:,k)==0,k)=DLk_unsorted(j_unsort,:);
% 
% %     fprintf('completed evaluation of surface density %d\n',k)
% end


end



function [lambda_int,lambda_surf,lambda_ext]=DLspectrum(p,u0,oblate)
    sp=(p+1)^2;
    ii = (1:sp)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;
    anm=factorial(nn-mm)./factorial(nn+mm).*(-1).^(mm).*(oblate.*(-1).*(u0.^2+1)+~oblate.*(u0.^2-1));
    L_obl=legendre_otc(p,1j.*u0,1,1,1);
    u0_pro=u0+oblate; % if i-th u0 for oblate, will add 1 to make sure it is >1. actual value won't matter.
    L_pro=legendre_otc(p,u0_pro,1,1,1);
    % L_pro=legendre_otc(p,u0,1,1,1);
    Pp=L_pro{1}; Qp=L_pro{2}; dPp=L_pro{3}; dQp=L_pro{4};
    Po=L_obl{1}; Qo=L_obl{2}; dPo=L_obl{3}; dQo=L_obl{4};

    lambda_int=anm.*(dQo*diag(oblate)+dQp*diag(~oblate));
    lambda_surf=anm./2.*((Po.*dQo+dPo.*Qo)*diag(oblate)+(Pp.*dQp+dPp.*Qp)*diag(~oblate));
    lambda_ext=anm.*(dPo*diag(oblate)+dPp*diag(~oblate));   
end

function [lambda_int,lambda_surf,lambda_ext]=DLspectrum_old(p,u0,oblate)
    sp=(p+1)^2;
    ii = (1:sp)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;
    if oblate
    	anm=factorial(nn-mm)./factorial(nn+mm).*(-1).^(mm+1).*(u0.^2+1);
        L=legendre_otc(p,1j.*u0,1,1,1);
    else
	    anm=factorial(nn-mm)./factorial(nn+mm).*(-1).^mm .*(u0.^2-1);
        L=legendre_otc(p,u0,1,1,1);
    end
    P=L{1}; Q=L{2}; dP=L{3}; dQ=L{4};

    lambda_int=anm.*dQ;
    lambda_surf=anm./2.*(P.*dQ+dP.*Q);
    lambda_ext=anm.*dP;   
end

function F=solid_harmonic(p,u0,u_x,oblate)
    F=1;
    if oblate
	    u_x = 1j.*u_x;
    end
    if abs(u_x) - u0 < -1e-14
        PQ=legendre_otc(p,u_x,0);
        P=PQ{1};
        F=P.';
    elseif abs(u_x) - u0 > 1e-14
        PQ=legendre_otc(p,u_x,1);
        Q=PQ{2};
        F=Q.';
    end
end

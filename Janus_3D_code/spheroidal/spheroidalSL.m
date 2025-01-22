function SL=spheroidalSL(params,X)
%--------------------------------------------------------------------%
% spheroidlSL computes the laplace Single Layer potential of a density 
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
%    (o) X = (x,y,z) coordinates of target evaluations. Is size [nt,3]
%    
% Returns SL of size [nt,ns] for self-evaluation or [nt,1] for provided X.
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


% Do basis transformations for all possible u0 values we might have.
% -----------------------------------------------------------------------
all_u0s=unique(u0);
Gshc = zeros(size(shc));
% Change to basis that diagonalizes the single layer operator

% for i=1:length(all_u0s)
%     this_u0=all_u0s(i);
%     this_if_obl=oblate(i);
% 
%     %Calculate basis transformation matrix G: Ynm -> Ynm / sqrt(u0^2-v^2
%     this_G = Gmatrix(p,this_u0,0,this_if_obl);
% 
%     % keep track of indices and shc of particles with the same u0
%     this_index = find(u0==this_u0); 
%     this_shc = reshape(shc(:,:,this_index),sp,[],1);
% 
%     %Apply basis transformation to spheroidal harmonic coefficients
%     this_Gshc = this_G\this_shc;
% 
%     Gshc(:,:,this_index) = reshape(this_Gshc,sp,[],length(this_index));
% end

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
end


% No target points given; do self evaluation
% -----------------------------------------------------------------------
if nargin==1

    %if X not given, calculate on-surface
    [theta_x,phi_x]=gl_grid(p);
    v_x=cos(theta_x);
    nt=length(theta_x);

    Y=zeros(nt,sp);
    for n=0:p  %loop over terms in spheroidal harmonic expansion
        Yn=Ynm(n,[],acos(v_x)',phi_x);
        Y(:,n^2+1:(n+1)^2)=Yn;
    end

    % Calculate spectra for interior, exterior, and surface
    [~,spectra_surf,~]=SLspectrum(p,u0,a,oblate);

    %reshape so that each column of DL_coefs corresponds to each function 
    % on each spheroid
    spectra_surf = reshape(spectra_surf,sp,1,ns);
    spectra_matrix = repmat(spectra_surf,1,nf,1);
    SL_coefs = spectra_matrix.*Gshc;
    SL_coefs = reshape(SL_coefs,sp,[],1);
    % display(SL_coefs);

    %multiply coefficients with spheroidal harmonics
    SL=Y*SL_coefs;

    % reshape to match original shape of sigma
    SL = reshape(SL,[],nf,ns);

    if isReal
        SL=real(SL);
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
    
    SL=cell(1,ns);

    % % Starting the actual computation
    %--------------------------------------------------------------------%

    % Calculate spectra for interior, exterior, and surface
    [spectra_int,spectra_surf,spectra_ext] = SLspectrum(p,u0,a,oblate);

    for k=1:ns  %loop over each spheroid surface we want to evaluate
        
        Xtk = Xt{k};
        [ntk,d]=size(Xtk);
        SLk = zeros(ntk,nf);

        if ntk > 0
            if d~=3
                error("Dimensions of target point array should be N x 3.")
            end
    
            spectra_regions={spectra_int(:,k),spectra_surf(:,k),spectra_ext(:,k)};

            % Convert targets to spheroidal coords
            S=cart2spheroidal(Xtk,a(k),oblate(k));
            
            u_x=S(:,1);
            % %%%%%%%%%%%%%%%%%%%
            % u_x=u0(k).*ones(size(S,1),1);
            % %%%%%%%%%%%%%%%%%%%%%

            if min(u_x)-1<1e-12 && oblate(k)==0
                fprintf("\n warning: u_x~1. check spheroid #%d\n ",k);
            end
            
            % % split up interior/surface/exterior
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
            SLregions=cell(3,1);
            
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
                        Yn=Ynm(n,[],real(acos(v_x_r)),phi_x_r); % v_x_r exceeds [-1,1] by 1e-8, but acos() returns imaginary values. Impose real values for Ynm.
                        Yr(:,n^2+1:(n+1)^2)=Yn;
                    end

                    % solid spheroidal harmonics
                    FYr=Fr.*Yr;

                    % Spheroidal harmonic coefficients of D[sigma]
                    spectra_matrix = repmat(spectra_regions{r},1,nf);
                    SLcoefs_r = spectra_matrix.*Gshc(:,:,k);
                    
                    % multiply spheroidal harmonic coefficients and
                    % harmonic functions together
                    SLregions{r}=FYr*SLcoefs_r;

                end
            end
    
            % Recombine all SLs from interior/exterior/surface
            SLk(u_x < u0(k)-1e-14,:) = SLregions{1};
            SLk(abs(u_x-u0(k))<=1e-14,:) = SLregions{2};
            SLk(u_x > u0(k)+1e-14,:) = SLregions{3};
            % SLk(u_x < u0(k),:) = SLregions{1};
            % SLk(u_x == u0(k),:) = SLregions{2};
            % SLk(u_x > u0(k),:) = SLregions{3};
    
            if isReal
                SLk=real(SLk);
            end
        end
           
        %evluation for particle k
        SL{k} = SLk;
        
    end

    if ~isa(X, "cell")
        SL = cell2mat(reshape(SL,1,1,ns));
    end

end
%--------------------------------------------------------------------%

end

function [lambda_int,lambda_surf,lambda_ext]=SLspectrum(p,u0,a,oblate)
    sp=(p+1)^2;
    ii = (1:sp)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;
    cnm=a.*factorial(nn-mm)./factorial(nn+mm).*(-1).^(mm).*(oblate.*1j.*sqrt(u0.^2+1)+~oblate.*sqrt(u0.^2-1));
    L_obl=legendre_otc(p,1j.*u0,1,1,1);
    u0_pro=u0+oblate; % if i-th u0 for oblate, will add 1 to make sure it is >1. actual value won't matter.
    L_pro=legendre_otc(p,u0_pro,1,1,1);
    Pp=L_pro{1}; Qp=L_pro{2}; 
    Po=L_obl{1}; Qo=L_obl{2}; 
    lambda_int=cnm.*(Qo*diag(oblate)+Qp*diag(~oblate));
    lambda_surf=cnm.*((Po.*Qo)*diag(oblate)+(Pp.*Qp)*diag(~oblate));
    lambda_ext=cnm.*(Po*diag(oblate)+Pp*diag(~oblate));   

    % if oblate
    %     cnm=1j.*a.*factorial(nn-mm)./factorial(nn+mm).*(-1).^(mm) .*sqrt(u0.^2+1);
    %     L=legendre_otc(p,1j.*u0,1);
    % else
    %     cnm=a.*factorial(nn-mm)./factorial(nn+mm).*(-1).^mm .*sqrt(u0.^2-1);
    %     L=legendre_otc(p,u0,1);
    % end
    % P=L{1}; Q=L{2};
    % 
    % lambda_int=cnm.*Q;
    % lambda_surf=cnm.*P.*Q;
    % lambda_ext=cnm.*P;   
end

function F=solid_harmonic(p,u0,u_x,oblate)
    F=1;
    if oblate
        % fprintf("oblate\n");
        u_x = 1j.*u_x;
    else
        % fprintf("prolate\n");
    end
    if abs(u_x)-u0 < -1e-14
        % try
        %     PQ=legendre_otc(p,u_x);
        % catch ME
        %     display(u_x);
        %     display(min(abs(u_x)));
        %     rethrow(ME);
        % end
        PQ=legendre_otc(p,u_x);
        P=PQ{1};
        F=P.';
    elseif abs(u_x) - u0 > 1e-14
        PQ=legendre_otc(p,u_x,1);
        Q=PQ{2};
        F=Q.';
    end
end

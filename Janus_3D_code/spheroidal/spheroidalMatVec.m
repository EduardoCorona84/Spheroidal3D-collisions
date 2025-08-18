function LP=spheroidalMatVec(params,potential,X,nu_vec)
%--------------------------------------------------------------------%
% spheroidalMatVec computes the Laplace layer potential of spheroids with 
% density 'sigma' evaluated either on the surfaces of all spheroids or at targets X. "far"
% particles are evaluated using smooth quadrature while near particles are 
% evaluated using spheroidal harmonics.
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
%       (o) sigma_coefficients = spherical harmonic coefficients of sigma
%       (o) matvec_eta = cutoff distance for near / far.
%    (o) potential = 'DL' for double layer, 'SL' for single layer
%    (o) X = target points at which to evaluate the potential. If not
%        provided, this function evaluates the potential at points on each
%        spheroid surface.
%    (o) nu_vec = the normal vectors associated with EACH target point (note
%        that this is only used for spheroidalSP and spheroidalDP).
%    
% Returns LP of size(sigma) for particle-to-particle evaluation or size [nrows(X),nf].
% 
%--------------------------------------------------------------------%

if isempty(params.sigma)
    error("No surface density given")
end
if isempty(params.u0)
    error("No surface parameter u_0 given")
end
if isempty(params.a)
    error("No scale (a) given")
end
if isempty(params.matvec_eta)
    error("No matvec_eta provided.");
end

sigma=params.sigma;
p=params.p;
if p < 2
    error('KernelEval (smooth quadrature) requires at least p=2.')
end
u0=params.u0;
a=params.a;
isReal=params.isReal;
[np,nf,ns]=size(sigma);
if_oblate = params.oblate;

if length(u0)==1
    u0=u0*ones(1,ns); 
end
if length(a)==1
    a=a*ones(1,ns); 
end

if nargin < 2
    error('Not enough input arguments. Must provide at least a SpheroidalParamters object and a potential.')

%particle-to-particle evaluation.
elseif nargin == 2 % Cannot be used for derivative terms (i.e. SP/DP)

    separation = params.separate_spheroids();
    
    % Do self evaluation of all particles
    % ------------------------------------------------------------
    if strcmp(potential,'DL')
        LP = spheroidalDL(params);
    elseif strcmp(potential, 'SL')
        LP = spheroidalSL(params);
    else
        error("invalid potential given. Use 'DL'/'SL' for double/single layer")
    end
    % ------------------------------------------------------------
    
    
    % Calculate effect from each particle on all the others
    for i=1:ns
    
        % Get target coordinates relative to self (particle i)
        [X,~]=params.get_X(i);
        ind=[1:i-1 i+1:ns];
        sep = separation(i,ind);
        sep_rep = repmat(sep, np,1);
        smooth_sep = sep_rep(:)==1; %indices where we separate smooth from spectral
    
        X_spectral = cell(1,ns);
        X_spectral{i} = X(~smooth_sep,:);
        X_smooth = X(smooth_sep,:);
    
        if strcmp(potential,'DL')
            LP_spectral_cell = spheroidalDL(params, X_spectral);
        elseif strcmp(potential, 'SL')
            LP_spectral_cell = spheroidalSL(params, X_spectral);
        end

        LP_smooth=[];
    
        if size(X_smooth,1) > 0 %at least one target is far
            fprintf("\n in smooth routine\n")
            % --------------------------------
            %  Do kernel eval for far targets
            
            %Kernel parameters
            if ~if_oblate
                Xself=prolate_spheroid_shape(p,u0(i),a(i));
            else
                Xself=oblate_spheroid_shape(p,u0(i),a(i));
            end
            Sns=SurfaceSph(Xself);
            
            pot=strcat(potential,'_L_3D');
            KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-8,2,400,1);
            KEparams.dim = 3;
            [~, gwt]=g_grid(p+1);
            wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
            wt = wt(:);
            Wns = Sns.geoProp.W; Wns= Wns.*wt;
            Nrns = reshape(Sns.geoProp.nor.to_array,[],3);
            KEparams.X = Xself;
            KEparams.nor = Nrns; 
            KEparams.W2 = Wns.';
            
            LP_smooth = Kernel_Eval(X_smooth,Xself,KEparams)*sigma(:,:,i);
        
        end

        LPi = zeros(np*(ns-1),nf);
        LPi(smooth_sep,:) = LP_smooth;
        LPi(~smooth_sep,:) = LP_spectral_cell{i};
        
        % reshape into matrix where each column is a spheroid (particle i
        % excluded)
        LPi=reshape(LPi,np,ns-1,nf);

        % Change LPi to shape np x nf x ns-1 to match with LP
        LPi=permute(LPi,[1 3 2]);

        % Add contribution from particle i on all other targets.
        LP(:,:,ind) = LP(:,:,ind) + LPi;
    end

elseif nargin > 2
    if (strcmp(potential, 'SP') || strcmp(potential, 'DP')) && nargin < 4
        error("No normal vector inputted.");
    end

    nt=size(X,1);

    [sep,Xt]=params.separate_targets(X);
    X_spectral = cell(1,ns);
    if strcmp(potential,'SP') || strcmp(potential, 'DP')
        nu_spectral = cell(1,ns);
        Nu_t = params.get_nu_targets(nu_vec);
    end
    
    % Separate target points and normal vectors into nearby and far
    for i=1:ns
        is_near = (sep(:,i) == 0);
        X_spectral{i} = Xt(is_near,:,i);
        if strcmp(potential, 'SP') || strcmp(potential, 'DP')
            nu_spectral{i} = Nu_t(is_near,:,i);
        end
    end

    %%%
    %%% Do "nearby" evaluation with spheroidal harmonics
    %%%
    if strcmp(potential,'DL')
        LP_spectral_cell = spheroidalDL(params,X_spectral);
    elseif strcmp(potential,'SL')
        LP_spectral_cell = spheroidalSL(params,X_spectral);
    elseif strcmp(potential,'SP')
        LP_spectral_cell = spheroidalSP(params, X_spectral, nu_spectral);
    elseif strcmp(potential, 'DP')
        LP_spectral_cell = spheroidalDP(params, X_spectral, nu_spectral);
    else
        error("Invalid potential given");
    end
    
    %%%
    %%% Do "far" evaluation with smooth quadrature
    %%%
    LP=zeros(nt,nf);
    
    % Calculate effect from each particle on each target point
    for i=1:ns

        % Add contributions from spectral method
        LP(sep(:,i)==0,:) = LP(sep(:,i)==0,:) + LP_spectral_cell{i}; 
    
        % Get targets for smooth quadrature
        X_smooth = Xt(sep(:,i)==1,:,i);
        
        if size(X_smooth,1) > 0 %at least one target is far
            % --------------------------------
            %  Do kernel eval for far targets
            
            %Kernel parameters
            if ~if_oblate(i)
                Xself=prolate_spheroid_shape(p,u0(i),a(i));
            else
                Xself=oblate_spheroid_shape(p,u0(i),a(i));
            end
            Sns=SurfaceSph(Xself);
            
            if strcmp(potential,'SP')
                pot='dSL_L_3D';
            elseif strcmp(potential, 'DP')
                pot='dDL_L_3D';
            else
                pot=strcat(potential,'_L_3D');
            end
            KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-8,2,400,1);
            KEparams.dim = 3;
            [~, gwt]=g_grid(p+1);
            wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
            wt = wt(:);
            Wns = Sns.geoProp.W; Wns= Wns.*wt;

            %%% Normal vector handling
            if strcmp(potential,'SP')
                target_norm_vecs = Nu_t(sep(:,i)==1,:,i);
                KEparams.nor = target_norm_vecs; 
            elseif strcmp(potential, 'DP')
                source_norm_vecs = reshape(Sns.geoProp.nor.to_array,[],3);
                target_norm_vecs = Nu_t(sep(:,i)==1,:);
                KEparams.nor = source_norm_vecs;
                KEparams.targnor = target_norm_vecs;
            else
                Nrns = reshape(Sns.geoProp.nor.to_array,[],3);
                KEparams.nor = Nrns;
            end
            KEparams.X = Xself;
            KEparams.W2 = Wns.';
            
            LP_smooth = Kernel_Eval(X_smooth,Xself,KEparams)*sigma(:,:,i);
    
            LP(sep(:,i)==1,:) = LP(sep(:,i)==1,:) + LP_smooth;
        end 
    end
end

if isReal
    LP = real(LP);
end

end

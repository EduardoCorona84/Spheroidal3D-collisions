function M = spheroidalMatVecKernel(params,potential,p,gradSP)

    if isempty(params.u0)
        error("No surface parameter u_0 given")
    end
    if isempty(params.a)
        error("No scale (a) given")
    end
    if isempty(params.matvec_eta)
        error("No matvec_eta provided.");
    end
    
    if nargin < 3
        p=pars.p;
    end
    
    if p < 2
        error('KernelEval (smooth quadrature) requires at least p=2.')
    end

    u0=params.u0;
    a=params.a;
    thetas=params.thetas;
    phis = params.phis;
    centers = params.centers;
    if_oblate = params.oblate;

    ns = size(centers,1);
    np = 2*p*(p+1);

    if length(thetas)==1
        thetas = thetas*ones(1,ns);
    end
    if length(phis)==1
        phis = phis*ones(1,ns);
    end
    if length(u0)==1
        u0=u0*ones(1,ns); 
    end
    if length(a)==1
        a=a*ones(1,ns); 
    end

    separation = params.separate_spheroids();

    % Do self evaluation of all particles (matrices on main block diagonal)
    % ------------------------------------------------------------
    IDparams = copy(params);
    IDparams.sigma = repmat(eye(np),1,1,ns);

    if strcmp(potential,'DL')
        LP = spheroidalDL(IDparams);
    elseif strcmp(potential, 'SL')
        LP = spheroidalSL(IDparams);
    elseif strcmp(potential, 'SP')
        LP = spheroidalSP(IDparams);
    else
        error("invalid potential given. Use 'DL'/'SL' for double/single layer")
    end
    % ------------------------------------------------------------
    if nargin==3
        LP = mat2cell(LP,np,np,ones(1,ns));
    else
        LP = mat2cell(LP,3*np,np,ones(1,ns));
    end
    M = blkdiag(LP{:});

    % Calculate effect from each particle on all the others
    if ns~=1
        for i=1:ns
            % Get target coordinates relative to self (particle i)
            [X,~]=IDparams.get_X(i);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get Nu from each target X.
            % Nu=zeros(size(X));
            % figure;
            [Nu,~]=IDparams.get_Norm_rot(p,i);
            % quiver3(X(:,1),X(:,2),X(:,3),Nu(:,1),Nu(:,2),Nu(:,3)); hold on;
            % for j=1:ns-1
            %     xind=(j-1)*np+1;
            %     quiver3(X(xind:xind+np-1,1),X(xind:xind+np-1,2),X(xind:xind+np-1,3),Nu(xind:xind+np-1,1),Nu(xind:xind+np-1,2),Nu(xind:xind+np-1,3)); hold on;
            %     surfsph=SurfaceSph(X(xind:xind+np-1,:));
            %     Nu(xind:xind+np-1,:) = reshape(surfsph.geoProp.nor.to_array,[],3);
            %     quiver3(X(xind:xind+np-1,1),X(xind:xind+np-1,2),X(xind:xind+np-1,3),Nu(xind:xind+np-1,1),Nu(xind:xind+np-1,2),Nu(xind:xind+np-1,3)); hold on;
            % end
            % hold off;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ind=[1:i-1 i+1:ns];
            sep = separation(i,ind);
            sep_rep = repmat(sep, np,1);
            smooth_sep = sep_rep(:)==1; %indices where we separate smooth from spectral

            X_spectral = cell(1,ns);
            X_spectral{i} = X(~smooth_sep,:);
            X_smooth = X(smooth_sep,:);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Nu_spectral = cell(1,ns);
            Nu_spectral{i} = Nu(~smooth_sep,:);
            Nu_smooth= Nu(smooth_sep,:);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            if strcmp(potential,'DL')
                LP_spectral_cell = spheroidalDL(IDparams, X_spectral);
            elseif strcmp(potential, 'SL')
                LP_spectral_cell = spheroidalSL(IDparams, X_spectral);
            elseif strcmp(potential, 'SP')
                LP_spectral_cell = spheroidalSP(IDparams, X_spectral,Nu_spectral);
            end
    
            LP_smooth=[];
        
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
                
                if strcmp(potential, 'SP')
                    pot='dSL_L_3D';
                else
                    pot=strcat(potential,'_L_3D');
                end
                KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-8,2,400,1);
                KEparams.dim = 3;
                [~, gwt]=g_grid(p+1);
                wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
                wt = wt(:);
                Wns = Sns.geoProp.W; Wns= Wns.*wt;
                if strcmp(potential,'DL')
                    Nrns = reshape(Sns.geoProp.nor.to_array,[],3);
                    KEparams.nor = Nrns; 
                elseif strcmp(potential,'SP')
                    % Strg1=SurfaceSph(X_smooth);
                    % Nrns = reshape(Strg1.geoProp.nor.to_array,[],3);
                    Nrns = Nu_smooth;
                    KEparams.nor = Nrns; 
                end
                KEparams.X = Xself;
                KEparams.W2 = Wns.';
                
                LP_smooth = Kernel_Eval(X_smooth,Xself,KEparams);
            
            end
    
            LPi = zeros(np*(ns-1),np);
            LPi(smooth_sep,:) = LP_smooth;
            LPi(~smooth_sep,:) = LP_spectral_cell{i};
    
            full_ind = ones(1,ns);
            full_ind(i)=0;
            full_ind = reshape(repmat(full_ind,np,1),[],1)==1;
    
            M(full_ind,(i-1)*np+1:i*np) = LPi;
            
        end
    end
end

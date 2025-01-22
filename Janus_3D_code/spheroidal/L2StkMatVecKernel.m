function M = L2StkMatVecKernel(pars,p)
    % create the Stokes single layer matrix from Laplace kernel matrices.

    if isempty(pars.u0)
        error("No surface parameter u_0 given")
    end
    if isempty(pars.a)
        error("No scale (a) given")
    end
    if isempty(pars.matvec_eta)
        error("No matvec_eta provided.");
    end
    
    if nargin < 2
        p=pars.p;
    end
    
    if p < 2
        error('KernelEval (smooth quadrature) requires at least p=2.')
    end

    u0=pars.u0;
    a=pars.a;
    centers = pars.centers;

    ns = size(centers,1);
    np = 2*p*(p+1);

    if length(u0)==1
        u0=u0*ones(1,ns); 
    end
    if length(a)==1
        a=a*ones(1,ns); 
    end

    separation = pars.separate_spheroids();

    % Do self evaluation of all particles (matrices on main block diagonal)
    % ------------------------------------------------------------
    IDparams = copy(pars);
    IDparams.sigma = repmat(eye(np),1,1,ns);

    % SL: np x np; 
    % M1 for each spheroid: 3np x 3np
    % [ SL 0  0 ][ sigma_x ] = [ SL[sigma_x] ]
    % [ 0  SL 0 ][ sigma_y ]   [ SL[sigma_y] ]
    % [ 0  0  SL][ sigma_z ]   [ SL[sigma_z] ] 
    SLM=spheroidalSL(IDparams);
    M1_mat = blkdiag(SLM,SLM,SLM);
    M1_mat = repmat(M1_mat,1,1,ns);
    M1_cells = mat2cell(M1_mat,3*np,3*np,ones(1,ns));

    % Surface SP but with arbitrary derivative vector, 
    % dSx: np x np; 
    % M2 for each spheroid: 3np x 3np
    % [ x*dSx y*dSx z*dSx ][ sigma_x ]
    % [ x*dSy y*dSy z*dSy ][ sigma_y ]
    % [ x*dSz y*dSz z*dSz ][ sigma_z ]

    % Form SPM = [ dSx ; dSy ; dSz ]
    SPM=zeros(np*3,np,ns);
    nu_x=repmat([1,0,0],np,1,ns);
    nu_y=repmat([0,1,0],np,1,ns);
    nu_z=repmat([0,0,1],np,1,ns);
    for ii=1:ns
        % Create separate system for each spheroid for self-to-self dS with
        % arbitrary normal to avoid calculating all to self effect.

        % LOCparams=SpheroidalParameters;
        % LOCparams.u0=IDparams.u0(ii);
        % LOCparams.a=IDparams.a(ii);
        % LOCparams.oblate=IDparams.oblate(ii);
        % LOCparams.thetas=IDparams.thetas(ii);
        % LOCparams.phis=IDparams.phis(ii);
        % LOCparams.centers=IDparams.centers(ii,:);
        % LOCparams.sigma=IDparams.sigma(:,:,ii);
        % LOCparams.get_shc();
        % X_self2 = LOCparams.get_X();
        % 
        % [LP_x,LP_y,LP_z] = spheroidalSP(LOCparams,X_self2,nu_x,nu_y,nu_z);
        [LP_x,LP_y,LP_z] = spheroidalSP(IDparams,[],nu_x,nu_y,nu_z);

        SPM((ii-1)*np+1:ii*np,:,ii)=LP_x;
        SPM(np+(ii-1)*np+1:np+ii*np,:,ii)=LP_y;
        SPM(np*2+(ii-1)*np+1:np*2+ii*np,:,ii)=LP_z;
    end

    % Multiply by [x y z]; target points self eval same as self points.
    [Xtrg,~]=IDparams.get_X();
    M2_cells=cell(1,ns);
    for ii=1:ns
        if isa(Xtrg,"cell")
            Xtrg_ii=Xtrg{ii};
        else
            Xtrg_ii=Xtrg(:,:,ii);
        end
        M2_cells{ii}=[repmat(Xtrg_ii(:,1),3,1).*SPM(:,:,ii),repmat(Xtrg(:,2),3,1).*SPM(:,:,ii),repmat(Xtrg(:,3),3,1).*SPM(:,:,ii)];
    end

    % SPM[Xsrc dot sigma]
    % dSx: np x np;
    % Xsrc_x: np x np;
    % [ dSx ]                        [ sigma_x ]
    % [ dSy ][ Xsrc_x Xsrc_y Xsrc_z ][ sigma_y ]
    % [ dSz ]                        [ sigma_z ]

    M3_cells=cell(1,ns);
    for ii=1:ns
        SPii=SPM(:,:,ii); % 3np x np
        if ~pars.oblate(ii)
            Xloc=prolate_spheroid_shape(pars.p,pars.u0(ii),pars.a(ii));
        else
            Xloc=oblate_spheroid_shape(pars.p,pars.u0(ii),pars.a(ii));
        end
        M3_cells{ii}=SPii*[diag(Xloc(:,1)),diag(Xloc(:,2)),diag(Xloc(:,3))];
    end

    M=blkdiag(M1_cells{:}) - blkdiag(M2_cells{:}) + blkdiag(M3_cells{:});
    M=1/2.*M;
    

    % Calculate effect from each particle on all the others
    if ns~=1
        for i=1:ns
            % Get target coordinates relative to self (particle i)
            [X,~]=IDparams.get_X(i);
            Ri=IDparams.Rmats(:,:,i);
            Nu_x_ii=(Ri * nu_x')';
            Nu_y_ii=(Ri * nu_y')';
            Nu_z_ii=(Ri * nu_z')';
            ind=[1:i-1 i+1:ns];
            sep = separation(i,ind);
            sep_rep = repmat(sep, np,1);
            smooth_sep = sep_rep(:)==1; %indices where we separate smooth from spectral

            X_spectral = cell(1,ns);
            X_spectral{i} = X(~smooth_sep,:);
            X_smooth = X(smooth_sep,:);
            Nu_spectral_x = cell(1,ns); Nu_spectral_y=cell(1,ns); Nu_spectral_z=cell(1,ns);
            Nu_spectral_x{i} = Nu_x_ii(~smooth_sep,:);
            Nu_spectral_y{i} = Nu_y_ii(~smooth_sep,:);
            Nu_spectral_z{i} = Nu_z_ii(~smooth_sep,:);

            SL_spectral_cell=spheroidalSL(IDparams,X_spectral);
            [SP_spectral_x_cell,SP_spectral_y_cell,SP_spectral_z_cell]=spheroidalSP(IDparams,X_spectral,Nu_spectral_x,Nu_spectral_y,Nu_spectral_z);

            % TODO: make SL_stk spectral matrix
            LP_spectral_cell=cell(1,ns);

            % ---------------------------------------------

            LP_smooth=[];
            if size(X_smooth)~=0
                if ~if_oblate(i)
                    Xself=prolate_spheroid_shape(p,u0(i),a(i));
                else
                    Xself=oblate_spheroid_shape(p,u0(i),a(i));
                end
                Sns=SurfaceSph(Xself);
    
                KEparams = Kernel_Eval_parameters('SL_Stk_3D',0,1,1,1,1e-8,2,400,1);
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
    
                LP_smooth = Kernel_Eval(Xtrg_ii,Xv,KEparams);

                % TODO: Need to either permute LP_Kernel for [sig_x;sig_y;sig_z]
                % input, or change Xv,Wv to make new kernel, and perform
                % tests again. 
                
            end

            % TODO: change below for dimension matching.
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
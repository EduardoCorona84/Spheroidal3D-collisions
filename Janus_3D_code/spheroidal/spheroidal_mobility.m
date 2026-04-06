function spheroidal_mobility(fname,Fparams,init)
    %{
    Entrypoint for the mobility solver for spheroidal suspension in Stokes flow.
    Code was adapted from the code for Janus particles.
     
    Inputs: 
    fname - (string) filename for experiment info
    
    Fparams - (struct) parameter struct for rigid body simulation, with fields:
    ----
    type            - (string) problem type (FTfun)
    typeMV          - (string) 'SSph' for (scalar) spheroidal harmonics, 'Rbs' for rotation based singular quad
    denseMV         - (bool) dense vs FMM for far-field
    num_timesteps   - (int)    number of timesteps
    dt              - (double) timestep length
    tdisc           - (string) time discretization 
        Implemented are "euler", "trapz", and "rk4" (and also "abash"?)
    comp            - (bool) compute intermediate quantities FT and VW (for debugging purposes?)

    parbd - (struct) struct with rigid body parameters:
        n3          - (int) number of rigid bodies n_b
        equ_radii   - (double) n_b x 1 array of equatorial radii for spheres/spheroids; advised to be set to 1 for efficiency purposes.
        polar_radii - (double) n_b x 1 array of polar radii for spheres/spheroids; obviously, this should be less than major_radii
            If one truly wants the u0/a values for the spheroids instead, see the utility functions -Brian.
            The convention for the codebase is that if equi_radius > polar_radius, it is an oblate. Otherwise,
            it's a prolate.
        shape_type  - (string) n_b x 1 array of strings: should be 'prolate', 'oblate', or 'sphere'.
            Spheres are not implemented for the time being.
        p           - (int) spheroidal harmonic order (bodies) 
        Ct          - (double) n_b x 3 array of centers 
        collision_eps - (double) collision buffer/tolerance
        mdist       - (double) collision buffer for body-body interactions
    
    parslv - (struct) linear solver parameters such as 
         prec     - (string) preconditioner type, '' for unprec, 'bkdiag'
                    (block diagonal), 'TT' (tensor train)
         precond_type - (string) 'bkdiag' or 'TT'
         solver   - 'gmres' (default)
         tol (tolerance), maxit (maximum restart cycles), rst (restart size), etc.
         Contact-LCP options:
            col_dense_max_pairs - when the number of frictionless contact
                        pairs is at most this value, explicitly build the
                        small projected contact matrix A = F^T M F rather
                        than using a matrix-free LCP operator.
            td_sparse_nn - when true, build a sparse contact-pair body-
                        block preconditioner for TD solves near or at
                        contact.
            td_sparse_nn_reg - diagonal regularization factor applied to
                        the sparse contact-pair TD operator before sparse
                        factorization.
            td_sparse_nn_droptol - ILU drop tolerance for the sparse
                        contact-pair TD preconditioner.
            td_sparse_nn_inverse - inverse method used on the sparse
                        contact-pair TD operator; should be 'ilu' or 'lu'.
            td_sparse_nn_reuse_local_state - when true, reuse the assembled
                        and inverted active-body local TD state within the
                        current timestep.
         Augmentation options:
            deflate - when true, near-contact/active-contact TD_main solves
                        attach a contact-aware coarse basis before the Krylov
                        solve.
            deflate_max_dim - maximum total dimension of that coarse basis.
            deflate_sval_tol - singular-value cutoff used when building the
                        contact-derived basis from each active contact block.
            deflate_near_contact - when true, use the nearest prospective
                        contact pairs to switch TD_main into the targeted
                        augmented solve before collision detection activates.
            deflate_near_contact_factor - near-contact switch threshold,
                        measured in multiples of collision_eps.
            deflate_probe_pairs - number of nearest candidate pairs checked
                        for that near-contact switch.
    
    Depending on the type of problem that is implemented
    Ffun, Tfun = @(t,C,q) with output of size 3 x n_b.
        Force and torque prescriptions; required if type is 'FTfun'.
        The parameters are: t for timestep, C for the center of the body, and
        q for the quaternion associated with the body (to represent orientation).
    
    init    - (string) optional filename to resume a simulation from last
    recorded timestep

    Plot options (Fparams):
        plotFlag          enable plotting (bool)
        plotEvery         update cadence in timesteps (int, >=1)
        plotView          fixed camera view [azimuth elevation] or 3-vector
        plotAxis          fixed axis limits [xmin xmax ymin ymax zmin zmax]
        plotTitle         title prefix (string)
        plotGrid          grid on/off (bool)
        plotColor         'sigma' | 'mu' | 'none' (default none)
        plotColorMode     'inf' (default) | 'l2'
        plotColorLimits   [cmin cmax]
        plotSurfaceAlpha  surface transparency in (0, 1)
        plotTrajectories  show center trajectories (bool)
        plotTrajMaxPoints cap trajectory length (int)
        plotTrajLineWidth trajectory line width (scalar)
        plotForceVectors  show per-body force vectors (bool, default false)
        plotForceColor    force vector RGB color ([r g b])
        plotVectorLineWidth overlay line width (default 1.5)
        plotForceMaxLength max force-vector length (default 0.8*diam)

    Optional background flow (Fparams.background_flow):
        enabled           enable background flow (bool, default false)
        U0                3-vector offset velocity in the lab frame
        A                 3x3 trace-free velocity-gradient matrix
    ----

    %%%
    %%% CODE ANNOTATIONS
    %%%
    Mt stores rotational information for n bodies at each time step
    Mt0 is the initial rotation setup

    VW represents the translational (v)/angular velocity (w) each timestep (as a 6 x num_body matrix)
    VW0 represents the initial velocities

    Xrp are the reference surface points used for the body discretization.

    Fparams.parbd.np is number of discretization points on the surface of a body (is the same for all bodies?)

    The code "xind = (1:np) + np*(k-1)" shows up often; it collects all of the discretization points for the body k.
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global timings;
    global ROTATIONAL_VELOCITY_TOL;     % The threshold that the norm of the angular velocity should meet in order
    ROTATIONAL_VELOCITY_TOL = 1e-10;    % to update the body's angular position.
     
    %%(0.0) Input validation
    if ~isempty(Fparams)
        % Generic property and parbd validation with defaults
        assert(isscalar(Fparams.Nt) && Fparams.Nt>=1,...
            'Fparams.Nt must be a positive integer.');
        assert(isfield(Fparams,'dt') && isnumeric(Fparams.dt) && isscalar(Fparams.dt) && Fparams.dt>0,...
            'Fparams.dt must be a strictly positive scalar.');
        assert(any(strcmp(Fparams.tdisc, {'euler','trapz','rk4','abash'})), ...
            'Fparams.tdisc must be one of: euler, trapz, rk4, abash.');
        assert(islogical(Fparams.denseMV) || ismember(Fparams.denseMV,[0,1]), ...
            'Fparams.denseMV must be true or false.');
        assert(islogical(Fparams.comp) || ismember(Fparams.comp,[0,1]), ...
            'Fparams.comp must be true or false.');

        % parbd
        req = {'p','Ct','equ_radii','polar_radii','mdist'};
        for k = 1:numel(req)
            assert(isfield(Fparams.parbd,req{k}), ...
                'Missing required field Fparams.parbd.%s', req{k});
        end

        assert(isnumeric(Fparams.parbd.Ct) && size(Fparams.parbd.Ct,2)==3, ...
            'Fparams.parbd.Ct must be an n3-by-3 double array.');
        assert(isscalar(Fparams.parbd.p) && Fparams.parbd.p>=2, ...
            'Fparams.parbd.p must be a positive integer greater than 2.');
        assert(numel(Fparams.parbd.equ_radii)==size(Fparams.parbd.Ct,1) && numel(Fparams.parbd.polar_radii)==size(Fparams.parbd.Ct,1), ...
            'Fparams.parbd.equ_radii and Fparams.parbd.polar_radii must have length equal to size(parbd.Ct,1).');
        assert(isscalar(Fparams.parbd.collision_eps)  && Fparams.parbd.collision_eps > 0, ...
            'Fparams.parbd.collision_eps must be a strictly positive scalar.');
        assert(isscalar(Fparams.parbd.mdist) && Fparams.parbd.mdist > 0, ...
            'Fparams.parbd.mdist must be a strictly positive scalar.');
    end

    if isfield(Fparams, 'plotFlag')
        plot_enabled = Fparams.plotFlag;
    else
        plot_enabled = false;
    end

    LOCAL_log_input_options(fname, Fparams, init);
    LOCAL_manage_td_sparse_nn_cache('reset', [], []);

    %%(0.1) (optional) Load data in init, initialize output arrays
    num_timesteps = Fparams.Nt; num_body=Fparams.parbd.n3; 
    tt = zeros(num_timesteps+1, 1); 
    sigma = cell(num_timesteps+1, 1); mu=sigma; U=sigma; VW=U; Xt=U; Ct=Xt; FT=mu; psi_Lap = mu; 
    Mt = cell(num_timesteps+1,num_body);
    Energy = zeros(1,num_timesteps);
    dt0 = Fparams.dt; 
    
    if isempty(init)
        init_flag = false; 
        
        % Evolution
        t=0; tt(1)=0;  
        Fparams.parslv.prev=[]; %initialize preconditioner params
        Mt0 = LOCAL_get_initial_rotations(Fparams, num_body);
        Mt(1,:) = Mt0;
    else % Start the simulation from a previously saved state
        init_flag = true;  
        
        % Load previous file 
        load(init);
        
        prev_timesteps = sum(tt>0);

        % Recover initial data at time t0 = tt(lid): 
        T0 = tt(prev_timesteps); Xt0 = Xt{prev_timesteps}; Ct0 = Ct{prev_timesteps}; Mt0 = Mt(prev_timesteps,:);
        VW0 = VW{prev_timesteps}; 

        Fparams.parslv.prec = []; 
        Fparams.parslv.prev = []; 

        % Re-initialize arrays (after load)
        tt = zeros(num_timesteps+1,1); 
        sigma = cell(num_timesteps+1,1); mu=sigma; U=sigma; VW=U; Xt=U; Ct=Xt; FT=mu; psi_Lap = mu;
        Mt = cell(num_timesteps+1,num_body);
        Energy = zeros(1,num_timesteps);

        t=T0; tt(1)=T0; 
        Xt{1}=Xt0; Ct{1}=Ct0; Mt(1,:) = Mt0; 
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %(0.2) Initialize timings and Fparams struct before simulation
    
    zN = zeros(num_timesteps,1); 

    timings = struct();

    timings.setup_surace    = 0;
    timings.setup_kernel    = 0;
    timings.incoming        = zN;

    timings.velocities          = struct();
    timings.velocities.solve    = zN;
    timings.velocities.apply    = zN;
    timings.velocities.vw       = zN;
    timings.velocities.col      = zN;
    timings.velocities.shell    = zN;
    timings.velocities.total    = zN;

    timings.advance = zN;

    timings.operator        = struct();
    timings.operator.surf   = zN;
    timings.operator.diag   = zN;
    timings.operator.offd   = zN;
    timings.operator.total  = zN;

    timings.total = zN;

    tic;
    Fparams = SpheroidalMS_initparams(Fparams);
    timings.setup_surace = toc;
    timedisc = Fparams.tdisc; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %(0.3) Initialize kernels (for MatVecs) and nullspace info
    tic; 
    Xref = Fparams.parbd.Xrp;
    Xrp = Xref;
    C0 = Fparams.parbd.C; 
    np = Fparams.parbd.np; normW = zeros(num_body,1);  
    if init_flag % Load state from memory
        for k=1:num_body
            normW(k) = norm(VW0(4:6, k));
        end 
    end
    
    % Get kernels/nullspace
    Kernels=[]; 
    [Kernels,Nullsp,Fparams,timings] = SpheroidalMS_UpdateOperators(Xrp,C0,Mt0,normW,Kernels,Fparams,timings,0); 
    timings.setup_kernel=toc;
    Xrp = Fparams.parbd.Xrp;
    C0 = Fparams.parbd.C;
    np = Fparams.parbd.np;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (0.4) Initialize collision info  
    [colevent, collist ,mindst, distances, closest_points_1, closest_points_2] = LOCAL_check_collision(C0, Mt0, Fparams);
    X2 = [];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (0.5) Initial state + plot setup
    if isempty(Xt{1})
        Xt{1} = Xrp;
        Ct{1} = Fparams.parbd.C;
    end

    plot_state = [];
    if Fparams.plotFlag
        plot_state = plot_init(Fparams, Xt{1}, np, num_body);
        if plot_state.traj_enable
            plot_state = plot_update_trajectory(plot_state, Ct{1});
        end
        plot_state = plot_update(plot_state, Xt{1}, Ct{1}, t, [], [], []);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (1.0) Timestepping (i.e. the main loop)
    
    % Linearize differential variational inequality each timestep
    for i=1:num_timesteps
        dt = dt0;
        if strcmp(timedisc,'euler')
            fprintf('\n ---------------------------------------------------------- \n')
            fprintf('\nExplicit euler step for timestep %d\n', i)
            [Xt{i+1},Mt(i+1,:),Ct{i+1},U{i},FT{i},sigma{i},mu{i},VW{i},Kernels,Nullsp,...
                Fparams,colevent,collist,closest_points_1,closest_points_2,dt,psi_Lap{i},Energy(i)] = ...
                LOCAL_euler_step(Xt{i},Xref,X2,Mt(i,:),Ct{i},Kernels,Nullsp,Fparams,colevent,collist, closest_points_1, closest_points_2, t, dt,i);
        elseif strcmp(timedisc,'trapz')
            error('Calls need to be updated.')
            % (1) Predictor step: 
            fprintf('\n ---------------------------------------------------------- \n')
            fprintf('\n (1) Trapezoidal, predictor step \n')
            [Xt1,Mt1,Ct1,U1,FT1,sigma1,mu1,VW1,Kernels,Nullsp,...
                Fparams,colevent,collist,dt] = ...
                LOCAL_euler_step(Xt{i},Xref,X2,Mt(i,:),Ct{i},Kernels,Nullsp,Fparams,colevent,collist,t,dt,i);
            
            % (2) Corrector step: 
            fprintf('\n ---------------------------------------------------------- \n')
            fprintf('\n (2) Trapezoidal, corrector step \n')
            % Get incoming force distribution: 
            tic; 
            [FT2,sigma2,VW2,psi_Lap,Energy] = LOCAL_get_incoming_Fc(Fparams,t+dt,dt,Kernels,Nullsp,Xt1,Sc); 
            fprintf('\n Time to compute incoming force: %e',toc);
            timings.incoming(i) = timings.incoming(i)+toc;  
            fprintf('\n Fluid Solve at time %.2f',t)
            tic; 
            [sigma2,mu2,U2,VW2] = LOCAL_compute_velocities(sigma2,VW2,Ct1,Kernels,Nullsp,Fparams,colevent,collist,i,dt); 
            timings.velocities.total(i) = timings.velocities.total(i) + timings.velocities.solve(i) + timings.velocities.apply(i) + timings.velocities.vw(i) + timings.velocities.col(i); 
            fprintf('\n Time to compute velocities / fluid solve: %e',timings.velocities.total(i)); 
            
            % Correct VW{i} as average of VW0 and VW1 (and associated quantities)
            VW{i}    = 0.5*(VW1+VW2); 
            sigma{i} = 0.5*(sigma1+sigma2); 
            mu{i}    = 0.5*(mu1+mu2);
            U{i}     = 0.5*(U1+U2); 
            FT{i}    = 0.5*(FT1+FT2); 
            
            [Xt{i+1},Mt(i+1,:),Ct{i+1},Kernels,Nullsp,Fparams,colevent,collist,dt] = ...
                LOCAL_advance_step(VW{i},mu{i},sigma{i},Xt{i},Xref,X2,Mt(i,:),Ct{i},Kernels,Fparams,dt,i);
            
        elseif strcmp(timedisc,'rk4')
            error('Calls need to be updated.')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
            % Explicit Runge-Kutta 4th order method
            % (1) First step, f1 = f(t0,u0)
            fprintf('\n ---------------------------------------------------------- \n')
            fprintf('\n (1) Runge-Kutta 4th order, t=t0, x=x(t0) \n')
            dt41 = 0.5*dt; 
            t41 = t;
            
            [Xt41,Mt41,Ct41,U41,FT41,sigma41,mu41,VW41,Ker41,Null41,...
                Fpar41,cev41,clst41,dt41] = ...
                LOCAL_euler_step(Xt{i},Xref,X2,Mt(i,:),Ct{i},Kernels,Nullsp,Fparams,colevent,collist,t41,dt41,i);
            
            % (2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('\n ---------------------------------------------------------- \n')
            fprintf('\n (2) Runge-Kutta 4th order, t_{1/2}=t0+(dt/2), x=x(t_0)+(dt/2)*v1 \n')
            t42 = t+dt41; 
            dt42 = 0.5*dt; 
            
            % Get incoming force distribution: 
            tic; 
            [FT42,sigma42,VW42,psi_Lap, Energy] = LOCAL_get_incoming_Fc(Fpar41,t42,dt42,Ker41,Null41,Xt41,Sc); 
            fprintf('\n Time to compute incoming force: %e',toc);
            timings.incoming(i) = timings.incoming(i)+toc;  
            fprintf('\n Fluid Solve at time %.2f',t)
            tic; 
            [sigma42,mu42,U42,VW42] = LOCAL_compute_velocities(sigma42,VW42,Ct41,Ker41,Null41,Fpar41,cev41,clst41,i,dt42); 
            timings.velocities.total(i) = timings.velocities.total(i) + timings.velocities.solve(i) + timings.velocities.apply(i) + timings.velocities.vw(i) + timings.velocities.col(i); 
            fprintf('\n Time to compute velocities / fluid solve: %e',timings.velocities.total(i));
            
            % Advance Xt42 ~ X(t) + (dt/2)*V42 
            [Xt42,Mt42,Ct42,Ker42,Null42,Fpar42,cev42,clst42,dt42] = ...
                LOCAL_advance_step(VW42,mu42,sigma42,Xt{i},Xref,X2,Mt(i,:),Ct{i},Kernels,Fparams,dt42,i);
            
            % (3) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('\n ---------------------------------------------------------- \n')
            fprintf('\n (3) Runge-Kutta 4th order, t_{1/2}=t0+0.5*dt, x=x(t_{1/2}) \n')
            t43 = t+dt42; 
            dt43 = dt; 
            
            % Get incoming force distribution: 
            tic; 
            [FT43,sigma43,VW43] = LOCAL_get_incoming_Fc(Fpar42,t43,dt43,Ker42,Null42,Xt42,Sc); 
            fprintf('\n Time to compute incoming force: %e',toc);
            timings.incoming(i) = timings.incoming(i)+toc;  
            fprintf('\n Fluid Solve at time %.2f',t)
            tic; 
            [sigma43,mu43,U43,VW43] = LOCAL_compute_velocities(sigma43,VW43,Ct42,Ker42,Null42,Fpar42,cev42,clst42,i,dt43); 
            timings.velocities.total(i) = timings.velocities.total(i) + timings.velocities.solve(i) + timings.velocities.apply(i) + timings.velocities.vw(i) + timings.velocities.col(i); 
            fprintf('\n Time to compute velocities / fluid solve: %e',timings.velocities.total(i));
            
            % Advance Xt43 ~ X(t) + (dt)*V43 
            [Xt43,Mt43,Ct43,Ker43,Null43,Fpar43,cev43,clst43,dt43] = ...
                LOCAL_advance_step(VW43,mu43,sigma43,Xt{i},Xref,X2,Mt(i,:),Ct{i},Kernels,Fparams,dt43,i);
            
            % (4) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('\n ---------------------------------------------------------- \n')
            fprintf('\n (4) Runge-Kutta 4th order, t_{1}=t0+dt, x=x(t_{1}) \n')
            t44 = t+dt43; 
            dt44 = dt; 
            % Get incoming force distribution: 
            tic; 
            [FT44,sigma44,VW44] = LOCAL_get_incoming_Fc(Fpar43,t44,dt44,Ker43,Null43,Xt43,Sc); 
            timings.incoming(i) = timings.incoming(i)+toc;
            fprintf('\n Fluid Solve at time %.2f \n',t)
            [sigma44,mu44,U44,VW44] = LOCAL_compute_velocities(sigma44,VW44,Ct43,Ker43,Null43,Fpar43,cev43,clst43,i,dt44); 
            
            timings.velocities.total(i) = timings.velocities.total(i) + timings.velocities.solve(i) + timings.velocities.apply(i) + timings.velocities.vw(i) + timings.velocities.col(i); 
            fprintf('\n Time to compute velocities / fluid solve %.2f \n',timings.velocities.total(i));
            
            % Average velocities and related quantities using RK4 weights
            VW{i}    = (1/6)*(VW41+2*VW42+2*VW43+VW44); 
            sigma{i} = (1/6)*(sigma41+2*sigma42+2*sigma43+sigma44); 
            mu{i}    = (1/6)*(mu41+2*mu42+2*mu43+mu44); 
            U{i}     = (1/6)*(U41+2*U42+2*U43+U44); 
            FT{i}    = (1/6)*(FT41+2*FT42+2*FT43+FT44);
            
            % Advance Xtp ~ X(t) + (dt/6)*(V41 + 2*V42 + 2*V43 + V44) 
            [Xt{i+1},Mt(i+1,:),Ct{i+1},Kernels,Nullsp,Fparams,colevent,collist,dt] = ...
                LOCAL_advance_step(VW{i},mu{i},sigma{i},Xt{i},Xref,X2,Mt(i,:),Ct{i},Kernels,Fparams,dt,i);
        end
        
        timings.total(i) = timings.velocities.total(i) + timings.operator.total(i) + ...
            timings.advance(i) + timings.incoming(i); 
        fprintf('\n Total computing time for timestep %d : %e ',i,timings.total(i))
        fprintf('\n -------------------------------------------------------------\n'); 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        t = t+dt;
        tt(i+1)=t;
        fprintf('\n dt: %2.2f\n',dt)

        if plot_enabled
            if plot_state.traj_enable
                plot_state = plot_update_trajectory(plot_state, Ct{i+1});
            end
            plot_state = plot_update(plot_state, Xt{i+1}, Ct{i+1}, t, sigma{i}, mu{i}, FT{i});
        end

        % Save progress every other step.
        save_data = LOCAL_save_data(tt, Xt, Mt, Ct, FT, sigma, mu, U, VW, psi_Lap, Energy);
        if mod(i,2)==1                                                                                                                                              
            save(fname,'-v7.3','-struct','save_data');
        else                                                                                                                                                        
            save([fname '2'],'-v7.3','-struct','save_data');
        end  
        save([fname '_profile'],'timings'); 
    end %% END for linearization
end

function Mt0 = LOCAL_get_initial_rotations(Fparams, num_body)
    Mt0 = repmat({eye(3)}, 1, num_body);

    if ~isfield(Fparams, 'initial_rotation_matrices') || isempty(Fparams.initial_rotation_matrices)
        return;
    end

    raw_rotations = Fparams.initial_rotation_matrices;
    if iscell(raw_rotations)
        if numel(raw_rotations) ~= num_body
            error('Fparams.initial_rotation_matrices must provide one 3x3 matrix per body.');
        end
        Mt0 = reshape(raw_rotations, 1, []);
    elseif isnumeric(raw_rotations) && isequal(size(raw_rotations), [3 3 num_body])
        for body_idx = 1:num_body
            Mt0{body_idx} = raw_rotations(:,:,body_idx);
        end
    elseif isnumeric(raw_rotations) && isequal(size(raw_rotations), [3 3]) && num_body == 1
        Mt0 = {raw_rotations};
    else
        error('Fparams.initial_rotation_matrices must be empty, a cell array, or a 3x3xn3 numeric array.');
    end
end
    
function data = LOCAL_save_data(tt, Xt, Mt, Ct, FT, sigma, mu, U, VW, psi_Lap, Energy)
    data = struct();
    data.tt = tt;
    data.Xt = Xt;
    data.Mt = Mt;
    data.Ct = Ct;
    data.FT = FT;
    data.sigma = sigma;
    data.mu = mu;
    data.U = U;
    data.VW = VW;
    data.psi_Lap = psi_Lap;
    data.Energy = Energy;
end

%% Mobility solver system code
function [Xtp,Mtp,Ctp,U,FT,sigma,mu,VW,Kernels,Nullsp,Fparams,colevent,collist,closest_points_1,closest_points_2,dt,psi_Lap,Energy] = ...
    LOCAL_euler_step(Xt,X0,X2,Mt,Ct,Kernels,Nullsp,Fparams,colevent,collist, closest_points_1, closest_points_2, t,dt,it)
    %{
    Performs a single forward-Euler step of the system of rigid-body particles.
    Uses the BIE operators and boundary information at time t to compute surface
    velocities and rigid-body motions, advances centers and orientations,
    updates operators for the new geometry at the next timestep.

    Inputs
    Xt       - (double np*n3 x 3)
        current surface points at time t (rotated)
    X0       - (double np*n3 x 3)
        reference (unrotated) surface points
    X2       - (double np2*n3 x 3)
        model surface points for collision checks?
        (may be empty/unused...?)
    Mt       - (1-by-n3 cell)
        array of 3x3 rotation matrices at time t
    Ct       - (double n3 x 3)
        centers of bodies at time t
    Kernels  - struct
        BIE operator MatVecs (e.g., TD, SD and diagonals)
    Nullsp   - struct
        nullspace operators (C, B, D, L)
    Fparams  - struct
        simulation parameters
    colevent - boolean
        collision status entering the step
    collist  - (integer k x 2)
        list of candidate colliding body index pairs
    t        - double
        current simulation time
    dt       - double
        timestep size
    it       - integer
        timestep index used for timings bookkeeping

    Outputs
    Xtp      - (double np*n3 x 3)
        surface points at time t+dt
    Mtp      - (1-by-n3 cell)
        array of 3x3 rotation matrices at time t+dt
    Ctp      - (double n3 x 3)
        centers at time t+dt
    U        - (double 3*np*n3 x 1)
        stacked surface velocity vector (interleaved)
    FT       - (double 6*n3 x 1)
        total forces/torques at time t
    sigma    - (double 3*np*n3 x 1)
        incident traction density (particular)
    mu       - (double 3*np*n3 x 1)
        scattered density from Fredholm solve
    VW       - (double 6 x n3)
        rigid-body velocities; 1:3 translational V, 4:6 angular W per body
    Kernels  - struct
        BIE operators updated for geometry at t+dt
    Nullsp   - struct
        nullspace operators updated for geometry at t+dt
    Fparams  - struct
        simulation parameters (with updated point clouds (?))
    colevent - boolean
        collision status after advancing centers
    collist  - (integer k x 2)
        updated candidate collision pairs
    dt       - double
        step size
    psi_Lap  - placeholder
    Energy   - placeholder
    %}

    psi_Lap = NaN; % ugly placeholder
    
    global timings; 
    %Sc = Fparams.parbd.Sc; 
    np = Fparams.parbd.np; 
    n3 = Fparams.parbd.n3; 
    
    MRot = @(wh,t) RotationMat(wh,t);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get incoming force distribution.
    tic;
    [FT,sigma,VW,Energy] = LOCAL_get_incoming_Fc(Fparams,t,dt,Kernels,Nullsp,Xt);
    fprintf('\n Time to compute incoming force: %e ',toc)
    timings.incoming(it) = timings.incoming(it) + toc;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\n Fluid Solve at time %.2f ',t)
    tic;
    [sigma,mu,U,VW] = LOCAL_compute_velocities(sigma,VW,Ct,Kernels,Nullsp,Fparams,colevent,collist,it,dt, Mt, closest_points_1, closest_points_2);
    timings.velocities.total(it) = timings.velocities.total(it) + timings.velocities.solve(it) + timings.velocities.apply(it) ...
        + timings.velocities.vw(it) + timings.velocities.col(it);
    fprintf('\n Time to compute velocities / fluid solve: %e',timings.velocities.total(it));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Advance center Ct
    tic; 
    Ctp = LOCAL_advance_center(Ct,dt,VW); 
    fprintf('\n Time to advance centers C(t): %e',toc); 
    timings.advance(it) = timings.advance(it) + toc; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Advance rotation matrix Mt and X
    tic; 
    [Mtp,Xtp,normW] = LOCAL_advance_rotation(MRot,VW,Mt,Xt,X0,dt,np,n3); 
    fprintf('\n Time to advance R(t) and X(t): %e',toc)
    timings.advance(it) = timings.advance(it) + toc;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check for collision after moving centers
    tic; 
    [Ctp, Mtp, Xtp, normW, dt, colevent, collist, closest_points_1, closest_points_2] ...
    = LOCAL_collision_info(VW, Ct, Ctp, Mt, Mtp, MRot, Xt, Xtp, X0, dt, normW, Fparams);
    fprintf('\n Time for collision detection: %e',toc);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Update operators
    fprintf('\nOperator update\n')
    [Kernels,Nullsp,Fparams,timings] = SpheroidalMS_UpdateOperators(Xtp,Ctp,Mtp,normW,Kernels,Fparams,timings,it); 
    timings.operator.total(it) = timings.operator.total(it) + timings.operator.surf(it) + timings.operator.diag(it) + timings.operator.offd(it);
    fprintf('\n Time to update surface and operators: %e',timings.operator.total(it));
end
    
function [FT, fM, VW, Energy] = LOCAL_get_incoming_Fc(Fparams,t,dt,Kernels,Nullsp,Xt)
    %{
    This function is used to calculate the data on the surfaces for the solve.
    %}
    Energy = 0;
    parslv = Fparams.parslv; 

    n3=Fparams.parbd.n3; p = Fparams.parbd.p; np = Fparams.parbd.np; 
    VW=[]; C = Fparams.parbd.C; 
    Bk = Nullsp.B; Ck = Nullsp.C;
    
    switch Fparams.type
    case 'FTfun' % Force and torque are given
        Ffun = Fparams.Ffun; 
        Tfun = Fparams.Tfun; 
        Ct = Fparams.parbd.C; 
        
        Force = Ffun(t,Ct);
        Torque = Tfun(t,Ct);

        % Rebuild [F_1; T_1; F_2; T_2; ...], matching Build_SpheroidalAuxMats.
        FT = reshape([Force; Torque], [], 1);
        fM = Bk'*FT;
    case 'MHD'
        KLD = Kernels.KLD; 
        SLD = Kernels.SLD; 
        fprintf('\n Magnetic potential Solve at time %.2f \n',t)

        %% Magnetic solve for potential \phi.
        % Build RHS (i.e. \eta H_0 \cdot n)
        H0 = Fparams.H0(:);
        Nr_all = Fparams.parbd.Nrp; % Normals are calculated in set_params
        rhs = Fparams.eta * (Nr_all * H0);  % (np*n3) x 1
        
        % Solve
        q_density = Lslv(KLD,rhs,parslv); 
        fprintf('\n Magnetic potential solve res = %1.4g \n', norm(Lapp(KLD,q_density)-rhs)); 

        %% Compute Maxwell stress, forces and torques

        % phi at Gamma (Continuous)
        phi = -Xt*H0 + Lapp(SLD,q_density); 

        % -1*\nabla \phi = -\nabla_{\Gamma} phi - phi_n n
        GradientOfPotential = @(phi, phi_n,S) -1*S.geoProp.Grad(phi) -1*vec3d([phi_n; phi_n; phi_n]).*S.geoProp.nor;
        
        % Maxwell stress dotted with normal: n \cdot (E \oprod E - 1/2 |E|^2 I)
        maxwell_traction = @(E,S) times(dot(E, S.geoProp.nor), E) - times(dot(E,E), S.geoProp.nor)/2;

        % dphi/dn at Gamma (this is derived from continuity at the interface)
        % Fparams.mur = \mu / \mu_0 (dimensionless quantity)
        phi_n_e = Fparams.mur/(1-Fparams.mur)*q_density; 
        phi_n_i = 1/(1-Fparams.mur)*q_density; 

        % Magnetic Field (exterior and interior)
        fM = zeros(3*np*n3,1); 
        H_i = cell(n3,1); H_e=H_i; 
        for j=1:n3
            indx=(1:np)+np*(j-1); 
            indv=(1:3*np)+3*np*(j-1); 

            % Rebuild surface (since we need access to geoProp.Grad)
            S = SurfaceSph(vec3d(Xt(indx,:)));

            phi_body = phi(indx);
            phi_body = phi_body(:);
            phi_n_i_body = phi_n_i(indx);
            phi_n_i_body = phi_n_i_body(:);
            phi_n_e_body = phi_n_e(indx);
            phi_n_e_body = phi_n_e_body(:);

            H_i{j} = GradientOfPotential(phi_body, phi_n_i_body, S);
            H_e{j} = GradientOfPotential(phi_body, phi_n_e_body, S);

            %Maxwell stress . normal (traction)
            traction = maxwell_traction(H_e{j}, S) - maxwell_traction(H_i{j}, S); 
            traction = real(reshape(traction.to_array,[],3))'; 
            fM(indv) = traction(:); 
        end

        fprintf('\n Magnetic forces and torques \n'); 
        FT = real(Ck*fM); 
        display(reshape(FT,6,n3))
    case 'JanusAmp'
        DLMODD=Kernels.DLMODD; dDLMODD=Kernels.dDLMODD;
        flabel=Fparams.SurfaceLabel;
        % Solves for density, psi & uses them to compute normal derivative
        if isa(DLMODD,'function_handle')
            K = @(V) DLMODD(V);
            [psi,~,rel_residual,I]=gmres(K, flabel,100,1e-6);
            fprintf('\n Janus Amph D+0.5I BIE solve error = %e', rel_residual); 
            
            phi = K(psi); 
            phi_n_e = dDLMODD(psi);
        else
            K = DLMODD;
            [psi,~,rel_residual,I]=gmres(K, flabel,100,1e-6);
            fprintf('\n Janus Amph D+0.5I BIE solve error = %e', rel_residual); 
            phi = K*psi;
            phi_n_e = dDLMODD*psi;
        end

        % Computes energy (not necessary for dynamics)
        W = Fparams.parbd.Wg(:);
        Energy = real(-W' * (phi.*phi_n_e));
        
        phisquared = phi.*phi;

        % \nabla phi = \nabla_{\Gamma} phi + phi_n n
        GradientOfPotential = @(phi, phi_n, S) S.geoProp.Grad(phi) + vec3d([phi_n; phi_n; phi_n]).*S.geoProp.nor;
        
        % -n \cdot (2 * E \oprod E - |E|^2 I) + lambda*u^2
        % This is Equation (47) in the paper (there are typos in the paper, so instead refer to the original paper it's referencing).
        % Recall that \lambda^2 = 1/p^2.
        % Fparams.gamma is the ratio of the amphiphilic to viscous pressure (denoted eta in the paper).
        maxwellSnor = @(E,S,phi2) (Fparams.gamma) * ( ...
                            times(dot(E,E), S.geoProp.nor) - 2*times(dot(E, S.geoProp.nor), E) ...
                            + (Fparams.lambda)*times(phi2,S.geoProp.nor) ...
                        );
        fM = zeros(3*np*n3,1); 
 
        for j=1:n3
            indx=(1:np)+np*(j-1); 
            indv=(1:3*np)+3*np*(j-1); 

            S = SurfaceSph(vec3d(Xt(indx,:)));
            H_e{j} = GradientOfPotential(phi(indx), phi_n_e(indx), S);

            %Maxwell stress . normal (traction)
            traction = maxwellSnor(H_e{j}, S, phisquared(indx)); 
            traction = real(reshape(traction.to_array,[],3))'; 
            fM(indv) = traction(:); 
        end

        fprintf('\n Forces and torques on amphiphillic Janus particles \n'); 
        FT = real(Ck*fM); 
        display(reshape(FT,6,n3))
    otherwise
        error('Type of problem passed in is not implemented.');
    end
end

function [sigma,mu,U,VW] = LOCAL_compute_velocities(sigma,VW,Ct,Kernels,Nullsp,Fparams,col,collist,i,dt, Mt, closest_points_1, closest_points_2)
    %{
    This does a mobility solve, and then analyzes whether any corrections are needed. If so,
    it does the correction using the collision resolution algorithm and returns the adjusted
    densities and velocities.

    Inputs

    Outputs
        sigma : density
        mu : 
        VW : rigid body velocities (translational/rotational velocities)
    %}
    global timings; 
    parslv = Fparams.parslv;
    tau = Fparams.parbd.tau;
    W = Fparams.parbd.W;
    num_body = Fparams.parbd.n3;
    rd = ones(num_body,1);
    rdt = repmat((rd.').^(-4),3,1);
    rdw = repmat((rd.').^(-2),3,1);
    vind = reshape(repmat(6*(0:num_body-1),3,1),1,[])+repmat((1:3),1,num_body);
    wind = reshape(repmat(6*(0:num_body-1),3,1),1,[])+repmat((4:6),1,num_body);
    ambient = SpheroidalMS_EvalBackgroundFlow(Fparams.parbd, Fparams.background_flow);
    parslv_td = LOCAL_build_td_main_solver_params(Kernels, Nullsp, Fparams, Ct, Mt, col, collist, ...
        closest_points_1, closest_points_2, i);

    % Do a generic fluid solve, and then adjust if collision occurs.
    if isempty(VW) || Fparams.comp
        % Fluid Solve
        % (1) U_inc=S[sigma] (particular solution given forces and torques)
        % sigma is the incoming traction distribution 
        
        % (2) U_sc=S[mu] ("scattered" field with zero forces and torques)
        % RHS -(aI+K)*sigma
        tic;
        % Note that Kernels.TD = 0.5I + K + L. So, we need to the
        % L[\sigma] term below to get rid of L.
        B = -ambient.traction + Nullsp.L*sigma - Lapp(Kernels.TD,sigma); 
        timings.velocities.apply(i) = 0.5*toc;
        
        % Solve Fredholm eq (aI + K + L)*mu = -(aI+K)*sigma for mu
        tic; 
        mu = Lslv(Kernels.TD,B,parslv_td);
        fprintf('\n Time for solve: %e',toc); 
        timings.velocities.solve(i) = toc;  

        tic; 
        U = LOCAL_rebuild_total_surface_velocity(Kernels, sigma, mu, ambient); 
        fprintf('\n Time for apply (of S) to compute U: %e',toc);
        timings.velocities.apply(i) = timings.velocities.apply(i) + 0.5*toc;

        tic;
        VW = LOCAL_extract_RBM(U,Nullsp,W,tau,rdw,rdt,vind,wind,num_body);
        fprintf('\n Time to compute V and W: %e',toc); 
        timings.velocities.vw(i) = toc; 

        errbs = norm(U-Nullsp.D'*VW(:))/norm(U);
        fprintf('\n Rigid body velocity error: %e',errbs)
    elseif col
        % If comp=false and we entered with collision info, still compute baseline solve.
        tic;
        B = -ambient.traction + Nullsp.L*sigma - Lapp(Kernels.TD,sigma);
        timings.velocities.apply(i) = toc;

        tic;
        mu = Lslv(Kernels.TD,B,parslv_td);
        timings.velocities.solve(i) = toc;

        U = LOCAL_rebuild_total_surface_velocity(Kernels, sigma, mu, ambient);
        VW = LOCAL_extract_RBM(U,Nullsp,W,tau,rdw,rdt,vind,wind,num_body);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Collision event
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if col
        tic;

        fprintf('\n-------------------------------------------------');
        fprintf('\n Computing contact force at collision sites: \n')
        display(collist(:,1:2)')
        fprintf('-------------------------------------------------\n');

        [F_c,mu_c,rho_c] = ...
            LOCAL_Compute_Contact_LCP(collist,Kernels,Nullsp,Fparams,Ct,VW,dt,Mt,closest_points_1,closest_points_2,i);

        F_c = reshape(F_c,6,[]);
        display(F_c(1:3,:));

        if ~isempty(mu_c)
            mu = mu + mu_c; 
            sigma = sigma + rho_c; 
            U = LOCAL_rebuild_total_surface_velocity(Kernels, sigma, mu, ambient); 
            VW = LOCAL_extract_RBM(U,Nullsp,W,tau,rdw,rdt,vind,wind,num_body);
        end

        fprintf('\n Time for contact force correction: %e',toc);
        timings.velocities.col(i) = toc;

        errbs = norm(U-Nullsp.D'*VW(:))/norm(U);
        fprintf('\n Rigid body velocity error after collision correction: %e',errbs)
    end

    VW = real(VW);
end

function U = LOCAL_rebuild_total_surface_velocity(Kernels, sigma, mu, ambient)
    U = ambient.velocity + Lapp(Kernels.SD, (mu + sigma));
end

function VW = LOCAL_extract_RBM(U,Nullsp,W,tau,rdw,rdt,vind,wind,num_body)
    CU = Nullsp.C*U;
    IU = CU(vind);
    WxI = CU(wind);

    VW = zeros(6,num_body);
    if ~iscell(W)
        VW(1:3,:) = (1/sum(W))*(rdw.*reshape(IU,3,num_body));
        VW(4:6,:) = tau\(rdt.*reshape(WxI,3,num_body));
    else
        IUv = reshape(IU,3,num_body);
        WxIv = reshape(WxI,3,num_body);
        for j=1:num_body
            VW(1:3,j) = (1/sum(W{j}))*IUv(:,j);
            VW(4:6,j) = tau{j}\WxIv(:,j);
        end
    end
end
    
function Ctp = LOCAL_advance_center(Ct,dt,VW)
    Ctp = Ct + dt*VW(1:3,:)';
end

function [Mtp,Xtp,normW] = LOCAL_advance_rotation(MRot, VW, Mt, Xt, X0, dt, np, num_body)
    normW = zeros(num_body,1); Mtp=Mt; Xtp=Xt;
    for k=1:num_body
        xind = (1:np)+np*(k-1); 
        normW(k) = norm(VW(4:6,k));

        if normW(k) > 1e-14 % Checks if significant enough to update.
            Mtp{k} = MRot(VW(4:6,k),dt)*Mt{k};          
            Xtp(xind,:) = X0(xind,:)*Mtp{k}';
        else   
            Mtp{k} = Mt{k}; 
            Xtp(xind,:) = Xt(xind,:);    
        end   
    end
end

%% Collision resolution
function [F_c,mu_c,rho_c,lcp_info] = ...
    LOCAL_Compute_Contact_LCP(collist,Kernels,Nullsp,Fparams,Ct,VW,dt,Mt,closest_points_1,closest_points_2,timestep_idx)
    %{
    Calculates the force and the resulting densities due to the collision.
    Note that this does not resolve the collision, only calculates the effect of it.

    Inputs:
    collist:
    closest_points_1:
        double 3 x number of collision pairs of points on body 1 of each colliding pair
    closest_points_2:
        double 3 x number of collision pairs of points on body 2 of each colliding pair
    Kernels
    Nullsp - 
    Fparams - parameters of mobility solver
    Ct - centers of bodies
    Mt       - (1-by-n3 cell)
        array of 3x3 rotation matrices at time t
    Ct       - (double n3 x 3) (will transpose this for column vectors)
    % centers are row vectors
    Ct = Ct';
    VW - translational/rotational velocities
    dt - timestep size
    Mt - rotation matrices of bodies
        sigma    - (double 3*np*n3 x 1)
        incident traction density (particular)
    mu       - (double 3*np*n3 x 1)
        scattered density from Fredholm solve
    VW       - (double 6 x n3
        velocities are column vectors
        rigid-body velocities; 1:3 translational V, 4:6 angular W per body
    closest_points_1 - closest points on body 1 of each colliding pair
    closest_points_2 - closest points on body 2 of each colliding pair

    Outputs:
    F_c -
    mu_c
    rho_c
    %}
    parslv = Fparams.parslv;
    n3 = Fparams.parbd.n3;
    collision_eps = Fparams.parbd.collision_eps;
    max_radii = max(Fparams.parbd.equ_radii, Fparams.parbd.polar_radii);

    body_1_idx = collist(:,1);
    body_2_idx = collist(:,2);
    num_contact_pairs = length(body_1_idx);
    if num_contact_pairs == 0
        mu_c = [];
        rho_c = [];
        F_c = [];
        lcp_info = struct('err', 0, 'iter', 0, 'flag', 0, 'msg', 'no-contact');
        return;
    end

    [distances, closest_points_1, closest_points_2] = LOCAL_spheroidal_distances(Ct, Mt, Fparams, collist);
    distances = distances(:); % Force into column vector

    TD = Kernels.TD;
    SD = Kernels.SD;
    Bk = Nullsp.B;
    Lk = Nullsp.L;
    Ak = Nullsp.A;

    F = LOCAL_build_contact_force_matrix(collist, closest_points_1, closest_points_2, Ct, Mt, Fparams.parbd);

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Build A = F^T * M * F %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    dense_contact_max_pairs = max(0, round(LOCAL_get_solver_option(parslv, 'col_dense_max_pairs', 4)));
    matfree = num_contact_pairs > dense_contact_max_pairs;

    % Separate outer solves from contact-LCP inner solves
    td_build_info = LOCAL_make_td_build_info(Fparams, Nullsp, timestep_idx);

    contact_parslv = parslv;
    % Contact-LCP inner TD solves seemed to be more robust with plain restarted GMRES
    % than with recycled variants.
    contact_parslv.solver = 'gmres';

    contact_parslv = LOCAL_configure_td_contact_preconditioner(contact_parslv, TD, ...
        size(Bk, 2), n3, collist, closest_points_1, closest_points_2, distances, 'Contact TD', ...
        td_build_info);
    if ~matfree
        fprintf('\n Contact LCP projected A: dense build for %d pairs ', num_contact_pairs);
        contact_parslv.tol = parslv.coltol;
        Amat = LOCAL_build_projected_contact_matrix(TD, SD, Lk, Bk.', Ak, F, contact_parslv);
    else
        contact_parslv.tol = parslv.coltol;
        rho_c = @(x) (Bk.')*(F*x);
        % Contact-LCP inner TD solves forced to plain GMRES
        Amat = @(x) real(F.'*(Ak*Lapp(SD,Lslv(TD,-Lapp(TD,rho_c(x))+Lk*rho_c(x),contact_parslv)+rho_c(x))));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Build b = (1/dt)*phi + F^T * V %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Recall that phib is the minimum separation function
    pair_scale = max(max_radii(body_1_idx), max_radii(body_2_idx));
    % We subtract by collision_eps*pair_scale here to account for the buffer between spheroids:
    % that is, if phi > 0, we are outside the buffer, and if phi < 0, we are inside the buffer, so we
    % are too close/overlap! This is needed since in the theoretical model \Phi assumes we can measure
    % distances and account for collisions exactly.
    phib = (1/dt) * (distances - collision_eps*pair_scale);
    bvec = phib + real((F.')*VW(:));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %LCP solve
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    max_iter=parslv.colmaxit; 
    tol_rel=parslv.col_tolrel;
    tol_abs=parslv.col_tolabs;
    profile=1;
    lam0 = zeros(size(bvec));

    lcp_flag = NaN;
    lcp_msg = '';
    if matfree
        switch parslv.colsolver
            case 'Newton'
                [lam, err, iter, lcp_flag, ~, lcp_msg] = ...
                    minmap_newton_matfree(Amat, bvec, lam0, max_iter, tol_rel, tol_abs, profile);
            case 'APGD'
                [lam, err, iter, lcp_flag, ~, lcp_msg] = ...
                    APGD_matfree(Amat, bvec, lam0, max_iter, tol_rel, tol_abs, profile);
            case 'BBPGD'
                [lam, err, iter, lcp_flag, ~, lcp_msg] = ...
                    BBPGD_matfree(Amat, bvec, lam0, max_iter, tol_rel, tol_abs, profile);
            otherwise
                error('Invalid LCP solver in params.');
        end
    else
        switch parslv.colsolver
            case 'Newton'
                [lam, err, iter, lcp_flag, ~, lcp_msg] = ...
                    minmap_newton(Amat, bvec, lam0, max_iter, tol_rel, tol_abs, profile);
            case 'APGD'
                [lam, err, iter, lcp_flag, ~, lcp_msg] = ...
                    APGD(Amat, bvec, lam0, max_iter, tol_rel, tol_abs, profile);
            case 'BBPGD'
                [lam, err, iter, lcp_flag, ~, lcp_msg] = ...
                    BBPGD(Amat, bvec, lam0, max_iter, tol_rel, tol_abs, profile);
            otherwise
                error('Invalid LCP solver in params.');
        end
    end

    lcp_info = struct('err', err, 'iter', iter, 'flag', lcp_flag, 'msg', lcp_msg);
    fprintf(['\n' parslv.colsolver ' LCP solution error = %e, iters = %d, flag = %d (%s) \n'], ...
        err, iter, lcp_flag, lcp_msg);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contact forces and modified densities
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if norm(lam) > 0 % Is at least one constraint active?
        % lam is the contact force magnitudes, so need to convert back to vector form
        % using the contact-force directions.
        F_c = F*lam;
        
        % Given the new data (i.e. force/torque) on surface, correct the densities
        % by doing another mobility solve.
        [mu_c,rho_c] = LOCAL_mobility_solve(TD, SD, Lk, Bk.', [], F_c, ...
            LOCAL_post_contact_solver_params(parslv, TD, size(Bk, 2), n3, ...
                collist, closest_points_1, closest_points_2, distances, td_build_info) ...
        );
    else % Else, there are no collision-related corrections that are needed.
        mu_c = [];
        rho_c = [];
        F_c = [];
    end
end

function Amat = LOCAL_build_projected_contact_matrix(TD, SD, Lk, BkT, Ak, F, contact_parslv)
    % Build the small dense projected contact matrix A = F' * M * F one
    % contact column at a time. This avoids forming the full mobility
    % matrix while still giving an exact projected operator for small
    % contact sets.
    num_contact_pairs = size(F, 2);
    Amat = zeros(num_contact_pairs, num_contact_pairs);
    for j = 1:num_contact_pairs
        rho_col = BkT * F(:, j);
        BIE_RHS = -Lapp(TD, rho_col) + Lk * rho_col;
        mu_col = Lslv(TD, BIE_RHS, contact_parslv);
        mobility_col = Ak * Lapp(SD, mu_col + rho_col);
        Amat(:, j) = real(F.' * mobility_col);
    end
end

function [Ctp, Mtp, Xtp, normW, dt, colevent, collist, closest_points_1, closest_points_2] ...
    = LOCAL_collision_info(VW, Ct, Ctp, Mt, Mtp, MRot, Xt, Xtp, X0, dt, normW, Fparams)
    %{
    This function ONLY does the following:
        (1) validate the proposed next state after explicit advancement,
        (2) and if validation fails, recomputes the proposed state for the next timestep
            by continually halving the timestep until a suitable spatial state is found.
        (3) and upon this, returns a valid spatial state (assuming maxbis is never achieved).

    Inputs
    Fparams - struct of params
    X2 - I think this is the surface points of the bodies (not clear why this matters for spheroidal collision)
    Ct - Centers of bodies at current time
    Ctp - Centers of the bodies at next time (not clear which is current and which is next)
    VW - translational/rotational velocities
    MRot - This is a function handle for a rotation matrix, takes in a timestep, and angular velocity vector (why do we need a timestep?)
    Mt - Has the rotational information of the bodies at each timestep
    dt - timestep size

    Outputs
    colevent - boolean, did a collision event happen?
    collist - list of colliding pairs
    dt - possibly modified timestep size
    Ctp - possibly modified centers of bodies at next time

    %}
    n3 = Fparams.parbd.n3;
    np = Fparams.parbd.np;
    collision_eps = Fparams.parbd.collision_eps;
    maxbis = 10; % Conservative bisection guard (mostly to handle strong forces)

    % Note that the following pattern is a do-while loop.

    % Validate the proposed spatial state for the next timestep
    [colevent,collist,mindst,distances,closest_points_1,closest_points_2] = ...
        LOCAL_check_collision(Ctp, Mtp, Fparams);
    cond = colevent && (mindst < 0.1*collision_eps);

    bis = 0;
    % Recompute candidate state with halved dt, then test if the candidate is too close:
    % i.e. mindst is less than some tolerance (0.1*collision_eps).
    while cond && bis < maxbis
        dt = dt/2;
        bis = bis + 1;
        Ctp = LOCAL_advance_center(Ct,dt,VW);
        [Mtp,Xtp,normW] = ...
            LOCAL_advance_rotation(MRot, VW, Mt, Xt, X0, dt, np, n3);

        [colevent,collist,mindst,distances,closest_points_1,closest_points_2] = ...
            LOCAL_check_collision(Ctp, Mtp, Fparams);
        cond = colevent && (mindst < 0.1*collision_eps);

        fprintf('\n bisection = %d: dt = %1.4e, mindst = %1.4e', ...
            bis, dt, mindst);
    end

    if cond
        fprintf('\n Warning: maximum number of bisections (%d) reached.', maxbis);
    end

    if colevent
        fprintf('\n Min pairwise relative distance after collision detection: %2.4f ',mindst);
    else
        collist = [];
        fprintf('\n Min pairwise relative distance: %2.4f',mindst);
    end
end
    
function [colevent,collist,mindst] = LOCAL_check_collision_sph(C,Fparams)
    % Calculates distance between spheres using their centers and radius.
    n3 = size(C,1); 
    
    max_radii = max(Fparams.parbd.equ_radii, Fparams.parbd.polar_radii);
    diam = 2 * max_radii; 
    collision_eps = Fparams.parbd.collision_eps;
    scaling_factor = 1.1;
    
    if n3 > 1
        % Compute center distances
        distC = get_distances_between_centers(C);
    
        % Find pairs within a buffer band of the circumscribed spheres.
        [ii,jj]=meshgrid(1:n3); 
        %this creates all combinations of indices (i, j) for i,j=1,...,n3 (all bodies)
        %finds all the indices whose distance is less than diam (r_i + r_j) + buffer which is relative as its in terms of 1.1*collision_eps*max(r_i, r_j)
        id = distC<=diam+scaling_factor*collision_eps*max_radii & ii<jj; 
        %This selects the indices where the 
        ip = ii(id); 
        jp = jj(id); 
    
        % Compute minimum relative distance between spheres
        reldist = (distC-diam)./max_radii;
        mindst=min(reshape(reldist(~eye(n3)),[],1));
    else
        ip=[]; jp=[]; 
        mindst=Inf; 
    end

    colevent = mindst < scaling_factor*collision_eps; 
    collist = [ip jp]; 
end

function [colevent, collist, mindst, distances, closest_points_1, closest_points_2] = ...
        LOCAL_check_collision(C, Mt, Fparams)
    % C - centers of bodies (This is an double array, n_b x 3)
    % Mt - rotation matrices of bodies at current time

    % Do a quick check using the spheroids' circumscribed spheres.
    [~, collist_sph, ~] = LOCAL_check_collision_sph(C,Fparams);
    collist = collist_sph;

    if isempty(collist)
        colevent = false;
        mindst = Inf; % Minimum distance
        distances = [];
        closest_points_1 = [];
        closest_points_2 = [];
        return;
    end
    
    % we are going to implement a not vectorized version for now, just to try and get this to work (and we will assume this is only for spheroids)

    [distances, closest_points_1, closest_points_2] = ...
        LOCAL_spheroidal_distances(C, Mt, Fparams, collist);

    % Compute refined distances per collision candidate pair.
    collision_eps = Fparams.parbd.collision_eps;
    dist = distances(:);

    % Set the new collision state using the finer scale.
    mindst = min(dist);
    active = dist < 1.1*collision_eps;
    collist = collist(active,:);
    distances = distances(active).';
    closest_points_1 = closest_points_1(:,active);
    closest_points_2 = closest_points_2(:,active);
    colevent = ~isempty(collist);
end

%% Utility functions
function LOCAL_log_input_options(fname, Fparams, init)
    function value_str = LOCAL_parse_string_param(value, empty_label)
        if nargin < 2
            empty_label = '<empty>';
        end

        if isempty(value)
            value_str = empty_label;
        else
            value_str = strtrim(value);
        end
    end

    fprintf('\n============================================================\n');
    fprintf('spheroidal_mobility input options\n');
    fprintf('fname: %s\n', LOCAL_parse_string_param(fname, '<empty>'));
    fprintf('init: %s\n', LOCAL_parse_string_param(init, '<none>'));
    fprintf('Fparams:\n');
    if isempty(Fparams)
        fprintf('  <empty>\n');
    else
        fprintf('%s', evalc('disp(Fparams)'));
        if isfield(Fparams, 'parslv')
            fprintf('Fparams.parslv:\n');
            fprintf('%s', evalc('disp(Fparams.parslv)'));
        end
    end
    fprintf('============================================================\n\n');
end

function correction_parslv = LOCAL_post_contact_solver_params(parslv, TD, problem_size, n3, contact_pairs, closest_points_1, closest_points_2, pair_distances, build_info)
    correction_parslv = parslv;
    correction_parslv.solver = 'gmres';
    correction_parslv.deflate = false;
    if isfield(correction_parslv, 'deflate_basis')
        correction_parslv.deflate_basis = [];
    end
    correction_parslv = LOCAL_configure_td_contact_preconditioner(correction_parslv, TD, ...
        problem_size, n3, contact_pairs, closest_points_1, closest_points_2, pair_distances, ...
        'Post-contact TD', build_info);
end

function parslv_td = LOCAL_build_td_main_solver_params(Kernels, Nullsp, Fparams, Ct, Mt, col, collist, closest_points_1, closest_points_2, timestep_idx)
    % Start from the supplied TD solver settings, then specialize them
    % for the current contact or near-contact geometry.
    parslv_td = Fparams.parslv;

    % Identify the contact pairs. This can be the active collision set or
    % a small near-contact probe set used to switch early.
    [contact_pairs, contact_points_1, contact_points_2, contact_distances] = ...
        LOCAL_select_td_main_contact_pairs(parslv_td, Ct, Mt, Fparams, col, collist, ...
            closest_points_1, closest_points_2);

    % Build and attach a coarse basis from the contact directions when deflation is enabled.
    parslv_td = LOCAL_attach_contact_deflation(parslv_td, Nullsp.B, Ct, Mt, Fparams, ...
        contact_pairs, contact_points_1, contact_points_2);

    % Attach the TD preconditioner that for the active/near-contact pairs.
    parslv_td = LOCAL_configure_td_contact_preconditioner(parslv_td, Kernels.TD, ...
        size(Nullsp.B, 2), Fparams.parbd.n3, contact_pairs, contact_points_1, contact_points_2, contact_distances, 'TD_main', ...
        LOCAL_make_td_build_info(Fparams, Nullsp, timestep_idx));
end

function build_info = LOCAL_make_td_build_info(Fparams, Nullsp, timestep_idx)
    build_info = struct( ...
        'parbd', Fparams.parbd, ...
        'Lk', Nullsp.L, ...
        'BkT', Nullsp.B.', ...
        'typeMV', Fparams.typeMV, ...
        'timestep_idx', timestep_idx ...
    );
end

function parslv = LOCAL_attach_contact_deflation(parslv, Bk, Ct, Mt, Fparams, collist, closest_points_1, closest_points_2)
    if ~LOCAL_get_solver_option(parslv, 'deflate', false) || isempty(collist)
        parslv.deflate_basis = [];
        return;
    end

    if isempty(closest_points_1) || isempty(closest_points_2)
        [~, closest_points_1, closest_points_2] = LOCAL_spheroidal_distances(Ct, Mt, Fparams, collist);
    end

    max_dim = LOCAL_get_solver_option(parslv, 'deflate_max_dim', min(8, size(collist, 1)));
    sval_tol = LOCAL_get_solver_option(parslv, 'deflate_sval_tol', 1e-2);
    % Always build the augmentation basis one connected contact component at
    % a time. The component builder uses a contact-frame mode matrix
    % (normal plus tangential relative-motion directions), not the physical
    % frictionless LCP force map.
    basis = LOCAL_build_component_augmented_contact_basis(Bk.', Ct, Mt, Fparams.parbd, ...
        collist, closest_points_1, closest_points_2, max_dim, sval_tol);
    parslv.deflate_basis = basis;
end

function [contact_pairs, closest_points_1, closest_points_2, pair_distances] = ...
        LOCAL_select_td_main_contact_pairs(parslv, Ct, Mt, Fparams, col, collist, input_points_1, input_points_2)
    contact_pairs = [];
    closest_points_1 = [];
    closest_points_2 = [];
    pair_distances = [];

    % Only necessary if the sparse NN preconditioner or deflation is enabled
    if ~(LOCAL_get_solver_option(parslv, 'deflate', false) || LOCAL_get_solver_option(parslv, 'td_sparse_nn', false))
        return;
    end

    collision_eps = Fparams.parbd.collision_eps;
    if col && ~isempty(collist)
        contact_pairs = collist;
        if isempty(input_points_1) || isempty(input_points_2)
            [distances, closest_points_1, closest_points_2] = LOCAL_spheroidal_distances(Ct, Mt, Fparams, collist);
        else
            closest_points_1 = input_points_1;
            closest_points_2 = input_points_2;
            distances = LOCAL_spheroidal_distances(Ct, Mt, Fparams, collist);
        end
        pair_distances = distances(:);
        return;
    end

    % If there is no active collision list yet, we treat a few very close pairs
    % as "contact-like".
    if ~LOCAL_get_solver_option(parslv, 'deflate_near_contact', true)
        return;
    end

    % Probe only a small candidate set, then keep those whose true gap is
    % within the threshold_ratio multiple of collision_eps.
    threshold_ratio = LOCAL_get_solver_option(parslv, 'deflate_near_contact_factor', 1.5);
    probe_pairs = LOCAL_select_probe_pairs(Ct, Fparams, ...
        LOCAL_get_solver_option(parslv, 'deflate_probe_pairs', 3));
    if isempty(probe_pairs)
        return;
    end

    % Refine the probe set with the actual spheroidal distance calculations.
    [probe_distances, probe_points_1, probe_points_2] = LOCAL_spheroidal_distances(Ct, Mt, Fparams, probe_pairs);
    near_rows = probe_distances(:) < threshold_ratio * collision_eps;
    if ~any(near_rows)
        return;
    end

    % Return only the probe pairs that are genuinely close enough.
    contact_pairs = probe_pairs(near_rows, :);
    closest_points_1 = probe_points_1(:, near_rows);
    closest_points_2 = probe_points_2(:, near_rows);
    pair_distances = probe_distances(near_rows);
end

function probe_pairs = LOCAL_select_probe_pairs(Ct, Fparams, max_pairs)
    % Choose a small set of body pairs worth checking with the full
    % spheroidal distance routine used by the near-contact code.
    probe_pairs = [];

    n3 = size(Ct, 1);
    if n3 < 2
        return;
    end

    [~, candidate_pairs, ~] = LOCAL_check_collision_sph(Ct, Fparams);
    if isempty(candidate_pairs)
        % Fall back to the globally nearest body centers.
        distC = get_distances_between_centers(Ct);
        distC(1:n3+1:end) = Inf;
        [ii, jj] = find(triu(true(n3), 1));
        center_dists = distC(sub2ind([n3 n3], ii, jj));
        [~, order] = sort(center_dists, 'ascend');
        keep = order(1:min(max_pairs, numel(order)));
        probe_pairs = [ii(keep) jj(keep)];
        return;
    end

    body_1_idx = candidate_pairs(:, 1);
    body_2_idx = candidate_pairs(:, 2);
    max_radii = max(Fparams.parbd.equ_radii, Fparams.parbd.polar_radii);
    distC = get_distances_between_centers(Ct);

    % First, do a pass through the distances through spherical distances
    diam = max_radii(body_1_idx) + max_radii(body_2_idx);
    rel_gap = (distC(sub2ind(size(distC), body_1_idx, body_2_idx)) - diam) ./ ...
        max(max_radii(body_1_idx), max_radii(body_2_idx));
    [~, order] = sort(rel_gap, 'ascend');
    keep = order(1:min(max_pairs, numel(order)));
    probe_pairs = candidate_pairs(keep, :);
end

function parslv = LOCAL_configure_td_contact_preconditioner(parslv, TD, problem_size, n3, contact_pairs, closest_points_1, closest_points_2, pair_distances, label, build_info)
    if nargin < 10
        build_info = [];
    end
    if isempty(contact_pairs)
        return;
    end

    if LOCAL_get_solver_option(parslv, 'td_sparse_nn', false)
        parslv = LOCAL_configure_td_sparse_nn_preconditioner(parslv, TD, problem_size, n3, contact_pairs, closest_points_1, closest_points_2, pair_distances, label, build_info);
        return;
    end
end

function parslv = LOCAL_configure_td_sparse_nn_preconditioner(parslv, TD, problem_size, n3, contact_pairs, closest_points_1, closest_points_2, pair_distances, label, build_info)
    % Decides which bodies belong in the sparse local block, asks the low-level
    % builder for the LU factorizations, and then installs the resulting apply function
    % into parslv.prec.
    if nargin < 10
        build_info = [];
    end

    sparse_inverse_method = LOCAL_get_solver_option(parslv, 'td_sparse_nn_inverse', 'ilu');
    fprintf('\n %s sparse NN prec: building inverse=%s for %d active pairs \n', ...
        label, sparse_inverse_method, size(contact_pairs, 1));

    build_tic = tic;

    % Converts contact pair list into a graph representation. Then,
    % send the graph to the local sparse operator builder, which assembles,
    % regularizes, and factors.
    [active_bodies, adjacency, edge_count] = LOCAL_build_td_sparse_nn_graph(contact_pairs, n3);
    [local_state, local_stats] = LOCAL_build_td_sparse_nn_local_state(TD, problem_size, n3, ...
        active_bodies, adjacency, parslv, build_info);

    % Build the final preconditioner state for LOCAL_apply_td_sparse_nn_preconditioner:
    % keep the cached LU factorizations for the active bodies, but also retain the
    % original block-diagonal preconditioner for every other body.
    % See LOCAL_apply_td_sparse_nn_preconditioner for more details.
    state = local_state;
    state.num_dofs = problem_size;
    state.base_prec = LOCAL_get_solver_option(parslv, 'prec', []);

    degrees = full(sum(adjacency, 2) - 1);
    stats = struct( ...
        'ok', true, ...
        'reason', '', ...
        'cache_status', local_stats.cache_status, ...
        'sparse_inverse_method', sparse_inverse_method, ...
        'num_active_bodies', numel(active_bodies), ...
        'num_edges', edge_count - numel(active_bodies), ...
        'avg_degree', mean(max(degrees, 0)), ...
        'num_block_entries', nnz(adjacency), ...
        'factor_nnz', local_stats.factor_nnz, ...
        'assembly_time', local_stats.assembly_time, ...
        'condest_time', local_stats.condest_time, ...
        'condest_A_local', local_stats.condest_A_local, ...
        'factor_time', local_stats.factor_time, ...
        'build_time', toc(build_tic) ...
    );

    % Hand off preconditioner to parslv and log all debug information about it.
    parslv.prec = @(X) LOCAL_apply_td_sparse_nn_preconditioner(state, X);
    fprintf(['\n %s sparse NN prec: inverse=%s cache=%s bodies=%d edges=%d ' ...
        'avg_degree=%1.2f nnz_blocks=%d local_w=1 ' ...
        'build_time=%e assembly_time=%e condest_time=%e condest=%e factor_time=%e factor_nnz=%d \n'], ...
        label, stats.sparse_inverse_method, stats.cache_status, ...
        stats.num_active_bodies, stats.num_edges, stats.avg_degree, stats.num_block_entries, stats.build_time, ...
        stats.assembly_time, stats.condest_time, stats.condest_A_local, ...
        stats.factor_time, stats.factor_nnz);
end

function [local_state, local_stats] = LOCAL_build_td_sparse_nn_local_state(TD, problem_size, n3, active_bodies, adjacency, parslv, build_info)
    % Assemble the active-body local operator, regularize it, estimate its
    % condition number, and does a LU factorization on it.
    local_state = [];
    sparse_inverse_method = LOCAL_get_solver_option(parslv, 'td_sparse_nn_inverse', 'ilu');
    local_stats = struct( ...
        'ok', false, ...
        'reason', '', ...
        'cache_status', 'new', ...
        'factor_nnz', 0, ...
        'assembly_time', 0, ...
        'condest_time', 0, ...
        'condest_A_local', NaN, ...
        'factor_time', 0 ...
    );

    Nb = problem_size / n3;
    active_dim = Nb * numel(active_bodies);
    reg_scale = LOCAL_get_solver_option(parslv, 'td_sparse_nn_reg', 1e-10);
    reuse_local_state = LOCAL_get_solver_option(parslv, 'td_sparse_nn_reuse_local_state', true);

    cache_key = [];
    if reuse_local_state
        cache_key = build_info.timestep_idx;
        cached_value = LOCAL_manage_td_sparse_nn_cache('get', cache_key, []);
        if ~isempty(cached_value)
            local_state = cached_value.local_state;
            local_stats = cached_value.stats;
            local_stats.cache_status = 'reused';
            local_stats.assembly_time = 0;
            local_stats.factor_time = 0;
            return;
        end
    end

    assembly_tic = tic;
    A_local = LOCAL_build_td_sparse_nn_pair_block_operator(TD, build_info, active_bodies, adjacency, Nb);
    assembly_time = toc(assembly_tic);

    diag_abs = abs(diag(A_local));
    scale = max(1, full(max(diag_abs)));
    
    if issparse(A_local)
        A_local = A_local + reg_scale * scale * speye(active_dim);
    else
        A_local = A_local + reg_scale * scale * eye(active_dim);
    end

    if LOCAL_get_solver_option(parslv, 'td_sparse_nn_estimate_cond', false)
        condest_tic = tic;
        local_stats.condest_A_local = condest(sparse(A_local));
        local_stats.condest_time = toc(condest_tic);
    end

    [inverse_state, inverse_stats] = LOCAL_build_td_sparse_nn_inverse(A_local, sparse_inverse_method, parslv);
    if ~inverse_stats.ok
        local_stats.reason = inverse_stats.reason;
        return;
    end

    local_state = struct( ...
        'active_dof_idx', LOCAL_get_body_idx(active_bodies, Nb), ...
        'A_local', A_local, ...
        'solver_type', inverse_state.solver_type, ...
        'P', inverse_state.P, ...
        'Q', inverse_state.Q, ...
        'R', inverse_state.R, ...
        'L', inverse_state.L, ...
        'U', inverse_state.U ...
    );
    local_stats.ok = true;
    local_stats.factor_nnz = inverse_stats.factor_nnz;
    local_stats.assembly_time = assembly_time;
    local_stats.factor_time = inverse_stats.factor_time;

    if reuse_local_state
        LOCAL_manage_td_sparse_nn_cache('set', cache_key, struct('local_state', local_state, 'stats', local_stats));
    end
end

function [active_bodies, adjacency, edge_count] = LOCAL_build_td_sparse_nn_graph(contact_pairs, n3)
    active_bodies = unique(contact_pairs(:));
    m = numel(active_bodies);
    adjacency = false(m, m);
    edge_count = 0;
    if m == 0
        return;
    end

    local_body_idx = zeros(n3, 1);
    local_body_idx(active_bodies) = 1:m;
    adjacency(1:(m + 1):m^2) = true;

    for pair_idx = 1:size(contact_pairs, 1)
        i_local = local_body_idx(contact_pairs(pair_idx, 1));
        j_local = local_body_idx(contact_pairs(pair_idx, 2));
        if i_local == 0 || j_local == 0 || i_local == j_local
            continue;
        end
        adjacency(i_local, j_local) = true;
        adjacency(j_local, i_local) = true;
    end

    adjacency = adjacency | adjacency.';
    edge_count = nnz(adjacency);
end

function A_local = LOCAL_build_td_sparse_nn_direct_operator(build_info, active_bodies, ~, ~)
    local_parbd = LOCAL_build_td_sparse_nn_local_parbd(build_info.parbd, active_bodies);
    local_dof_idx = LOCAL_get_body_idx(active_bodies, build_info.parbd.Nb);
    L_local = build_info.Lk(local_dof_idx, local_dof_idx);
    TD_local = SpheroidalMS_MatVec([], L_local, build_info.typeMV, local_parbd, ...
        local_parbd.kerd, 0.5, 'TSL_Stk_3D');
    A_local = real(TD_local);
end

function local_parbd = LOCAL_build_td_sparse_nn_local_parbd(parbd, active_bodies)
    % Only consider the active body geometry
    local_parbd = SpheroidalMS_set_params( ...
        equ_radii = parbd.equ_radii(active_bodies), ...
        polar_radii = parbd.polar_radii(active_bodies), ...
        p = parbd.p, ...
        C = parbd.C(active_bodies, :), ...
        collision_eps = parbd.collision_eps, ...
        mdist = parbd.mdist, ...
        doAna = parbd.doAna, ...
        flag_pot = "TSL_Stk_3D", ...
        kerd = parbd.kerd, ...
        dense = true, ...
        bodydist = parbd.bodydist, ...
        MRot = {parbd.MRot{active_bodies}}, ...
        tsl_dealiasing = parbd.tsl_dealiasing, ...
        tsl_dealiasing_pad = parbd.tsl_dealiasing_pad ...
    );
end

function A_local = LOCAL_build_td_sparse_nn_pair_block_operator(TD, build_info, active_bodies, adjacency, Nb)
    % Actually construct the preconditioner. If the full dense TD operator is
    % already available, the active-body block can be extracted directly
    % instead of being reassembled pair-by-pair.
    if isnumeric(TD) && ~issparse(TD)
        active_dof_idx = LOCAL_get_body_idx(active_bodies, Nb);
        A_local = sparse(real(TD(active_dof_idx, active_dof_idx)));
        return;
    end

    m = numel(active_bodies);
    active_dim = Nb * m;
    A_local = spalloc(active_dim, active_dim, max(1, nnz(adjacency)) * Nb * Nb);

    % Build diagonal self blocks through the dense self-eval path.
    for body_local = 1:m
        local_idx = (1:Nb) + Nb * (body_local - 1);
        A_self = LOCAL_build_td_self_block_prec(build_info.parbd, build_info.Lk, Nb, active_bodies(body_local));
        A_local(local_idx, local_idx) = sparse(real(A_self));
    end

    % Build each undirected off-diagonal body pair, then place them in the matrix.
    for src_local = 1:m
        src_idx = (1:Nb) + Nb * (src_local - 1);
        for dst_local = (src_local + 1):m
            if ~adjacency(dst_local, src_local)
                continue;
            end

            pair_dense = LOCAL_build_td_sparse_nn_direct_operator( ...
                build_info, [active_bodies(dst_local); active_bodies(src_local)], [], []);
            dst_idx = (1:Nb) + Nb * (dst_local - 1);

            A_local(dst_idx, src_idx) = sparse(real(pair_dense(1:Nb, (Nb + 1):(2 * Nb))));
            A_local(src_idx, dst_idx) = sparse(real(pair_dense((Nb + 1):(2 * Nb), 1:Nb)));
        end
    end
end

function [inverse_state, inverse_stats] = LOCAL_build_td_sparse_nn_inverse(A_local, sparse_inverse_method, parslv)
    inverse_state = struct('solver_type', '', 'P', [], 'Q', [], 'R', [], 'L', [], 'U', []);
    inverse_stats = struct( ...
        'ok', false, ...
        'reason', '', ...
        'factor_time', 0, ...
        'factor_nnz', 0 ...
    );

    factor_tic = tic;
    if strcmp(sparse_inverse_method, 'lu')
        [L, U, P, Q, R] = lu(sparse(A_local));
        inverse_state.solver_type = 'sparse-lu';
        inverse_state.P = P;
        inverse_state.Q = Q;
        inverse_state.R = R;
    else % ILU
        ilu_setup = struct( ...
            'type', 'crout', ...
            'droptol', LOCAL_get_solver_option(parslv, 'td_sparse_nn_droptol', 1e-3), ...
            'udiag', 1 ...
        );
        [L, U] = ilu(sparse(A_local), ilu_setup);
        inverse_state.solver_type = 'ilu';
    end
    inverse_state.L = L;
    inverse_state.U = U;
    inverse_stats.ok = true;
    inverse_stats.factor_time = toc(factor_tic);
    inverse_stats.factor_nnz = nnz(L) + nnz(U);
end

function Y = LOCAL_apply_td_sparse_nn_preconditioner(state, X)
    Y = LOCAL_apply_preconditioner(state.base_prec, X); % Apply block-diagonal preconditioner first
    if isempty(state.active_dof_idx) % No active bodies
        return;
    end

    % Only apply sparse NN preconditioner on active bodies (this also includes the block-diagonal)
    rhs = X(state.active_dof_idx, :);
    local_action = LOCAL_solve_td_sparse_nn_local_system(state, rhs);
    Y(state.active_dof_idx, :) = local_action;
end

function local_correction = LOCAL_solve_td_sparse_nn_local_system(state, local_residual)
    % Solve the active-body block using the cached LU factorization.
    % This replaces the base preconditioner action on the active DoFs with a
    % more accurate preconditioner for A_local * x = local_residual.
    switch state.solver_type
        case 'dense-lu'
            % Dense LU path: state.P is the row-pivot vector, so first
            % permute the residual into pivoted row order, then apply L/U (or U^{-1}L^{-1}).
            local_correction = state.U \ (state.L \ local_residual(state.P, :));
        case 'sparse-lu'
            % MATLAB sparse LU returns factors satisfying
            % P * (R \ A_local) * Q = L * U. Undo that ordering to recover the
            % solution of A_local * x = local_residual on this local block.
            local_correction = state.Q * (state.U \ (state.L \ (state.P * (state.R \ local_residual))));
        otherwise
            % ILU stores only the incomplete triangular factors, so this acts
            % as an approximate local inverse via forward/back substitution.
            local_correction = state.U \ (state.L \ local_residual);
    end
end

function dof_idx = LOCAL_get_body_idx(body_ids, Nb)
    % If each body owns Nb consecutive rows, then body k contributes
    % rows ((k-1)*Nb + 1) : (k*Nb).
    body_ids = body_ids(:);
    dof_idx = zeros(Nb * numel(body_ids), 1);
    block_start_idx = 1;
    for body_idx = 1:numel(body_ids)
        % Global row block for this body in the full TD system.
        rows = (1:Nb) + Nb*(body_ids(body_idx) - 1);
        % Store those rows contiguously so callers can index all active-body
        % DoFs with one vector, preserving the body order in body_ids.
        dof_idx(block_start_idx:(block_start_idx + Nb - 1)) = rows(:);
        block_start_idx = block_start_idx + Nb;
    end
end

function basis = LOCAL_build_augmented_contact_basis(BkT, F, max_dim, sval_tol)
    % Build a coarse basis from the dominant singular directions of a
    % contact-derived mode matrix F, then convert those directions into the
    % density space through an application of B^T.
    basis = [];
    if isempty(F) || max_dim <= 0
        return;
    end

    [Ub, S, ~] = svd(F, 'econ');
    sing = diag(S);

    % Keep only the dominant singular directions, then orthonormalize after
    % mapping them into the density space (i.e. after left-applying B^T).
    eligible_count = find(sing >= sval_tol * sing(1), 1, 'last');
    if isempty(eligible_count)
        eligible_count = 1;
    end
    keep = 1:min(max_dim, eligible_count);
    basis = orth(real(BkT * Ub(:, keep)));
end

function basis = LOCAL_build_component_augmented_contact_basis(BkT, Ct, Mt, parbd, collist, closest_points_1, closest_points_2, max_dim, sval_tol)
    % Build the deflation basis one connected contact component at a time,
    % then merge those component-local bases into one global coarse space.
    basis = [];
    if isempty(collist) || max_dim <= 0
        return;
    end

    components = LOCAL_contact_components(parbd.n3, collist);
    component_basis = zeros(size(BkT, 1), 0);
    for comp_idx = 1:numel(components)
        rows = components{comp_idx};
        % Build a solver-only contact-mode matrix for this connected
        % component. In addition to the normal contact direction, include
        % the two tangent directions to hopefully mediate sliding problems.
        F_comp = LOCAL_build_augmented_contact_mode_matrix(collist(rows, :), closest_points_1(:, rows), ...
            closest_points_2(:, rows), Ct, Mt, parbd);
        basis_comp = LOCAL_build_augmented_contact_basis(BkT, F_comp, size(F_comp, 2), sval_tol);
        component_basis = [component_basis, basis_comp];
    end

    % Re-orthonormalize after concatenation, then truncate to the requested
    % global dimension cap.
    basis = orth(real(component_basis));
    if size(basis, 2) > max_dim
        basis = basis(:, 1:max_dim);
    end
end

function components = LOCAL_contact_components(n3, collist)
    % Partition the contact-pair list into connected components. Two rows
    % belong to the same component if they are linked through shared bodies.
    components = {};
    if isempty(collist)
        return;
    end

    unassigned = true(size(collist, 1), 1);
    while any(unassigned)
        % Start from one unassigned contact row and grow the component until
        % no new rows touch any currently active body. Note that this is the
        % flood-fill algorithm.
        seed = find(unassigned, 1);
        active_rows = false(size(unassigned));
        active_bodies = false(n3, 1);
        active_bodies(collist(seed, :)) = true;

        changed = true;
        while changed
            % Pull in every remaining row that touches the current body set.
            hits = unassigned & (active_bodies(collist(:, 1)) | active_bodies(collist(:, 2)));
            changed = any(hits);
            if ~changed
                break;
            end
            active_rows = active_rows | hits;
            unassigned(hits) = false;
            bodies = unique(collist(active_rows, :));
            active_bodies(bodies) = true;
        end

        components{end+1} = find(active_rows);
    end
end

function [Mf,rho_c] = LOCAL_mobility_solve(TD, SD, Lk, BMR, VNS, f, parslv)
    % The VNS parameter remains baffling to me, so we shall keep it for now...
    % Apply contact-force map and recover correction densities.
    rho_c = BMR*f;
    BIE_RHS = -Lapp(TD,rho_c) + Lk*rho_c;
    mu_c = Lslv(TD,BIE_RHS,parslv);

    if ~isempty(VNS)
        Mf = VNS*Lapp(SD,mu_c+rho_c);
    else
        Mf = mu_c;
    end
end

function y = Lapp(A,x)
    % Left-apply the matrix A to the vector x.
    if isnumeric(A)
        y=A*x;
    else
        y=real(A(x));
    end
end

function x = Lslv(A,b,parslv)
    % Linear solve for Ax = b, with parameters given in parslv.
    prec = LOCAL_get_solver_option(parslv, 'prec', []);
    silent = logical(LOCAL_get_solver_option(parslv, 'silent', false));

    if isempty(prec)
        pr = [];
    else
        pr = prec;
    end

    % Dense/direct path
    if isnumeric(A)
        x = A\b;
        return;
    end

    A_op = LOCAL_get_matvec_func(A);
    x = zeros(size(b));
    use_deflation = isfield(parslv, 'deflate_basis') && ~isempty(parslv.deflate_basis);
    for i=1:size(b,2)
        if use_deflation
            [x(:,i),flag,rs,total_iters,it,deflation_info] = ...
                LOCAL_run_augmented_solver(A,b(:,i),pr,parslv);
            if ~silent
                fprintf('\n gmres %d flag=%d outer=%d inner=%d total=%d basis=%d coarse=%1.4g relres=%1.4g \n', ...
                    i, flag, it(1), it(2), total_iters, deflation_info.basis_dim, ...
                    deflation_info.coarse_relres, rs);
            end
        else
            restart = parslv.rst;
            tol = parslv.tol;
            maxit = parslv.maxit;
            [x(:,i),flag,rs,it] = gmres(A_op,b(:,i),restart,tol,maxit,pr);
            total_iters = LOCAL_total_gmres_iters(it, restart);
            if ~silent
                fprintf('\n gmres %d flag=%d outer=%d inner=%d total=%d relres=%1.4g \n', ...
                    i, flag, it(1), it(2), total_iters, rs);
            end
        end
    end
end

function value = LOCAL_get_solver_option(parslv, field_names, default_value)
    if ischar(field_names) || (isstring(field_names) && isscalar(field_names))
        field_names = {char(field_names)};
    end

    for j = 1:numel(field_names)
        field_name = field_names{j};
        if isfield(parslv, field_name) && ~isempty(parslv.(field_name))
            value = parslv.(field_name);
            return;
        end
    end

    value = default_value;
end

function [x, flag, relres, total_iters, iter, info] = LOCAL_run_augmented_solver(A, b, prec, parslv)
    A_op = LOCAL_get_matvec_func(A);

    solve_opts = struct( ...
        'rst', parslv.rst, ...
        'tol', parslv.tol, ...
        'maxit', parslv.maxit, ...
        'prec', prec, ...
        'deflate_basis', parslv.deflate_basis, ...
        'krylov_solver', 'gmres' ...
    );

    [x, flag, relres, iter, info] = SpheroidalMS_augmented_GMRES(A_op, b, solve_opts);
    total_iters = LOCAL_total_gmres_iters(iter, parslv.rst);
end

function value = LOCAL_manage_td_sparse_nn_cache(action, key, value_in)
    persistent keys values

    if isempty(keys)
        keys = [];
        values = {};
    end

    if nargin < 3
        value_in = [];
    end

    switch action
        case 'reset'
            keys = [];
            values = {};
            value = [];
        case 'get'
            idx = find(keys == key, 1);
            if isempty(idx)
                value = [];
            else
                value = values{idx};
            end
        case 'set'
            idx = find(keys == key, 1);
            if isempty(idx)
                keys(end+1) = key;
                values{end+1} = value_in;
            else
                values{idx} = value_in;
            end
            value = [];
        case 'remove'
            idx = find(keys == key, 1);
            if ~isempty(idx)
                keys(idx) = [];
                values(idx) = [];
            end
            value = [];
        otherwise
            error('Invalid sparse-NN cache action "%s".', action);
    end
end

function A_op = LOCAL_get_matvec_func(A)
    if isnumeric(A)
        A_op = @(x) A*x;
    else
        A_op = @(x) real(A(x));
    end
end

function F = LOCAL_build_contact_force_matrix(collist, closest_points_1, closest_points_2, Ct, Mt, parbd)
    n3 = parbd.n3;
    num_contact_pairs = size(collist, 1);
    F = zeros(6 * n3, num_contact_pairs);
    if num_contact_pairs == 0 || isempty(closest_points_1) || isempty(closest_points_2)
        return;
    end

    quadratic_forms = LOCAL_build_contact_quadratic_forms(Ct, Mt, parbd);

    body_1_idx = collist(:, 1);
    body_2_idx = collist(:, 2);
    for pair_idx = 1:num_contact_pairs
        translational_indi = (1:3) + 6 * (body_1_idx(pair_idx) - 1);
        translational_indj = (1:3) + 6 * (body_2_idx(pair_idx) - 1);
        rotational_indi = (4:6) + 6 * (body_1_idx(pair_idx) - 1);
        rotational_indj = (4:6) + 6 * (body_2_idx(pair_idx) - 1);

        cp_i = closest_points_1(:, pair_idx);
        cp_j = closest_points_2(:, pair_idx);

        body_j = body_2_idx(pair_idx);
        nvec = LOCAL_contact_normal_from_quadratic_form(quadratic_forms, body_j, cp_j, Ct);

        r_i = cp_i - Ct(body_1_idx(pair_idx), :).';
        r_j = cp_j - Ct(body_2_idx(pair_idx), :).';

        F(translational_indi, pair_idx) = nvec;
        F(translational_indj, pair_idx) = -nvec;
        F(rotational_indi, pair_idx) = cross(r_i, nvec);
        F(rotational_indj, pair_idx) = -cross(r_j, nvec);
    end
end

function F = LOCAL_build_augmented_contact_mode_matrix(collist, closest_points_1, closest_points_2, Ct, Mt, parbd)
    % Build solver-only contact-frame modes for the augmentation basis.
    % Each contact contributes:
    %   1. normal relative-motion mode
    %   2. first tangential relative-motion mode
    %   3. second tangential relative-motion mode
    %
    % The purpose of this is to help GMRES resolve sliding and near-contact modes.
    n3 = parbd.n3;
    num_contact_pairs = size(collist, 1);
    F = zeros(6 * n3, 3 * num_contact_pairs);
    if num_contact_pairs == 0 || isempty(closest_points_1) || isempty(closest_points_2)
        return;
    end

    quadratic_forms = LOCAL_build_contact_quadratic_forms(Ct, Mt, parbd);
    body_1_idx = collist(:, 1);
    body_2_idx = collist(:, 2);
    for pair_idx = 1:num_contact_pairs
        cp_i = closest_points_1(:, pair_idx);
        cp_j = closest_points_2(:, pair_idx);
        body_i = body_1_idx(pair_idx);
        body_j = body_2_idx(pair_idx);

        nvec = LOCAL_contact_normal_from_quadratic_form(quadratic_forms, body_j, cp_j, Ct);
        [t1, t2] = LOCAL_contact_tangent_frame(nvec);
        r_i = cp_i - Ct(body_i, :).';
        r_j = cp_j - Ct(body_j, :).';

        col_offset = 3 * (pair_idx - 1);
        F = LOCAL_store_contact_direction_mode(F, col_offset + 1, body_i, body_j, r_i, r_j, nvec);
        F = LOCAL_store_contact_direction_mode(F, col_offset + 2, body_i, body_j, r_i, r_j, t1);
        F = LOCAL_store_contact_direction_mode(F, col_offset + 3, body_i, body_j, r_i, r_j, t2);
    end
end

function quadratic_forms = LOCAL_build_contact_quadratic_forms(Ct, Mt, parbd)
    n3 = parbd.n3;
    quadratic_forms = zeros(3, 3, n3);
    for body_idx = 1:n3
        spheroid_params = LOCAL_get_spheroid_params(body_idx, Ct, Mt, parbd);
        D = diag([spheroid_params.a spheroid_params.b spheroid_params.c].^-2);
        quadratic_forms(:, :, body_idx) = spheroid_params.R * D * (spheroid_params.R.');
    end
end

function nvec = LOCAL_contact_normal_from_quadratic_form(quadratic_forms, body_idx, contact_point, Ct)
    grad = quadratic_forms(:, :, body_idx) * (contact_point - Ct(body_idx, :).');
    nvec = grad / norm(grad);
end

function [t1, t2] = LOCAL_contact_tangent_frame(nvec)
    nvec = nvec / norm(nvec);

    [~, axis_idx] = min(abs(nvec));
    reference_axis = zeros(3,1);
    reference_axis(axis_idx) = 1;

    t1 = cross(nvec, reference_axis);
    t1 = t1 / norm(t1);

    t2 = cross(nvec, t1);
    t2 = t2 / norm(t2);
end

function F = LOCAL_store_contact_direction_mode(F, col_idx, body_i, body_j, r_i, r_j, dvec)
    translational_indi = (1:3) + 6 * (body_i - 1);
    translational_indj = (1:3) + 6 * (body_j - 1);
    rotational_indi = (4:6) + 6 * (body_i - 1);
    rotational_indj = (4:6) + 6 * (body_j - 1);

    F(translational_indi, col_idx) = dvec;
    F(translational_indj, col_idx) = -dvec;
    F(rotational_indi, col_idx) = cross(r_i, dvec);
    F(rotational_indj, col_idx) = -cross(r_j, dvec);
end

function ITSSDd = LOCAL_build_td_inverse_blocks_prec(parbd, Lk, Nb, n3)
    Iblock = eye(Nb);
    ITSSDd = cell(n3,1);

    for k = 1:n3
        Akk = LOCAL_build_td_self_block_prec(parbd, Lk, Nb, k);
        ITSSDd{k} = Akk \ Iblock;
    end
end

function Akk = LOCAL_build_td_self_block_prec(parbd, Lk, ~, body_idx)
    Akk = SpheroidalMS_BuildTDSelfBlock(parbd, Lk, body_idx);
end

function Y = LOCAL_apply_preconditioner(prec, X)
    if isempty(prec)
        Y = X;
        return;
    end

    if isnumeric(prec)
        Y = prec\X;
    else
        Y = prec(X);
    end
end

function total_iters = LOCAL_total_gmres_iters(it, restart)
    if it(1) <= 0
        total_iters = it(2);
    else
        total_iters = (it(1)-1)*restart + it(2);
    end
end

function [distances, closest_points_1, closest_points_2] = LOCAL_spheroidal_distances(C, Mt, Fparams, collist)

    %switch statement for different distance algorithms.
    %right now we will just implement two options, moving balls and GJK signed volumes accelerated
    switch lower(Fparams.parbd.bodydist.algo)
        case 'moving_balls'
            distance_algo = @moving_balls_pair;
        case 'gjk signed volumes accelerated'
            distance_algo = @GJK_signed_volumes_accelerated_pair;
        otherwise
            error('Unknown spheroidal distance algorithm %s', Fparams.parbd.distance.algo);
    end

    %row vector of the pairwise distances
    distances = zeros(1, size(collist, 1));
    %row vector of the closest points of the 1st particle in the collision pair
    closest_points_1 = zeros(3, size(collist, 1));
    %row vector of the closest points of the 2nd particle in the collision pair
    closest_points_2 = zeros(3, size(collist, 1));
    % These are column vectors

    for collision_pair = 1:size(distances, 2)
        curr_index_1 = collist(collision_pair, 1);
        curr_index_2 = collist(collision_pair, 2);

        spheroid_1_params = LOCAL_get_spheroid_params(curr_index_1, C, Mt, Fparams.parbd);
        spheroid_2_params = LOCAL_get_spheroid_params(curr_index_2, C, Mt, Fparams.parbd);
        [closest_points_1(:, collision_pair), closest_points_2(:, collision_pair), distances(collision_pair)] = ...
            distance_algo(spheroid_1_params, spheroid_2_params, Fparams.parbd.bodydist.tol, Fparams.parbd.bodydist.max_iter);
    end

end

function spheroid_params = LOCAL_get_spheroid_params(body_idx, C, Mt, parbd)
    spheroid_params.C = C(body_idx, :);
    spheroid_params.R = Mt{body_idx};
    switch parbd.shape_type(body_idx)
        case "oblate"
            spheroid_params.a = parbd.polar_radii(body_idx);
            spheroid_params.c = parbd.equ_radii(body_idx);
        case "prolate"
            spheroid_params.a = parbd.equ_radii(body_idx);
            spheroid_params.c = parbd.polar_radii(body_idx);
        case "sphere"
            spheroid_params.a = parbd.equ_radii(body_idx);
            spheroid_params.c = parbd.equ_radii(body_idx);
        otherwise
            error('Unknown shape type.');
    end
    spheroid_params.b = spheroid_params.a;
end


function M = RotationMat(wh,t)
    %{
    An implementation of Rodrigues' rotation matrix formula.

    Inputs
    wh - (double) angular velocity vector
    t  - (double) timestep

    Outputs
    M - (double) 3x3 rotation matrix
    %}
    nwh = norm(wh); 
    t = nwh*t; 
    wh = wh./nwh; 
    
    M = [
        1-(wh(2)^2+wh(3)^2)*(1-cos(t)) , wh(2)*wh(1)*(1-cos(t))-wh(3)*sin(t) , wh(1)*wh(3)*(1-cos(t))+wh(2)*sin(t);...
        wh(1)*wh(2)*(1-cos(t))+wh(3)*sin(t),1-(wh(1)^2+wh(3)^2)*(1-cos(t)),wh(2)*wh(3)*(1-cos(t))-wh(1)*sin(t);...
        wh(1)*wh(3)*(1-cos(t))-wh(2)*sin(t),wh(2)*wh(3)*(1-cos(t))+wh(1)*sin(t),1-(wh(2)^2+wh(1)^2)*(1-cos(t))
    ];
end

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
        out         - (bool) external vs internal evaluation (set to 1) 
    
    parslv - (struct) linear solver parameters such as 
         prec     - (string) preconditioner type, '' for unprec, 'bkdiag'
                    (block diagonal), 'TT' (tensor train)
         precond_type - (string) 'bkdiag' or 'TT'
         solver   - gmres, pcg, bicg, etc. 
         tol (tolerance), maxit (maximum iterations), rst (restart), etc.
    
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
    ----

    %%%
    %%% CODE ANNOTATIONS
    %%%
    Mt stores rotational information for n bodies at each time step
    Mt0 is the initial rotation setup

    VW represents the translational (v)/angular velocity (w) each timestep (as a 6 x num_body matrix)
    VW0 represents the initial velocities

    Xrp are the tracking points needed for collision? Doesn't seem to be used either way.

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
        req = {'p','Ct','equ_radii','polar_radii','mdist','out'};
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
        assert(islogical(Fparams.parbd.out) || ismember(Fparams.parbd.out,[0,1]), ...
            'Fparams.parbd.out must be true or false.');
    end

    if isfield(Fparams, 'plotFlag')
        plot_enabled = Fparams.plotFlag;
    else
        plot_enabled = false;
    end

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
     
        for k=1:num_body
            Mt{1,k} = eye(3);   
        end
        
        Mt0 = Mt(1,:);
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
    Xrp = Fparams.parbd.Xrp; C0 = Fparams.parbd.C; 
    np = Fparams.parbd.np; normW = zeros(num_body,1);  
    if init_flag % Load state from memory
        for k=1:num_body
            normW(k) = norm(VW0(4:6, k));
            % Collect all discretization points for each body
            xind = (1:np) + np*(k-1);
            % Apply rotation (from loaded state) to each discretization point
            Xrp(xind,:) = Xrp(xind,:)*Mt0{k}';
        end 
    end
    
    % Get kernels/nullspace
    Kernels=[]; 
    [Kernels,Nullsp,Fparams,timings] = SpheroidalMS_UpdateOperators(Xrp,C0,Mt0,normW,Kernels,Fparams,timings,0); 
    timings.setup_kernel=toc;

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
                LOCAL_euler_step(Xt{i},Xt{1},X2,Mt(i,:),Ct{i},Kernels,Nullsp,Fparams,colevent,collist, closest_points_1, closest_points_2, t, dt,i);
        elseif strcmp(timedisc,'trapz')
            error('Calls need to be updated.')
            % (1) Predictor step: 
            fprintf('\n ---------------------------------------------------------- \n')
            fprintf('\n (1) Trapezoidal, predictor step \n')
            [Xt1,Mt1,Ct1,U1,FT1,sigma1,mu1,VW1,Kernels,Nullsp,...
                Fparams,colevent,collist,dt] = ...
                LOCAL_euler_step(Xt{i},Xt{1},X2,Mt(i,:),Ct{i},Kernels,Nullsp,Fparams,colevent,collist,t,dt,i);
            
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
                LOCAL_advance_step(VW{i},mu{i},sigma{i},Xt{i},Xt{1},X2,Mt(i,:),Ct{i},Kernels,Fparams,dt,i);
            
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
                LOCAL_euler_step(Xt{i},Xt{1},X2,Mt(i,:),Ct{i},Kernels,Nullsp,Fparams,colevent,collist,t41,dt41,i);
            
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
                LOCAL_advance_step(VW42,mu42,sigma42,Xt{i},Xt{1},X2,Mt(i,:),Ct{i},Kernels,Fparams,dt42,i);
            
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
                LOCAL_advance_step(VW43,mu43,sigma43,Xt{i},Xt{1},X2,Mt(i,:),Ct{i},Kernels,Fparams,dt43,i);
            
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
                LOCAL_advance_step(VW{i},mu{i},sigma{i},Xt{i},Xt{1},X2,Mt(i,:),Ct{i},Kernels,Fparams,dt,i);
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
    fprintf('\n Operator update')
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
    lprec=[]; 
    acc = parslv.tol; rst = parslv.rst; maxit = parslv.maxit;

    n3=Fparams.parbd.n3; p = Fparams.parbd.p; np = Fparams.parbd.np; 
    VW=[]; C = Fparams.parbd.C; 
    Bk = Nullsp.B; Ck = Nullsp.C; Lk = Nullsp.L; 
    
    switch Fparams.type
    case 'FTfun' % Force and torque are given
        Ffun = Fparams.Ffun; 
        Tfun = Fparams.Tfun; 
        Ct = Fparams.parbd.C; 
        
        Force = Ffun(t,Ct);  
        Torque = Tfun(t,Ct);
        FT = [Force;Torque]; FT = FT(:); 
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
        error('Need to implement kernels and make a new projection matrix.');
        SLMODD=Kernels.SLMODD; dSLMODD=Kernels.dSLMODD;
        DLMODD=Kernels.DLMODD; dDLMODD=Kernels.dDLMODD;
        flabel=Fparams.SurfaceLabel;
        % Solves for density, psi & uses them to compute normal derivative
        if isa(SLMODD,'function_handle')
            K = @(V) SLMODD(V) + DLMODD(V);
            P = make_projection(Fparams.parmod.p,1);
            Proj = @(V) reshape(P*reshape(V,[],n3),[],1);

            % 
            K_full_rank = @(V) Proj(K(Proj(V))) + V - Proj(V);    
            [psi,~,~,I]=gmres(K_full_rank, Proj(flabel),100,1e-6);
            fprintf('\n Janus Amph S+D BIE solve error = %e, %e', ...
                norm(K(psi)-Proj(flabel))/norm(flabel), norm(K_full_rank(psi)-Proj(flabel))/norm(flabel)); 
            
            phi = K(psi); 
            phi_n_e=dSLMODD(psi) + dDLMODD(psi);
        else
            K=SLMODD+DLMODD;
            P=make_projection(Fparams.parmod.p,n3);
            K_full_rank=P*K*P+(eye(np*n3)-P);
            [psi,~,~,I]=gmres(K_full_rank, P*(flabel),100,1e-6);
            %display(cond(K_full_rank));
            Kond = cond(K_full_rank); 

            fprintf('\n Janus Amph S+D BIE solve error = %e, %e',norm(K*psi-P*flabel)/norm(flabel), norm(K_full_rank*psi-P*flabel)/norm(flabel)); 
            phi = K*psi; 
            phi_n_e=dSLMODD*psi + dDLMODD*psi;
        end

        % Computes energy (not necessary for further dynamics)
        W = repmat(Fparams.parmod.W,n3, 1);
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
        B = Nullsp.L*sigma-Lapp(Kernels.TD,sigma); 
        timings.velocities.apply(i) = 0.5*toc;
        
        % Solve Fredholm eq (aI + K + L)*mu = -(aI+K)*sigma for mu
        tic; 
        mu = Lslv(Kernels.TD,B,parslv);
        fprintf('\n Time for solve: %e',toc); 
        timings.velocities.solve(i) = toc;  

        tic; 
        U = Lapp(Kernels.SD,(mu+sigma)); 
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
        B = Nullsp.L*sigma - Lapp(Kernels.TD,sigma);
        timings.velocities.apply(i) = toc;

        tic;
        mu = Lslv(Kernels.TD,B,parslv);
        timings.velocities.solve(i) = toc;

        U = Lapp(Kernels.SD,(mu+sigma));
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
            LOCAL_Compute_Contact_LCP(collist,Kernels,Nullsp,Fparams,Ct,VW,dt,Mt,closest_points_1,closest_points_2);

        F_c = reshape(F_c,6,[]);
        display(F_c(1:3,:));

        if ~isempty(mu_c)
            mu = mu + mu_c; 
            sigma = sigma + rho_c; 
            U = Lapp(Kernels.SD,(mu+sigma)); 
            VW = LOCAL_extract_RBM(U,Nullsp,W,tau,rdw,rdt,vind,wind,num_body);
        end

        fprintf('\n Time for contact force correction: %e',toc);
        timings.velocities.col(i) = toc;

        errbs = norm(U-Nullsp.D'*VW(:))/norm(U);
        fprintf('\n Rigid body velocity error after collision correction: %e',errbs)
    end

    VW = real(VW);
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
function [F_c,mu_c,rho_c] = ...
    LOCAL_Compute_Contact_LCP(collist,Kernels,Nullsp,Fparams,Ct,VW,dt,Mt,closest_points_1,closest_points_2)
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
    tol = parslv.tol;
    collision_eps = Fparams.parbd.collision_eps;
    max_radii = max(Fparams.parbd.equ_radii, Fparams.parbd.polar_radii);

    body_1_idx = collist(:,1);
    body_2_idx = collist(:,2);
    num_contact_pairs = length(body_1_idx);
    if num_contact_pairs == 0
        mu_c = [];
        rho_c = [];
        F_c = [];
        return;
    end

    % TODO: get rid of this...
    [distances, closest_points_1, closest_points_2] = LOCAL_spheroidal_distances(Ct, Mt, Fparams, collist);
    distances = distances(:); % Force into column vector

    TD = Kernels.TD;
    SD = Kernels.SD;
    Bk = Nullsp.B;
    Lk = Nullsp.L;
    Ak = Nullsp.A;

    % f(x) = (x-C)^T A (x-C) - 1 = 0, and \nabla f = 2*A*(x-C).
    % TODO: cache this.
    quadratic_forms = zeros(3,3,n3);
    for body_idx = 1:n3
        spheroid_params = LOCAL_get_spheroid_params(body_idx, Ct, Mt, Fparams.parbd);
        D = diag([spheroid_params.a spheroid_params.b spheroid_params.c].^-2);
        quadratic_forms(:,:,body_idx) = spheroid_params.R*D*(spheroid_params.R.');
    end

    % Build contact-direction matrix F
    F = zeros(6*n3,num_contact_pairs); 
    for k = 1:num_contact_pairs
        translational_indi = (1:3)+6*(body_1_idx(k)-1);
        translational_indj = (1:3)+6*(body_2_idx(k)-1);

        rotational_indi = (4:6)+6*(body_1_idx(k)-1);
        rotational_indj = (4:6)+6*(body_2_idx(k)-1);

        % Use the analytic spheroid gradient to build the contact normal.
        cp_i = closest_points_1(:,k);
        cp_j = closest_points_2(:,k);

        % body_i = body_1_idx(k);
        body_j = body_2_idx(k);

        % Convention: contact normal vector points from body 2 to body 1, so we 
        % use the outward normal of body 2 at cp_j.
        grad_j = quadratic_forms(:,:,body_j)*(cp_j - Ct(body_j,:).');
        nvec = grad_j/norm(grad_j);

        r_i = cp_i - Ct(body_1_idx(k),:).';
        r_j = cp_j - Ct(body_2_idx(k),:).';

        F(translational_indi,k) = nvec;
        F(translational_indj,k) = -nvec;
        F(rotational_indi,k) = cross(r_i, nvec);
        F(rotational_indj,k) = -cross(r_j, nvec);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Build A = F^T * M * F %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    matfree = ~Fparams.denseMV || num_contact_pairs > 1;

    if ~matfree
        rho_c = (Bk.')*F; % This is rho_c -- the incident field.
        BIE_RHS = -Lapp(TD,rho_c)+Lk*rho_c; % (-0.5I-T)*rho_c
        mu_c = Lslv(TD,BIE_RHS,parslv); % \mu_c = (0.5I + T + L)^{-1}(-0.5I - T)*rho_c

        % Setup LCP: x \perp A*x + b (dense build of Amat = F^T M F)
        % Lapp(SD, MuNS+Bf) represents the mobility solve, and conversion of density to velocity.
        % Ak*Lapp(SD,MuNS+Bf)) extracts the rigid body motion from the velocity (i.e. the mobility matrix).
        % So, Amat = F^T (A S (0.5I + T + L)^{-1}(-0.5 - T)B^T) F, where whatever is in the paranthesis is
        %   the mobility matrix \mathcal{M} in the equation \mathcal{U} = \mathcal{M}F.
        Amat = real(F.'*(Ak*Lapp(SD,mu_c+rho_c)));
    else % What is happening here...?
        parslv.tol = parslv.coltol;
        rho_c = @(x) (Bk.')*(F*x);
        Amat = @(x) real(F.'*(Ak*Lapp(SD,Lslv(TD,-Lapp(TD,rho_c(x))+Lk*rho_c(x),parslv)+rho_c(x))));
        parslv.tol = tol;
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
    tol_rel=parslv.col_tolrel; %1e-6; 
    tol_abs=parslv.col_tolabs; %1e-9; 
    profile=1;

    if matfree
        switch parslv.colsolver
            case 'Newton'
                [lam, err, iter, ~, ~, ~] = ...
                    minmap_newton_matfree(Amat, bvec, zeros(size(bvec)), max_iter, tol_rel, tol_abs, profile);
            case 'APGD'
                [lam, err, iter, ~, ~, ~] = ...
                    APGD_matfree(Amat, bvec, zeros(size(bvec)), max_iter, tol_rel, tol_abs, profile);
            case 'BBPGD'
                [lam, err, iter, ~, ~, ~] = ...
                    BBPGD_matfree(Amat, bvec, zeros(size(bvec)), max_iter, tol_rel, tol_abs, profile);
            otherwise
                error('Invalid LCP solver in params.');
        end
    else
        switch parslv.colsolver
            case 'Newton'
                [lam, err, iter, ~, ~, ~] = ...
                    minmap_newton(Amat, bvec, zeros(size(bvec)), max_iter, tol_rel, tol_abs, profile);
            case 'APGD'
                [lam, err, iter, ~, ~, ~] = ...
                    APGD(Amat, bvec, zeros(size(bvec)), max_iter, tol_rel, tol_abs, profile);
            case 'BBPGD'
                [lam, err, iter, ~, ~, ~] = ...
                    BBPGD(Amat, bvec, zeros(size(bvec)), max_iter, tol_rel, tol_abs, profile);
            otherwise
                error('Invalid LCP solver in params.');
        end
    end

    fprintf(['\n' parslv.colsolver ' LCP solution error = %e, iters = %d \n'],err,iter);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contact forces and modified densities
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if norm(lam) > 0 % Is at least one constraint active?
        % lam is the contact force magnitudes, so need to convert back to vector form
        % using the contact-force directions.
        F_c = F*lam;
        
        % Given the new data (i.e. force/torque) on surface, correct the densities
        % by doing another mobility solve.
        [mu_c,rho_c] = LOCAL_mobility_solve(TD, SD, Lk, Bk.', [], F_c, parslv);
    else % Else, there is no collision-related corrections that are needed.
        mu_c = [];
        rho_c = [];
        F_c = [];
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
    %{
    Calculates distance between spheres using their centers and radius.

    Inputs
        -

    Outputs
    %}
    n3 = size(C,1); 
    
    max_radii = max(Fparams.parbd.equ_radii, Fparams.parbd.polar_radii);
    diam = 2 * max_radii; 
    collision_eps = Fparams.parbd.collision_eps;  
    
    if n3 > 1
        % Compute center distances
        distC = get_distances_between_centers(C);
    
        % Find pairs for which (C_i-C-j) <= (r_i+r_j)+1.1*collision_eps*max(r_i,r_j)
        [ii,jj]=meshgrid(1:n3); 
        %this creates all combinations of indices (i, j) for i,j=1,...,n3 (all bodies)
        %finds all the indices whose distance is less than diam (r_i + r_j) + buffer which is relative as its in terms of 1.1*collision_eps*max(r_i, r_j)
        id = distC<=diam+1.1*collision_eps*max_radii & ii<jj; 
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

    colevent = mindst < 1.1*collision_eps; 
    collist = [ip jp]; 
end

function [colevent, collist, mindst, distances, closest_points_1, closest_points_2] = ...
        LOCAL_check_collision(C, Mt, Fparams)
    % C - centers of bodies (This is an double array, n_b x 3)
    % Mt - rotation matrices of bodies at current time

    % Do a quick check using the spheroids' circumscribed spheres.
    [colevent,collist,~] = LOCAL_check_collision_sph(C,Fparams);
    if ~colevent
        colevent = false;
        mindst = Inf; % Minimum distance
        distances = [];
        closest_points_1 = [];
        closest_points_2 = [];
        return;
    end
    
    % we are going to implement a not vectorized version for now, just to try and get this to work (and we will assume this is only for spheroids)

    [distances, closest_points_1, closest_points_2] = LOCAL_spheroidal_distances(C, Mt, Fparams, collist);

    % Compute refined distances per collision candidate pair.
    collision_eps = Fparams.parbd.collision_eps;
    dist = distances(:);

    % Set the new collision state using the finer scale.
    mindst = min(dist);
    active = dist < 1.1*collision_eps; % Match the 1.1 factor from above.
    collist = collist(active,:);
    distances = distances(active).';
    closest_points_1 = closest_points_1(:,active);
    closest_points_2 = closest_points_2(:,active);
    colevent = ~isempty(collist);
end

%% Utility functions
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
    %{
    Linear solve for Ax = b, with parameters given in parslv.
    %}
    prec = parslv.prec; 
    
    if nargin<6
        if ~isempty(prec) 
            pr=prec;  
        else
            pr=[]; 
        end
    end
    
    if isnumeric(A)
        x=A\b; 
    else
        x = zeros(size(b)); 
        for i=1:size(b,2)
            if nargin==2
                [x(:,i),~,rs,it]=gmres(A,b(:,i),1,1e-6,200,pr);
            else
                [x(:,i),~,rs,it]=gmres(A,b(:,i),parslv.rst,parslv.tol,parslv.maxit,pr);
            end
            fprintf('\n gmres %d iters=%d, res=%1.4g \n',i,prod(it),rs); 
        end
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

function q = rotation2quaternion(R)
    % Convert rotation matrix to unit quaternion
    % Given by Gemini, I need to verify this

    % numerical stability purposes
    
    v_w = 1 + trace(R);
    v_x = 1 + R(1,1) - R(2,2) - R(3,3);
    v_y = 1 - R(1,1) + R(2,2) - R(3,3);
    v_z = 1 - R(1,1) - R(2,2) + R(3,3);

    [~, max_index] = max([v_w, v_x, v_y, v_z]);
    switch max_index
        case 1
            qw = 0.5 * sqrt(v_w);
            qx = (R(3,2) - R(2,3)) / (4 * qw);
            qy = (R(1,3) - R(3,1)) / (4 * qw);
            qz = (R(2,1) - R(1,2)) / (4 * qw);
        case 2
            qx = 0.5 * sqrt(v_x);
            qw = (R(3,2) - R(2,3)) / (4 * qx);
            qy = (R(1,2) + R(2,1)) / (4 * qx);
            qz = (R(1,3) + R(3,1)) / (4 * qx);
        case 3
            qy = 0.5 * sqrt(v_y);
            qw = (R(1,3) - R(3,1)) / (4 * qy);
            qx = (R(1,2) + R(2,1)) / (4 * qy);
            qz = (R(2,3) + R(3,2)) / (4 * qy);
        case 4
            qz = 0.5 * sqrt(v_z);  
            qw = (R(2,1) - R(1,2)) / (4 * qz);
            qx = (R(1,3) + R(3,1)) / (4 * qz);
            qy = (R(2,3) + R(3,2)) / (4 * qz);
    end
    q = [qw; qx; qy; qz];
end

function R = quaternion2rotation(q)
    % Convert unit quaternion to rotation matrix
    % Given by Gemini, I need to verify this

    qw = q(1);
    qx = q(2);
    qy = q(3);
    qz = q(4);

    R = [
        1 - 2*(qy^2 + qz^2), 2*(qx*qy - qw*qz), 2*(qx*qz + qw*qy);
        2*(qx*qy + qw*qz), 1 - 2*(qx^2 + qz^2), 2*(qy*qz - qw*qx);
        2*(qx*qz - qw*qy), 2*(qy*qz + qw*qx), 1 - 2*(qx^2 + qy^2)
    ];
end

function psi = construct_global_psi(q)
    % This constructs the psi matrix that interfaces between angular velocity and derivative of the rotational configuration in quaternion form
    % This assumes angular velocity is given in global frame.
    P = [0 -q(4) q(3);
         q(4) 0 -q(2);
         -q(3) q(2) 0];

    psi = (1/2).*[-q(2:4).' ; 
                  q(1).*eye(3) - P];
end

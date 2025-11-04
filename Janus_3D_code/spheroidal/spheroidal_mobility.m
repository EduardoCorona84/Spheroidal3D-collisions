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
            The convention for the codebase is that if equi_radius > polar_radius, it is a oblate. Otherwise,
            it's an prolate.
        shape_type  - (string) n_b x 1 array of strings: should be 'prolate', 'oblate', or 'sphere'.
            Spheres are not implemented for the time being.
        p           - (int) spheroidal harmonic order (bodies) 
        Ct          - (double) n_b x 3 array of centers 
        eps         - (double) epsilon buffer (collision dist)
        mdist       - (double) collision buffer for body-body interactions
        out         - (bool) external vs internal evaluation (set to 1) 
    
    parslv - (struct) linear solver parameters such as 
         prec     - (string) preconditioner type, '' for unprec, 'bkdiag'
                    (block diagonal), 'TT' (tensor train)
         prtype   - (string) 'bkdiag' or 'TT'
         solver   - gmres, pcg, bicg, etc. 
         tol (tolerance), maxit (maximum iterations), rst (restart), etc.
    
    Depending on the type of problem that is implemented
    Ffun, Tfun = @(t,C,q) with output of size 3 x n_b.
        Force and torque prescriptions; required if type is 'FTfun'.
        The parameters are: t for timestep, C for the center of the body, and
        q for the quaternion associated with the body (to represent orientation).
    
    init    - (string) optional filename to resume a simulation from last
    recorded timestep
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
        % Generic property validation

        % parbd validation
        assert()
    end

    %%(0.1) (optional) Load data in init, initialize output arrays
    num_timesteps = Fparams.Nt; num_body=Fparams.parbd.n3; 
    tt = zeros(num_timesteps+1, 1); 
    sigma = cell(num_timesteps+1, 1); mu=sigma; U=sigma; VW=U; Xt=U; Ct=Xt; FT=mu; psi_Lap = mu; 
    Mt = cell(num_timesteps+1,num_body);
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
        sigma = cell(num_timesteps+1,1); mu=sigma; U=sigma; VW=U; Xt=U; Ct=Xt; FT=mu;
        Mt = cell(num_timesteps+1,num_body);
            
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
    timedisc = Fparams.tdisc; Sc = Fparams.parbd.Sc; 
    
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
    %%%%%%%% NEED TO UPDATE THIS FOR SPHEROIDS
    [colevent,collist,~,~] = LOCAL_check_collision_sph(C0,Fparams);
    
    if ~strcmp(Fparams.parbd.Shape,'') % For non-spherical shapes...
        % Seems like the purpose of this is to create a point cloud on the surface of each body.
        Sc2 = SurfaceSph(rad*shape_gallery(2*p,Fparams.parbd.Shape)); % rad (i.e. radius) is undefined...

        % Model surface pts
        % The dimension of X2 is (np2 * num_body) * 3
        % Each row contains the Cartesian coordinates of the points on the point clouds
        X2 = repmat(reshape(Sc2.cart.to_array,[],3),num_body,1);
        np2 = 2*(2*p)*(2*p+1);    
    else  
       X2=[];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (0.5) Timestepping (i.e. the main loop)
    Xt{1}=Xrp; Ct{1}=Fparams.parbd.C;
    
    % Linearize differential variational inequality each timestep
    for i=1:num_timesteps
        dt = dt0;
        if strcmp(timedisc,'euler')
            fprintf('\n ---------------------------------------------------------- \n')
            fprintf('\n (1) Explicit euler step \n')
            [Xt{i+1},Mt(i+1,:),Ct{i+1},U{i},FT{i},sigma{i},mu{i},VW{i},Kernels,Nullsp,...
                Fparams,colevent,collist,dt,psi_Lap{i},Energy(i)] = ...
                LOCAL_euler_step(Xt{i},Xt{1},X2,Mt(i,:),Ct{i},Kernels,Nullsp,Fparams,colevent,collist,t,dt,i);

        elseif strcmp(timedisc,'trapz')
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
        fprintf('\n dt: %2.2f ',dt)

        % Save progress every other step.
        if mod(i,2)==1                                                                                                                                              
            save(fname,'-v7.3','tt','Xt','Mt','Ct','FT','sigma','mu','U','VW','psi_Lap','Energy');                                                                                         
        else                                                                                                                                                        
            save([fname '2'],'-v7.3','tt','Xt','Mt','Ct','FT','sigma','mu','U','VW','psi_Lap','Energy');                                                                                   
        end  
        save([fname '_profile'],'timings'); 
    end %% END for linearization
end
    
%% Mobility solver system code
function [Xtp,Mtp,Ctp,U,FT,sigma,mu,VW,Kernels,Nullsp,Fparams,colevent,collist,dt,psi_Lap,Energy] = LOCAL_euler_step(Xt,X0,X2,Mt,Ct,Kernels,Nullsp,Fparams,colevent,collist,t,dt,it)
    %{
    Given the BIE matrices at timestep t and the known boundary conditions, calculate 
    the velocities at timestep t+1 and advance all bodies with the new velocities.

    Inputs
    

    Outputs
    Xtp - 
    Mtp -
    Ctp -
    U -
    FT -
    sigma -
    mu -
    VW -
    Nullsp
    Fparams - (struct) parameters of problem
    colevent - (boolean) did a collision event happen in advancing the timestep?
    collist
    dt
    psi_Lap - unused variable
    Energy - unused variable
    %}
    
    global timings; 
    Sc = Fparams.parbd.Sc; 
    np = Fparams.parbd.np; 
    n3 = Fparams.parbd.n3; 
    
    MRot = @(wh,t) RotationMat(wh,t);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % Get incoming force distribution: 
    tic; 
    [FT,sigma,VW,psi_Lap,Energy] = LOCAL_get_incoming_Fc(Fparams,t,dt,Kernels,Nullsp,Xt,Sc); 
    Ct
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
    % Check for collision after moving centers
    tic; 
    [colevent,collist,dt,Ctp, closest_points_1, closest_points_2] ...
    = LOCAL_collision_info(Fparams,X2,Ct,Ctp,VW,MRot,Mt,dt);
    fprintf('\n Time for collision detection: %e',toc);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Advance rotation matrix Mt and X
    tic; 
    [Mtp,Xtp,normW] = LOCAL_advance_rotation(MRot,VW,Mt,Xt,X0,dt,np,n3); 
    fprintf('\n Time to advance R(t) and X(t): %e',toc)
    timings.advance(it) = timings.advance(it) + toc; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Update Sc and operators
    fprintf('\n Surface and operator update')
    [Kernels,Nullsp,Fparams,timings] = RBS_Update_Operators(Xtp,Ctp,Mtp,normW,Kernels,Fparams,timings,it);
    timings.operator.total(it) = timings.operator.total(it) + timings.operator.surf(it) + timings.operator.diag(it) + timings.operator.offd(it);
    fprintf('\n Time to update surface and operators: %e',timings.operator.total(it));
end
    
function [Xtp,Mtp,Ctp,Kernels,Nullsp,Fparams,colevent,collist,dt] = LOCAL_advance_step(VW,mu,sigma,Xt,X0,X2,Mt,Ct,Kernels,Fparams,dt,it)
    %{

    %}
    global timings;
    
    np = Fparams.parbd.np; 
    n3 = Fparams.parbd.n3; 
    
    MRot = @(wh,t) RotationMat(wh,t);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Advance center Ct 
    tic; 
    Ctp = LOCAL_advance_center(Ct,dt,VW); 
    fprintf('\n Time to advance centers C(t): %e',toc); 
    timings.advance(it) = timings.advance(it) + toc; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check for collision after moving centers
    tic; 
    [colevent,collist,dt,Ctp, closest_points_1, closest_points_2] ...
    = LOCAL_collision_info(Fparams,X2,Ct,Ctp,VW,MRot,Mt,dt);
    fprintf('\n Time for collision detection: %e',toc);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Advance rotation matrix Mt and X
    tic; 
    [Mtp,Xtp,normW] = LOCAL_advance_rotation(MRot,VW,Mt,Xt,X0,dt,np,n3); 
    fprintf('\n Time to advance R(t) and X(t): %e',toc)
    timings.advance(it) = timings.advance(it) + toc; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update operators
    fprintf('\n Surface and operator update')
    [Kernels,Nullsp,Fparams,timings] ...
    = RBS_Update_Operators(Xtp,Ctp,Mtp,normW,Kernels,Fparams,timings,it);
    timings.operator.total(it) = timings.operator.total(it) + timings.operator.surf(it) + timings.operator.diag(it) + timings.operator.offd(it);
    fprintf('\n Time to update operators: %e',timings.operator.total(it));
end
    
function [FT, fM, VW, Energy] = LOCAL_get_incoming_Fc(Fparams,t,dt,Kernels,Nullsp,Xt,Sc)
    %{
        This functions mostly seems to be for debugging and not too relevant for the
        actual mobility solver.
    %}
    Energy = 0;
    parslv = Fparams.parslv; 
    lprec=[]; 
    acc = parslv.tol; rst = parslv.rst; maxit = parslv.maxit;
    
    %W = Fparams.parbd.W; 
    n3=Fparams.parbd.n3; p = Fparams.parbd.p; np = Fparams.parbd.np; 
    VW=[]; C = Fparams.parbd.C; 
    Bk = Nullsp.B; Ck = Nullsp.C; Lk = Nullsp.L; 
    
    switch Fparams.type
    case 'FTfun' 
        Ffun = Fparams.Ffun; 
        Tfun = Fparams.Tfun; 
        Ct = Fparams.parbd.C; 
        
        Force = Ffun(t,Ct);  
        Torque = Tfun(t,Ct);
        FT = [Force;Torque]; FT = FT(:); 
        fM = Bk'*FT;
    end
end

function [sigma,mu,U,VW] = LOCAL_compute_velocities(sigma,VW,Ct,Kernels,Nullsp,Fparams,col,collist,i,dt, Mt, closest_points_1, closest_points_2)
    %{
    Computes rigid body velocities and advance centroids and rotation matrices.

    Inputs

    Outputs
        sigma : density (not at all necessary!)
        mu : 
        VW : rigid body velocities (translational/rotational velocities)
    %}
    global timings; 
    parslv = Fparams.parslv; 
    rd = Fparams.parbd.rd; tau = Fparams.parbd.tau; W = Fparams.parbd.W; 
    num_body = Fparams.parbd.n3; 
    rdt = repmat((rd.').^(-4),3,1); 
    rdw = repmat((rd.').^(-2),3,1);
    vind = reshape(repmat(6*(0:num_body-1),3,1),1,[])+repmat((1:3),1,num_body);
    wind = reshape(repmat(6*(0:num_body-1),3,1),1,[])+repmat((4:6),1,num_body);
    
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

        % U = U_inc + U_sc
        tic; 
        U = Lapp(Kernels.SD,(mu+sigma)); 
        fprintf('\n Time for apply (of S) to compute U: %e',toc);  
        timings.velocities.apply(i) = timings.velocities.apply(i) + 0.5*toc;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute (V,W) and advance Ct, Xt and Mt
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic;
        % The operator C here seems rather mysterious.
        % We need to rebuild the rigid body motions since we have updated it.
        CU = Nullsp.C*U;
        IU = CU(vind); 
        WxI = CU(wind); 
            
        % U(S_k) = V_k + W_k x (X(S_k) - C_k)
        VW = zeros(6,num_body); 
        if ~iscell(W)
            VW(1:3,:) = (1/sum(W))*(rdw.*reshape(IU,3,num_body)); 
            VW(4:6,:) = tau\(rdt.*reshape(WxI,3,num_body)); 
        else
            IUv = reshape(IU,3,num_body); WxIv = reshape(WxI,3,num_body); 
            for j=1:num_body
                VW(1:3,j) = (1/sum(W{j}))*IUv(:,j); 
                VW(4:6,j) = tau{j}\WxIv(:,j);
            end
        end
        fprintf('\n Time to compute V and W: %e',toc); 
        timings.velocities.vw(i) = toc; 
        
        errbs = norm(U-Nullsp.D'*VW(:))/norm(U);
        fprintf('\n Rigid body velocity error: %e',errbs)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Collision event
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif col
        tic; 
        B = Nullsp.L*sigma-Lapp(Kernels.TD,sigma);
        timings.velocities.apply(i) = toc;
        
        % Solve Fredholm eq TD*mu = B
        tic; 
        mu = Lslv(Kernels.TD,B,parslv); 
        timings.velocities.solve(i) = toc;

        tic;  
        fprintf('\n-------------------------------------------------');
        fprintf('\n Computing contact force at collision sites: \n')
        display(collist(:,1:2)')
        fprintf('-------------------------------------------------\n');
        
        % Compute contact force and force distribution updates
        %we might have to compute distances before the first pass
        [F_c,mu_c,rho_c] = LOCAL_Compute_Contact_LCP(collist,Kernels,Nullsp,Fparams,Ct,VW,dt, Mt,closest_points_1, closest_points_2);
        
        F_c = reshape(F_c,6,[]);
        display(F_c(1:3,:));
        
        if ~isempty(mu_c)
            %Update sigma, mu, U and VW
            mu = mu + mu_c; 
            sigma = sigma + rho_c; 
            U = Lapp(Kernels.SD,(mu+sigma)); 
            
            CU  = Nullsp.C*U; 
            IU  = CU(vind); 
            WxI = CU(wind); 
        
            if ~iscell(W)
                VW(1:3,:) = (1/sum(W))*(rdw.*reshape(IU,3,num_body)); 
                VW(4:6,:) = tau\(rdt.*reshape(WxI,3,num_body)); 
            else
                IUv = reshape(IU,3,num_body); WxIv = reshape(WxI,3,num_body); 
                for j=1:num_body
                    VW(1:3,j) = (1/sum(W{j}))*IUv(:,j); 
                    VW(4:6,j) = tau{j}\WxIv(:,j);
                end
            end
        end
    
        fprintf('\n Time for contact force correction: %e',toc);  
        timings.velocities.col(i) = toc;
            
        % Check rigid body velocity
        errbs = norm(U-Nullsp.D'*VW(:))/norm(U);
        fprintf('\n Rigid body velocity error after collision correction: %e',errbs)
    end
    
    VW = real(VW); 
end
    
function Ctp = LOCAL_advance_center(Ct,dt,VW)
    Ctp = Ct + dt*VW(1:3,:)';
end

function [Mtp,Xtp,normW] = LOCAL_advance_rotation(MRot, VW, Mt, Xt, X0, dt, np, num_body)
    normW = zeros(num_body,1); Mtp=Mt; Xtp=Xt;
    for k=1:num_body
        xind = (1:np)+np*(k-1); 
        normW(k) = norm(VW(4:6,k));

        if normW(k) > ROTATIONAL_VELOCITY_TOL % Checks if significant enough to update.
            Mtp{k} = MRot(VW(4:6,k),dt)*Mt{k};          
            Xtp(xind,:) = X0(xind,:)*Mtp{k}';
        else   
            Mtp{k} = Mt{k}; 
            Xtp(xind,:) = Xt(xind,:);    
        end   
    end
end

%% Collision resolution
function [F_c,mu_c,rho_c] = LOCAL_Compute_Contact_LCP(collist,Kernels,Nullsp,Fparams,Ct,VW,dt, Mt, distances, closest_points_1, closest_points_2)
    %{
    Calculates the force and the resulting densities due to the collision.
    Note that this does not resolve the collision, only calculates the effect of it.

    Inputs:
    collist
    Kernels
    Nullsp - 
    Fparams - parameters of mobility solver
    Ct - centers of bodies
    VW - translational/rotational velocities
    dt - timestep size
    Mt - rotation matrices of bodies
    closest_points_1 - closest points on body 1 of each colliding pair
    closest_points_2 - closest points on body 2 of each colliding pair

    Outputs:
    F_c -
    mu_c
    rho_c
    %}
    parslv = Fparams.parslv; 
    rd = Fparams.parbd.rd; 
    n3 = length(rd); 
    diam = Fparams.parbd.diam; %diam(i,j) = r_i + r_j
    mxrd = Fparams.parbd.mxrd; %max(r_i,r_j)
    eps = Fparams.parbd.eps;
    Nb = Fparams.parbd.Nb; 
    tol = parslv.tol; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %{
    ip = collist(:,1); jp = collist(:,2); 
    numFS = 0; 
    %}
    %

    TD = Kernels.TD; SD = Kernels.SD; 
    Bk = Nullsp.B; Ck = Nullsp.C; Lk = Nullsp.L; 
    numF = length(ip); 



    % Compute vectors and normal vectors for pairs
    %this is under the assumption that we have spheres. This will need to be modified for spheroids
    

    R = Ct(ip,:)-Ct(jp,:);       %Ci - Cj numF x 3
    NR = sqrt(sum(R.*R,2));      %|Ci-Cj| numF x 1 
    Rhat = repmat(1./NR,1,3).*R; %eij = (Ci - Cj)/|Ci-Cj|
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Build A
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F = zeros(6*n3,numF+numFS); 
    %I'm not sure what numFS is for. This looks like its constructing the "D" matrix from my notes. 

    %We also probably need the "quaternion matrix", the matrix that converts the angular velocity to the time derivative of the configuration. In spherical case this is just the identity, but we will have more than this.
    %NOTE: In this example we will be using rotation matrices, not quaternions.

    %for all the colision pairs, we fill in the normal directions into the F matrix
    %We also need to do this for torques, which is not done yet
    %The indexing of F is giving by the indexing of the collision pairs.
    %F^T is size: the number of collision pairs by 6 times the number of particles (in the collision list/event)
    %iterate through all collision pairs (order given by collist)
    for k=1:numF
        %normal (translational) forces
        translational_indi = (1:3)+6*(ip(k)-1);
        translational_indj = (1:3)+6*(jp(k)-1);

        rotational_indi = (4:6)+6*(ip(k)-1);
        rotational_indj = (4:6)+6*(jp(k)-1);
        %The normal vector is (by convention) pointing from particle 2 to particle 1 (j to i). The signs need to be consistent with this convention.

        %compute normal vector for collision pairs
        if distances(k) < separation_tol
            % if the particles are close enough, we will use the normal of the second particle in the collision pair
            ellipsoid_mat = Mt{jp(k)}*diag([Fparams.parbd.polar_radii(jp(k)) Fparams.parbd.equi_radii(jp(k)) Fparams.parbd.equi_radii(jp(k))].^(-2))*Mt{jp(k)}';
            %make sure dimensionality is correct here (this is something for me to do in general)
            normal = ellipsoid_mat*(closest_points_2(:,k)'-Ct(jp(k),:)');
            
        else
            normal = (closest_points_2(:,k)'-closest_points_1(:,k)')/distances(k);
        end

        %compute torques for collision pair 
        local_coordinates_1 = closest_points_1(:,k)'-Ct(ip(k),:)';
        local_coordinates_2 = closest_points_2(:,k)'-Ct(jp(k),:)';
        torque_direction_1 = -cross(local_coordinates_1, normal);
        torque_direction_2 = cross(local_coordinates_2, normal);

        F(translational_indi, k) = -normal;
        F(translational_indj, k) = normal;

        F(rotational_indi, k) = torque_direction_1;
        F(rotational_indj, k) = torque_direction_2;
    end

    % We need some way to store the rotational part of the configuration as a vector. 
    %By default, things are stored as rotation matrices. Rotation matrices can be stored as a vector, but I'm not sure what the "interfacing" matrix/operator would be from a vectorized rotation matrix to the angular velocity.
    %Because I have already done all of this with quaternions, I think I am going to just convert to quaternions for now and in the future we can explore alternative approaches.


    %I'm not sure what this is doing 
    %{
    for k=numF+1:numF+numFS
    indi = (1:3)+6*(ipsh(k-numF)-1);
    F(indi,k) = -Ct(ipsh(k-numF),:)./norm(Ct(ipsh(k-numF),:));
    end
    %}

    % A is built using the dense or matfree mobility matrix. Can be accelerated
    % by employing only self interaction (block-diagonal) for TD and SD. 
    % Different criteria can be added as needed
    matfree = ~Fparams.denseMV || numF+numFS > 1; 
    bkdiag=false;

    if ~matfree
        Bf = (Bk.')*F; 
        %(3) (-0.5I-K)*rho_c
        MNS = -Lapp(TD,Bf)+Lk*Bf;
        %(4) 3x3 MNS=VNS*S(mu_c+rho_c)  
        MuNS = Lslv(TD,MNS,parslv); 

        % Setup LCP x perp A*x + b (dense build of Amat = F^T M F)
        Amat = real(F.'*(Ck*Lapp(SD,MuNS+Bf))); 
    else
        parslv.tol = parslv.coltol; 
        Bf = @(x) (Bk.')*(F*x);
        if bkdiag
            S0 = @(x) reshape(Kernels.SSD0*(repmat(rd.',Nb,size(x,2)).*reshape(x,Nb,n3*size(x,2))),[],size(x,2));
            IT0 = @(x) reshape(Kernels.ITSSD0*reshape(x,Nb,n3*size(x,2)),[],size(x,2));
            Amat = @(x) real(F.'*(Ck*(S0(-IT0(Lapp(TD,Bf(x))+Lk*Bf(x))+Bf(x)))));
        else 
            Amat = @(x) real(F.'*(Ck*Lapp(SD,Lslv(TD,-Lapp(TD,Bf(x))+Lk*Bf(x),parslv)+Bf(x))));
        end
        
        parslv.tol = tol; 
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Build constant vector b
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute (1/dt)*phi
    %we will also need the quanternion part here as well (I think)
    phib = zeros(numF+numFS,1); 
    if numF>0
        phib(1:numF) = (1/dt)*(NR-diam(ip+n3*(jp-1))-eps*mxrd(ip+n3*(jp-1))); 
    end

    %b_k = (1/dt)*phi_k + F.'V_k
    bvec = phib + real((F.')*VW(:));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %LCP solve
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %LCP params
    max_iter=parslv.colmaxit; 
    tol_rel=parslv.col_tolrel; %1e-6; 
    tol_abs=parslv.col_tolabs; %1e-9; 
    profile=1;

    if matfree
        switch parslv.colsolver
            case 'Newton'
            % solve LCP using minmap Newton (matfree)
            [lam ,err ,iter, ~, ~, ~] = ...
            minmap_newton_matfree(Amat, bvec, zeros(size(bvec)), max_iter, tol_rel, tol_abs, profile );
            case 'APGD'
            % solve LCP using Accelerated PGD
            [lam ,err ,iter, ~, ~, ~] = ...
            APGD_matfree(Amat, bvec, zeros(size(bvec)), max_iter, tol_rel, tol_abs, profile );    
            case 'BBPGD'
            % solve LCP using Barzilai Borwein PGD
            [lam ,err ,iter, ~, ~, ~] = ...
            BBPGD_matfree(Amat, bvec, zeros(size(bvec)), max_iter, tol_rel, tol_abs, profile );    
            otherwise
            % solve LCP using Barzilai Borwein PGD
            [lam ,err ,iter, ~, ~, ~] = ...
            BBPGD_matfree(Amat, bvec, zeros(size(bvec)), max_iter, tol_rel, tol_abs, profile );    
        end
    else
        switch parslv.colsolver
            case 'Newton'
            % solve LCP using minmap Newton  
            [lam ,err ,iter, ~, ~, ~] = ...
            minmap_newton(Amat, bvec, zeros(size(bvec)), max_iter, tol_rel, tol_abs, profile );
            case 'APGD'
            % solve LCP using Accelerated PGD
            [lam ,err ,iter, ~, ~, ~] = ...
            APGD(Amat, bvec, zeros(size(bvec)), max_iter, tol_rel, tol_abs, profile );    
            case 'BBPGD'
            % solve LCP using Barzilai Borwein PGD
            [lam ,err ,iter, ~, ~, ~] = ...
            BBPGD(Amat, bvec, zeros(size(bvec)), max_iter, tol_rel, tol_abs, profile );    
            otherwise
            % solve LCP using Barzilai Borwein PGD
            [lam ,err ,iter, ~, ~, ~] = ...
            BBPGD(Amat, bvec, zeros(size(bvec)), max_iter, tol_rel, tol_abs, profile );    
        end
    end

    fprintf(['\n minmap ' parslv.colsolver ' LCP solution error = %e, iters = %d \n'],err,iter);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contact forces and modified densities
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if norm(lam)>0
        % Contact forces / torques
        F_c = F*lam; 
        % Obtain rho and mu densities
        [mu_c,rho_c] = Lapp_ctmat(TD,SD,Lk,Bk.',[],F_c,parslv);
    else
        mu_c=[]; rho_c=[]; F_c=[];  
    end
end

function [colevent,collist, dt,Ctp, closest_points_1, closest_points_2,] = LOCAL_collision_info(Fparams,X2,Ct,Ctp,VW,MRot,Mt,dt)
    %{
    
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
    n3 = size(Ct,1); Shape=Fparams.parbd.Shape; 
    eps = Fparams.parbd.eps;  
    
    % Check for collision between spheres
    [colevent,collist,mindst,mindstsh] = LOCAL_check_collision_sph(Ctp,Fparams);
    
    % Finer collision detection for non-spheres
    if ~strcmp(Shape,'')
        np2 = size(X2,1)/n3;
        [colevent,collist,mindst,closest_points_1,closest_points_2]=LOCAL_check_collision(collist,eps,Ctp,X2,MRot,Mt,VW,dt,np2); 
    end
    
    if colevent
        % If mindst<<eps, or <0, we need to adjust timestep
        bis=0;
        maxbis=3;
        
        cond = mindst < 0.1*eps;

        while cond && bis<=maxbis
        
            % Recompute dt and Ct{i+1} to avoid collision
            dt = dt/2;
            bis=bis+1; 
            Ctp = LOCAL_advance_center(Ct,dt,VW); 
        
            %Check for collision between spheres (or sphere envelopes)
            [colevent,collist,mindst,mindstsh] = LOCAL_check_collision_sph(Ctp,Fparams);
    
            % Finer collision detection for non-spheres
            if ~strcmp(Shape,'')
                %Xip and Xjp are not used here, but I think they should be used for the closest points. We should change LOCAL_collision_info to return these as these will be needed in the optimization problem formulation.
                [colevent,collist,mindst,Xip,Xjp]=LOCAL_check_collision(collist,eps,Ctp,X2,MRot,Mt,VW,dt,np2); 
            end
            
            fprintf('\n bisection = %d: dt = %1.4e, mindst = %1.4e, mindstsh = %1.4e',bis,dt,mindst,mindstsh);
            
            % update condition
            cond = mindst < 0.1*eps;
        end
        
        fprintf('\n Min pairwise relative distance after collision detection: %2.4f ',mindst);
    else
        collist=[]; 
        fprintf('\n Min pairwise relative distance: %2.4f',mindst);
    end
end
    
function [colevent,collist,mindst,mindstsh] = LOCAL_check_collision_sph(C,Fparams)
    %{
    Calculates distance between spheres using their centers and radius.

    Inputs
        -

    Outputs
    %}
    n3 = size(C,1); 
    
    rd = Fparams.parbd.rd; 
    mxrd = Fparams.parbd.mxrd; 
    diam = Fparams.parbd.diam; 
    eps = Fparams.parbd.eps;  
    
    if n3 > 1
        % Compute center distances
        distC = LOCAL_CenterDistances(C);
    
        % Find pairs for which (C_i-C-j) <= (r_i+r_j)+1.1*eps*max(r_i,r_j)
        [ii,jj]=meshgrid(1:n3); 
        %this creates all combinations of indices (i, j) for i,j=1,...,n3 (all bodies)
        %finds all the indices whose distance is less than diam (r_i + r_j) + buffer which is relative as its in terms of 1.1*eps*max(r_i, r_j)
        id = distC<=diam+1.1*eps*mxrd & ii<jj; 
        %This selects the indices where the 
        ip = ii(id); 
        jp = jj(id); 
    
        % Compute minimum relative distance between spheres
        mindst=min(reshape((distC-diam)./mxrd,[],1));
    else
        ip=[]; jp=[]; 
        mindst=Inf; 
    end

    colevent = mindst < 1.1*eps; 
    mindstsh = Inf;
    collist = [ip jp]; 
end
    
function [colevent,collist,mindst,Xip,Xjp]=LOCAL_check_collision(collist,eps,C,X,MRot,Mt,VW,dt,np)
    %collist - list of potential colliding pairs assuming they are spheres (for this we would use the largest radii?, assuming spheroids)
    % eps - an epsilon buffer for collision detection
    % C - centers of bodies (This is an double array, n_b x 3)
    % X - surface points of bodies (don't think this is relevant for spheroids)
    % MRot - function handle for rotation matrix takes in a timestep, and angular velocity vector (why do we need a timestep?)
    % Mt - rotation matrices of bodies at current time
    % VW translational/rotational velocities
    % dt - timestep size
    % np - number of surface points per body (don't think is relevant for spheroids)

    %add a switch statment for different shapes (this will just be for spheroids for now)
    %in the spheroid case, we will check if collevent occurs if we modify LOCAL_check_collision_sph and LOCAL_collision_info (for computing distances we can use the max radii for each spheroid)
    
    % we are going to implement a not vectorized version for now, just to try and get this to work (and we will assume this is only for spheroids)

    switch(Fparams.parbd.Shape)
        case {'oblate', 'prolate'}
            [distances, closest_points_1, closest_points_2] = LOCAL_spheroidal_distances(C, MRot, collist, Fparams);
    end

    %now we determine which distances are less than the epsilon buffer
    indices = distances < %some sort of criteria, probably need spheroid specific
    collist = collist(indices, :);
    colevent = ~isempty(collist);
    mindst = min(distances);
    



end

%% Utility functions
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
    
function den = LOCAL_CenterDistances(C)
    %{
    Calculates pairwise Euclidean distance between points in R^3.

    Inputs
    C - (double) n_b x 1 centers of bodies

    Outputs
    den - distances between bodies
    %}
    [Y_g1,  X_g1] = meshgrid(C(:,1), C(:,1));
    [Y_g2,  X_g2] = meshgrid(C(:,2), C(:,2));
    [Y_g3,  X_g3] = meshgrid(C(:,3), C(:,3));
    d1 = (X_g1 - Y_g1); d2 = (X_g2 - Y_g2); d3 = (X_g3 - Y_g3); 
    den = sqrt(d1.^2 + d2.^2 + d3.^2);

    if any(den(:) == 0)
        error('Zero distance in LOCAL_CenterDistance; handling of this is not implemented.');
    end
end

function [distances, closest_points_1, closest_points_2] = LOCAL_spheroidal_distances(collist, C, Mt, Fparams)

    %switch statement for different distance algorithms.
    %right now we will just implement two options, moving balls and GJK signed volumes accelerated
    switch lower(Fparams.parbd.distance.algo)
        case 'moving balls'
            distance_algo = @moving_balls_pair;
        case 'gjk signed volumes accelerated'
            distance_algo = @GJK_signed_volumes_accelerated_pair;
        otherwise
            error('Unknown spheroidal distance algorithm %s', Fparams.parbd.distance.algo);
    end

    %row vector of the pairwise distances
    distances = zeros(1, size(collist, 2));
    %row vector of the closest points of the 1st particle in the collision pair
    closest_points_1 = zeros(3, size(collist, 2));
    %row vector of the closest points of the 2nd particle in the collision pair
    closest_points_2 = zeros(3, size(collist, 2));

    for collision_pair = 1:size(distances, 2)
        %this depends on how I implement the distance algorithms, which params are passed, the params are just placeholders for now, but we really just need the centers, the shapes (and rotation matrices), tolerance and iters

        %create temporary structs to hold the spheroid params for the 2 spheroids in the collision pair
        spheroid_1_params.C = C(collist(collision_pair, 1), :);
        spheroid_1_params.R = Mt{collist(collision_pair, 1)};
        spheroid_1_params.a = Fparams.parbd.polar_radii(collist(collision_pair, 1));
        spheroid_1_params.b = Fparams.parbd.equ_radii(collist(collision_pair, 1));
        spheroid_1_params.c = spheroid_1_params.b;

        %spheroid 2
        spheroid_2_params.C = C(collist(collision_pair, 2), :);
        spheroid_2_params.R = Mt{collist(collision_pair, 2)};
        spheroid_2_params.a = Fparams.parbd.polar_radii(collist(collision_pair, 2));
        spheroid_2_params.b = Fparams.parbd.equ_radii(collist(collision_pair, 2));
        spheroid_2_params.c = spheroid_2_params.b;

        [closest_points_1(:, collision_pair), closest_points_2(:, collision_pair), distances(collision_pair)] = distance_algo(spheroid_1_params, spheroid_2_params, Fparams.parbd.distance.tol, Fparams.parbd.distance.max_iters);
    end

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
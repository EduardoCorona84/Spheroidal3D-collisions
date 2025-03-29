function mobility_solver(fname,Fparams,init)
    %{
    Entrypoint for the mobility solver for spheroidal suspension in Stokes flow.
    Code was adapted from the code for Janus particles.
     
    Inputs: 
    fname - (string) filename for experiment info
    
    Fparams - (struct) parameter struct for rigid body simulation, with fields:
    
    type            - (string) mobility problem type (FTfun)
    denseMV         - (bool) dense vs FMM for far-field
    num_timesteps   - (int)    number of timesteps
    dt              - (double) timestep length
    timedisc        - (string) time discretization (i.e. "euler", "trapz", or "rk4")
    comp            - (bool) compute intermediate quantities FT and VW
    lambda          - (double) parameter for modified laplace case
    e_north/e_south - (double) Janus particle relative permittivity 

    parbd - (struct) struct with rigid body parameters:
        Shape       - (string) rigid body shape, '' is default for sphere
        n3          - (int) number of rigid bodies n_b
        rd          - (double) radius (monodisperse) or array of n_b radii (polydisperse) of bodies
        p           - (int)    spherical harmonic order (bodies) 
        Ct          - (double) n_b x 3 array of centers 
        eps         - (double) epsilon buffer (collision dist)
        mdist       - (double) collision buffer for body-body interactions
        out         - (bool) external vs internal evaluation (set to 1) 
    
    parslv - (struct) linear solver parameters such as 
         prec     - (string) preconditioner type, '' for unprec, 'bkdiag'
                    (block diagonal), 'TT' (tensor train)
         solver   - gmres, pcg, bicg, etc. 
         tol (tolerance), maxit (maximum iterations), rst (restart), etc.
           
    parsh - (struct) optional shell geometry parameters: 
        psh - (int)    spherical harmonic order
        shrd - (double) shell radius
        mdsh - (double) min relative distance to origin for near-sing (<1)
        epsh - (double) collision buffer with geometry
        out   - (bool) external vs internal evaluation (set to 0)
     
    The default 'FTfun' (force and torque prescription) requires functions 
    Ffun,Tfun = @(t,C) with output of size 3 x n_b. 
    
    init    - (string) optional filename to resume a simulation from last
    recorded timestep

    %%%
    %%% CODE ANNOTATIONS
    %%%
    Mt stores rotational information for n bodies at each time step
    Mt0 is the initial rotation setup

    VW represents the translational (v)/angular velocity (w) each timestep (as a 6 x num_body matrix)
    VW0 represents the initial velocities

    Xrp (now renamed to X_ref_pts) are the tracking points needed for collision? Doesn't seem to be used either way.

    Fparams.parbd.np is number of discretization points on the surface of a body (is the same for all bodies?)

    The code "xind = (1:np) + np*(k-1)" shows up often; it collects all of the discretization points for the body k.
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global timings;
    global ROTATIONAL_VELOCITY_TOL;     % The threshold that the norm of the angular velocity should meet in order
    ROTATIONAL_VELOCITY_TOL = 1e-10;    % to update the body's angular position.
     
    %(0.1) (optional) Load data in init, initialize output arrays
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
    Fparams = MS_initparams(Fparams);
    timings.setup_surace = toc;
    timedisc = Fparams.tdisc; Sc = Fparams.parbd.Sc; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %(0.3) Initialize kernels (for MatVecs) and nullspace info
    tic; 
    Xrp = Fparams.parbd.Xrp; C0 = Fparams.parbd.C; 
    np = Fparams.parbd.np; normW = zeros(num_body,1);  
    if init_flag
        for k=1:num_body
            normW(k) = norm(VW0(4:6, k));
            % Collect all discretization points for each body
            xind = (1:np) + np*(k-1);
            % Apply rotation (from loaded state) to each discretization point
            Xrp(xind,:) = Xrp(xind,:)*Mt0{k}';
        end 
    end
    
    Kernels=[]; 
    [Kernels,Nullsp,Fparams,timings] = RBS_Update_Operators(Xrp,C0,Mt0,normW,Kernels,Fparams,timings,0); 
    timings.setup_kernel=toc;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (0.4) Initialize collision info  
    [colevent,collist,~,~] = LOCAL_check_collision_sph(C0,Fparams);
    
    if ~strcmp(Fparams.parbd.Shape,'') % For non-spherical shapes...
        % Seems like the purpose of this is to create a point cloud on the surface of each body.
        Sc2 = SurfaceSph(rad*shape_gallery(2*p,Fparams.parbd.Shape)); % rad is undefined...

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
    end
end
    
function [Xtp,Mtp,Ctp,U,FT,sigma,mu,VW,Kernels,Nullsp,Fparams,colevent,collist,dt,psi_Lap,Energy] = LOCAL_euler_step(Xt,X0,X2,Mt,Ct,Kernels,Nullsp,Fparams,colevent,collist,t,dt,it)
    global timings; 
    Sc = Fparams.parbd.Sc; 
    diam = Fparams.parbd.diam; 
    mxrd = Fparams.parbd.mxrd; 
    Shape = Fparams.parbd.Shape; 
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
    [sigma,mu,U,VW] = LOCAL_compute_velocities(sigma,VW,Ct,Kernels,Nullsp,Fparams,colevent,collist,it,dt); 
    timings.velocities.total(it) = timings.velocities.total(it) + timings.velocities.solve(it) + timings.velocities.apply(it) ...
        + timings.velocities.vw(it) + timings.velocities.col(it); 
    fprintf('\n Time to compute velocities / fluid solve: %e',timings.velocities.total(it)); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Advance center Ct 
    tic; 
    Ctp = LOCAL_advance_center(Ct,dt,VW,Fparams); 
    fprintf('\n Time to advance centers C(t): %e',toc); 
    timings.advance(it) = timings.advance(it) + toc; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check for collision after moving centers
    tic; 
    [colevent,collist,dt,Ctp] ...
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
    if isfield(Fparams,'parsh')
    tic; 
    Fparams = LOCAL_compute_shell_velocity(mu+sigma,Fparams);
    fprintf('\n Time to compute boundary correction: %e',toc)
    timings.velocities.shell(it) = toc;
    timings.velocities.total(it) = timings.velocities.total(it) + timings.velocities.shell(it); 
    end
    
    [Kernels,Nullsp,Fparams,timings] ...
    = RBS_Update_Operators(Xtp,Ctp,Mtp,normW,Kernels,Fparams,timings,it);
    timings.operator.total(it) = timings.operator.total(it) + timings.operator.surf(it) + timings.operator.diag(it) + timings.operator.offd(it);
    fprintf('\n Time to update surface and operators: %e',timings.operator.total(it));  
    
end
    
function [Xtp,Mtp,Ctp,Kernels,Nullsp,Fparams,colevent,collist,dt] = LOCAL_advance_step(VW,mu,sigma,Xt,X0,X2,Mt,Ct,Kernels,Fparams,dt,it)
    global timings;  
    
    np = Fparams.parbd.np; 
    n3 = Fparams.parbd.n3; 
    
    MRot = @(wh,t) RotationMat(wh,t);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Advance center Ct 
    tic; 
    Ctp = LOCAL_advance_center(Ct,dt,VW,Fparams); 
    fprintf('\n Time to advance centers C(t): %e',toc); 
    timings.advance(it) = timings.advance(it) + toc; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check for collision after moving centers
    tic; 
    [colevent,collist,dt,Ctp] ...
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
    if isfield(Fparams,'parsh')
    tic; 
    Fparams = LOCAL_compute_shell_velocity(mu+sigma,Fparams);
    fprintf('\n Time to compute boundary correction: %e',toc)
    timings.velocities.shell(it) = toc;
    timings.velocities.total(it) = timings.velocities.total(it) + timings.velocities.shell(it); 
    end
    
    [Kernels,Nullsp,Fparams,timings] ...
    = RBS_Update_Operators(Xtp,Ctp,Mtp,normW,Kernels,Fparams,timings,it);
    timings.operator.total(it) = timings.operator.total(it) + timings.operator.surf(it) + timings.operator.diag(it) + timings.operator.offd(it);
    fprintf('\n Time to update surface and operators: %e',timings.operator.total(it));
end
    
    function [FT, fM, VW, Energy] = LOCAL_get_incoming_Fc(Fparams,t,dt,Kernels,Nullsp,Xt,Sc)
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
    
function y = Lapp(A,x)
    % Left-apply the matrix A to the vector x.
    if isnumeric(A)
        y=A*x; 
    else
        y=real(A(x)); 
    end
end
    
    function x = Lslv(A,b,parslv)
        % Linear solve for Ax = b, with metaparameters given in parslv.
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
            %b = parslv.Pr(b); 
            %A = @(x) parslv.Pr(A(x)) + 0.5*(x-parslv.Pr(x)); 
            
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
    
    function den = LOCAL_CenterDistance(C)
        % Calculates pairwise Euclidean distance between points in R^3.
        [Y_g1,  X_g1] = meshgrid(C(:,1), C(:,1));
        [Y_g2,  X_g2] = meshgrid(C(:,2), C(:,2));
        [Y_g3,  X_g3] = meshgrid(C(:,3), C(:,3));
        d1 = (X_g1 - Y_g1); d2 = (X_g2 - Y_g2); d3 = (X_g3 - Y_g3); 
        den = sqrt(d1.^2 + d2.^2 + d3.^2);

        if any(den(:) == 0)
            error('Zero distance in LOCAL_CenterDistance; handling of this is not implemented.');
        end
    end
    
    function [F_c,mu_c,rho_c] = LOCAL_Compute_Contact_LCP(collist,Kernels,Nullsp,Fparams,Ct,VW,dt)
        error("LOCAL_Compute_Contact_LCP is not implemented for spheroids.");
    end
    
    function [sigma,mu,U,VW] = LOCAL_compute_velocities(sigma,VW,Ct,Kernels,Nullsp,Fparams,col,collist,i,dt)
        %global rst maxit acc timings rd; 
        global timings; 
        parslv = Fparams.parslv; 
        rd = Fparams.parbd.rd; tau = Fparams.parbd.tau; W = Fparams.parbd.W; 
        num_body = Fparams.parbd.n3; 
        rdt = repmat((rd.').^(-4),3,1); 
        rdw = repmat((rd.').^(-2),3,1);  
        shflg = isfield(Fparams,'parsh'); 
        vind = reshape(repmat(6*(0:num_body-1),3,1),1,[])+repmat((1:3),1,num_body);
        wind = reshape(repmat(6*(0:num_body-1),3,1),1,[])+repmat((4:6),1,num_body);
        
        % If inside shell, add traction from boundary correction to rhs
        Fparams.Tshell = []; 
        if shflg
            if ~isempty(Fparams.parsh.shellden)
                fprintf('\n Evaluation of traction correction field from shell:')
                Fparams.parsh.flag_pot = 'TSL_Stk_3D'; 
                Xtrg = Fparams.parbd.Xp; Nrtrg = Fparams.parbd.Nrp; 
                tic; 
                Fparams.Tshell = real(VSh_MatVec_RB_trg(Fparams.parsh.shellden,Xtrg,Nrtrg,Fparams.parsh));
                fprintf('\n Time for Traction computation from boundary: %e',toc); 
                
                parshS = Fparams.parsh; parshS.flag_pot = 'SL_Stk_3D'; 
                Ush = real(VSh_MatVec_RB_trg(parshS.shellden,Xtrg,Nrtrg,parshS)); 
            end
        end
        
        if isempty(VW) || Fparams.comp
            % Fluid Solve
            % (1) U_inc=S[sigma] (particular solution given forces and torques)
            % Sigma is the incoming traction distribution 
            
            % (2) U_sc=S[mu] ("scattered" field with zero forces and torques)
            % RHS -(aI+K)*sigma
            tic; 
            B = Nullsp.L*sigma-Lapp(Kernels.TD,sigma);
            timings.velocities.apply(i) = 0.5*toc;
            
            if ~isempty(Fparams.Tshell)
                B = B - Fparams.Tshell; 
            end
            
            % Solve Fredholm eq TD*mu = B
            tic; 
            mu = Lslv(Kernels.TD,B,parslv);
            fprintf('\n Time for solve: %e',toc); 
            timings.velocities.solve(i) = toc;  

            % U = U_inc + U_sc
            tic; 
            U = Lapp(Kernels.SD,(mu+sigma)); 
            fprintf('\n Time for apply (of S) to compute U: %e',toc);  
            timings.velocities.apply(i) = timings.velocities.apply(i) + 0.5*toc;
            
            if ~isempty(Fparams.Tshell)
                U = U + Ush; 
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Compute (V,W) and advance Ct, Xt and Mt
            tic; 
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
            
            elseif col || ~isempty(Fparams.Tshell) 
                tic; 
                B = Nullsp.L*sigma-Lapp(Kernels.TD,sigma);
                timings.velocities.apply(i) = toc;
                
                % If inside shell, add traction from boundary correction to rhs
                if ~isempty(Fparams.Tshell)
                    B = B - Fparams.Tshell; 
                end
                
                % Solve Fredholm eq TD*mu = B
                tic; 
                mu = Lslv(Kernels.TD,B,parslv); 
                timings.velocities.solve(i) = toc;
                
                if ~col
                U   = Lapp(Kernels.SD,(mu+sigma)) + Ush; 
            
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
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Collision event
        if col
        tic; 
                
            fprintf('\n-------------------------------------------------');
            fprintf('\n Computing contact force at collision sites: \n')
            display(collist(:,1:2)')
            fprintf('-------------------------------------------------\n');
            
            % Compute contact force and force distribution updates
            [F_c,mu_c,rho_c] = LOCAL_Compute_Contact_LCP(collist,Kernels,Nullsp,Fparams,Ct,VW,dt);
            
            F_c = reshape(F_c,6,[]);
            display(F_c(1:3,:));
            
            if ~isempty(mu_c)
                %Update sigma, mu, U and VW
                mu   = mu    + mu_c; 
                sigma= sigma + rho_c; 
                U   = Lapp(Kernels.SD,(mu+sigma)); 
                
                if ~isempty(Fparams.Tshell)
                    U = U + Ush; 
                end
                
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
    
    function params = LOCAL_compute_shell_velocity(mu,params)
        shf = @(q) VshAna([q(1:3:end);q(2:3:end);q(3:3:end)],'VW'); 
        
        fprintf('\n Evaluation of flow on shell and boundary correction:')
        Xtrg = params.parsh.Xp; Nrtrg = params.parsh.Nrp; 
        Usqr = real(VSh_MatVec_RB_trg(mu,Xtrg,Nrtrg,params.parbd));
        Uh = shf(Usqr); 
        sigmah = -params.parsh.eigI.*Uh; 
        sigma  = real(VshSyn(sigmah,'VW'));
        sigma = reshape(reshape(sigma,[],3).',[],1);
        
        % Set params for shell density apply
        params.parsh.U = Usqr; 
        params.parsh.Uh = Uh; 
        params.parsh.shellden = sigma; 
        params.parsh.Vh = sigmah; 
        
        params.parsh.flag_pot = 'SL_Stk_3D'; params.parsh.a = 0; 
        
        Uinf = Usqr + VSh_MatVec_RB2(sigma,[],params.parsh);  
        fprintf('\n Check flow at boundary = 0: ||Uinf||_2 = %0.5g , ||Uinf||_inf = %0.5g',norm(Uinf),max(abs(Uinf))) 
        if max(abs(Uinf))>1e-3
        plot(log10(abs(shf(Uinf)))); pause(0.0001);   
        end
        
        end
        
        function Ctp = LOCAL_advance_center(Ct,dt,VW,Fparams)
        %global type; 
        
        if strcmp(Fparams.type,'Purcellxy')
            th = Fparams.tht + dt*VW(6,:);
            th = th - 2*pi*(floor(th./(2*pi))); 
            Fparams.tht = th; 
            Ctp = zeros(3,3); 
            Ctp(:,1:2) = [Rer(Fparams.tht(1));Rer(Fparams.tht(2));Rer(Fparams.tht(3))];  
        else
            Ctp = Ct + dt*VW(1:3,:)';
        end
    end
    
    function M = RotationMat(wh,t)
        nwh = norm(wh); 
        t = nwh*t; 
        wh = wh./nwh; 
        
        M = [
            1-(wh(2)^2+wh(3)^2)*(1-cos(t)) , wh(2)*wh(1)*(1-cos(t))-wh(3)*sin(t) , wh(1)*wh(3)*(1-cos(t))+wh(2)*sin(t);...
            wh(1)*wh(2)*(1-cos(t))+wh(3)*sin(t),1-(wh(1)^2+wh(3)^2)*(1-cos(t)),wh(2)*wh(3)*(1-cos(t))-wh(1)*sin(t);...
            wh(1)*wh(3)*(1-cos(t))-wh(2)*sin(t),wh(2)*wh(3)*(1-cos(t))+wh(1)*sin(t),1-(wh(2)^2+wh(1)^2)*(1-cos(t))
        ];
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
    
    function [colevent,collist,dt,Ctp] = LOCAL_collision_info(Fparams,X2,Ct,Ctp,VW,MRot,Mt,dt)
        %global diam rd mxrd sheps;  
        n3 = size(Ct,1); Shape=Fparams.parbd.Shape; 
        eps = Fparams.parbd.eps;  
        
        %Check for collision between spheres (or sphere envelopes)
        [colevent,collist,mindst,mindstsh] = LOCAL_check_collision_sph(Ctp,Fparams);
        
        % Finer collision detection for non-spheres
        if ~strcmp(Shape,'')
            np2 = size(X2,1)/n3;
            [colevent,collist,mindst,Xip,Xjp]=LOCAL_check_collision(collist,eps,Ctp,X2,MRot,Mt,VW,dt,np2); 
        end
        
        if colevent
            
            % If mindst<<eps, or <0, we need to adjust timestep
            bis=0; maxbis=3; shell = isfield(Fparams,'parsh');
            
            cond = mindst < 0.1*eps;
            if shell
                sheps = Fparams.parsh.eps;
                cond = cond || mindstsh < 0.1*sheps;
            end
            
            while cond && bis<=maxbis
            
                % Recompute dt and Ct{i+1} to avoid collision
                dt = dt/2;
                bis=bis+1; 
                Ctp = LOCAL_advance_center(Ct,dt,VW,Fparams); 
            
                %Check for collision between spheres (or sphere envelopes)
                [colevent,collist,mindst,mindstsh] = LOCAL_check_collision_sph(Ctp,Fparams);
        
                % Finer collision detection for non-spheres
                if ~strcmp(Shape,'')
                    [colevent,collist,mindst,Xip,Xjp]=LOCAL_check_collision(collist,eps,Ctp,X2,MRot,Mt,VW,dt,np2); 
                end
                
                fprintf('\n bisection = %d: dt = %1.4e, mindst = %1.4e, mindstsh = %1.4e',bis,dt,mindst,mindstsh);
                
                % update condition
                cond = mindst < 0.1*eps;
                if shell
                    sheps = Fparams.parsh.eps;
                    cond = cond || mindstsh < 0.1*sheps;
                end
            end
            
            fprintf('\n Min pairwise relative distance after collision detection: %2.4f ',mindst);
            if shell
                fprintf('\n Min distance to geometry after collision detection: %2.4f ',mindstsh);
            end 
            
        else
            collist=[]; 
            fprintf('\n Min pairwise relative distance: %2.4f',mindst);
            if isfield(Fparams,'parsh')
                fprintf('\n Min relative distance to geometry: %2.4f, absolute distance: %2.4f ',mindstsh,mindstsh*Fparams.parsh.rd);
            end
        end
    end
    
    function [colevent,collist,mindst,mindstsh] = LOCAL_check_collision_sph(C,Fparams)
        num_body = size(C,1); 
        
        rd = Fparams.parbd.rd; 
        mxrd = Fparams.parbd.mxrd; 
        diam = Fparams.parbd.diam; 
        eps = Fparams.parbd.eps;  
        
        if num_body > 1
            % Compute center distances
            distC = LOCAL_CenterDistance(C);
        
            % Find pairs for which (C_i-C-j) <= (r_i+r_j)+1.1*eps*max(r_i,r_j)
            [ii,jj]=meshgrid(1:num_body); 
            id = distC<=diam+1.1*eps*mxrd & ii<jj; 
            ip = ii(id); 
            jp = jj(id); 
        
            % Compute minimum relative distance between spheres
            mindst=min(reshape((distC-diam)./mxrd,[],1));
        else
            ip=[]; jp=[]; 
            mindst=Inf; 
        end
        
        % If there is a spherical shell, compute signed distance to boundary
        if isfield(Fparams,'parsh')
            rdsh = Fparams.parsh.rd; 
            sheps = Fparams.parsh.eps;
                
            NC = sqrt(sum(C.*C,2)); 
            distSh = rdsh - rd - NC; %(R-r) - ||C_i||
            mindstsh = min(distSh)/Fparams.parsh.rd; 
            iish = find(distSh <= 1.1*sheps*rdsh); 
            jjsh = (num_body+1)*ones(size(iish)); %j=n3+1 -> collision with boundary
            
            colevent = mindst < 1.1*eps || mindstsh < 1.1*sheps;
            collist = [[ip ; iish] [jp ; jjsh]];
        else
            colevent = mindst < 1.1*eps; 
            mindstsh = Inf;
            collist = [ip jp]; 
        end
    end
    
    function [colevent,collist,mindst,Xip,Xjp]=LOCAL_check_collision(collist,eps,C,X,MRot,Mt,VW,dt,np)
        
end

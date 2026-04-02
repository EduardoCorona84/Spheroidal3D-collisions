function Fparams = SpheroidalMS_initparams(Fparams)
%{
Helper code to initialize additional parameters for the mobility solver.
This is a wrapper function; see SpheroidalMS_set_params for more helpful
comments on the required inputs.
%}

disp('Initializing parameters...');

% Handling of variables for specific problem types
if isfield(Fparams, 'type')
switch Fparams.type
   case 'MHD'
      % MHD parameters for Laplace BIE
      Fparams.eta = (Fparams.mur-1)/(Fparams.mur+1); 
end
end

Fparams.background_flow = LOCAL_process_background_flow(Fparams);

% Params struct for rigid bodies
p = Fparams.parbd.p;
C = Fparams.parbd.Ct;
equ_radii = Fparams.parbd.equ_radii;
polar_radii = Fparams.parbd.polar_radii;
doAna=1; 
kerd=3;
    collision_eps = Fparams.parbd.collision_eps;
mdist=Fparams.parbd.mdist;
dense = Fparams.denseMV;
bodydist = Fparams.parbd.bodydist;
if isfield(Fparams.parbd, 'tsl_dealiasing')
    tsl_dealiasing = Fparams.parbd.tsl_dealiasing;
else
    tsl_dealiasing = true;
end
if isfield(Fparams.parbd, 'tsl_dealiasing_pad')
    tsl_dealiasing_pad = Fparams.parbd.tsl_dealiasing_pad;
else
    tsl_dealiasing_pad = 4;
end
if isfield(Fparams.parbd, 'tsl_backend') && ~isempty(Fparams.parbd.tsl_backend)
    tsl_backend = Fparams.parbd.tsl_backend;
else
    tsl_backend = 'spheroidal';
end

Fparams.parbd = SpheroidalMS_set_params(...
    equ_radii=equ_radii, ...
    polar_radii=polar_radii, ...
    p=p, ...
    C=C, ...
    collision_eps=collision_eps, ...
    mdist=mdist, ...
    doAna=doAna, ...
    flag_pot='SL_Stk_3D', ...
    kerd=kerd, ...
    dense=dense, ...
    bodydist=bodydist, ...
    tsl_backend=tsl_backend, ...
    tsl_dealiasing=tsl_dealiasing, ...
    tsl_dealiasing_pad=tsl_dealiasing_pad ...
);

% Time discretization (adams bashford multistep)
if strcmp(Fparams.tdisc,'abash')
   Fparams.abk = {1,[3/2 -1/2],[23/12 -4/3 5/12], [55/24 -59/24 37/24 -3/8]}; 
   Fparams.tdisc = 'abash'; 
   
   if Fparams.order > 2
      Fparams.tdisc2 = 'rk4';
   else
      Fparams.tdisc2 = 'trapz'; 
   end 
end

end %% END OF MAIN FUNCTION

function background_flow = LOCAL_process_background_flow(Fparams)
    background_flow = struct( ...
        'enabled', false, ...
        'U0', zeros(1,3), ...
        'A', zeros(3,3) ...
    );

    if ~isfield(Fparams, 'background_flow') || isempty(Fparams.background_flow)
        return;
    end

    background_flow.enabled = Fparams.background_flow.enabled;
    background_flow.U0 = Fparams.background_flow.U0;
    background_flow.A = Fparams.background_flow.A;

    trace_tol = 1e-12;
    assert(abs(trace(background_flow.A)) <= trace_tol, ...
        'Fparams.background_flow.A must have trace 0 to represent incompressible flow.');
end

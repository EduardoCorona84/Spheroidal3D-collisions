function Fparams = SpheroidalMS_initparams(Fparams)
%{
Helper code to initialize additional parameters for the mobility solver.
This is a wrapper function; see SpheroidalMS_set_params for more helpful
comments on the required inputs.
%}

disp('Initializing parameters...');

% Params struct for rigid bodies
p = Fparams.parbd.p;
C = Fparams.parbd.Ct;
equ_radii = Fparams.parbd.equ_radii;
polar_radii = Fparams.parbd.polar_radii;
doAna=1; 
kerd=3;
eps = Fparams.parbd.eps;
mdist=Fparams.parbd.mdist;
out=Fparams.parbd.out;
dense = Fparams.denseMV;
bodydist = Fparams.parbd.bodydist;

Fparams.parbd = SpheroidalMS_set_params(...
    equ_radii=equ_radii, ...
    polar_radii=polar_radii, ...
    p=p, ...
    C=C, ...
    eps=eps, ...
    mdist=mdist, ...
    out=out, ...
    doAna=doAna, ...
    flag_pot='SL_Stk_3D', ...
    kerd=kerd, ...
    dense=dense, ...
    bodydist=bodydist ...
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
    
    
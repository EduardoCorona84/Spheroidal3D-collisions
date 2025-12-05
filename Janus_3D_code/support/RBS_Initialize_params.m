function Fparams = RBS_Initialize_params(Fparams)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extra parameters depending on type
Fparams.eta=1; 
switch Fparams.type
    case 'MHD'
        % MHD parameters for Laplace BIE
        Fparams.eta = (Fparams.mur-1)/(Fparams.mur+1); 
    case 'Purcellxy'
        fprintf('\n Generating C from initial angles (Purcell)\n'); 
        Fparams.parbd.Ct = zeros(3,3); 
        Rer = @(th) Fparams.R*[cos(th) sin(th)]; 
        Fparams.parbd.Ct(:,1:2) = [Rer(Fparams.tht(1));Rer(Fparams.tht(2));Rer(Fparams.tht(3))];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% params struct for rigid bodies
p = Fparams.parbd.p; C = Fparams.parbd.Ct; rd = Fparams.parbd.rd; doAna=1; 
kerd=3; mdist=Fparams.parbd.mdist; out=Fparams.parbd.out; Shape=Fparams.parbd.Shape;  
% Model Surface
Sc = cell(p,1); 
Sc{p} = SurfaceSph(shape_gallery(p,Shape));

Fparams.parbd = RBS_set_params(p,C,rd,'SL_Stk_3D',Shape,Sc,kerd,Fparams.denseMV,doAna,mdist,Fparams.parbd.eps,out); 
if isfield(Fparams, 'lofi')
    lofi_p = Fparams.lofi.p;
    lofi_Sc = cell(lofi_p,1); 
    lofi_Sc{lofi_p} = SurfaceSph(shape_gallery(lofi_p,Shape));
    Fparams.lofi = RBS_set_params(lofi_p,C,rd,'SL_Stk_3D',Shape,lofi_Sc,kerd,Fparams.denseMV,doAna,mdist,Fparams.lofi.eps,out); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize shell boundary (if it exists)
if isfield(Fparams,'parsh')
   %Eigenvalues for Inverse Operator (at outer shell) 
   %interior flow
   outsh=0; psh = Fparams.parsh.psh; shrd = Fparams.parsh.shrd; 
   mdsh = Fparams.parsh.mdist; kerd=3; 
   [SV,SW,SX] = RBS_get_eigen('SMat',outsh); %VSh eigenvalues for Single Layer (first kind eq)
   spsh=(psh+1)^2; 
   ii = (1:spsh)'; 
   nn = floor(sqrt(ii-1));
   eigS = [SV(nn);SW(nn);SX(nn)]; 
   eigIS = (1/shrd)./(eigS+(eigS==0)) - (1/shrd).*(eigS==0); 
   
   % Update parameters for kernel eval 
   Scsh{psh,1} = SurfaceSph(shrd*shape_gallery(psh,'')); doAn = 0; outsh=0; 
   Fparams.parsh = RBS_set_params(psh,[0 0 0],shrd,'TSL_Stk_3D','',Scsh,kerd,Fparams.denseMV,doAn,mdsh,Fparams.parsh.eps,outsh);  

   %Initial outer flow is zero (pass empty density, potential and eigenvalues)
   Fparams.parsh.shellden = []; Fparams.parsh.Vh = []; Fparams.parsh.eigI = eigIS;  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rotation-based singular quad if non-spherical geometry 
if ~strcmp(Fparams.parbd.Shape,'')
   Fparams.typeMV='Rbs';  
end

% time discretization (adams bashford multistep)
if strcmp(Fparams.tdisc,'abash')
   Fparams.abk = {1,[3/2 -1/2],[23/12 -4/3 5/12], [55/24 -59/24 37/24 -3/8]}; 
   Fparams.tdisc = 'abash'; 
   
   if Fparams.order>2
      Fparams.tdisc2='rk4';
   else
      Fparams.tdisc2='trapz'; 
   end 
end

if ~ isfield(Fparams, 'loadIntermediate')
    Fparams.loadIntermediate = false;
end

if ~ isfield(Fparams, 'plotFlag')
    Fparams.plotFlag = false;
end
end


function [SPx,SPy,SPz]=spheroidalMatVec_gradSP(params,X)
%--------------------------------------------------------------------%
% spheroidalMatVec for cartesian gradient of S' operator
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
%    (o) X = target points at which to evaluate the potential. If not
%        provided, this function evaluates the potential at points on each
%        spheroid surface.
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

nu_x=repmat([1,0,0],size(X,1),1,ns);
nu_y=repmat([0,1,0],size(X,1),1,ns);
nu_z=repmat([0,0,1],size(X,1),1,ns);

% start gradient calculation
nt=size(X,1);

[sep,Xt]=params.separate_targets(X);
X_spectral = cell(1,ns);
nu_x_spectral= cell(1,ns); nu_y_spectral=nu_x_spectral; nu_z_spectral=nu_x_spectral;

for i=1:ns
    X_spectral{i} = Xt(sep(:,i)==0,:,i);
    nu_x_spectral{i} = nu_x(sep(:,i)==0,:,i);
    nu_y_spectral{i} = nu_y(sep(:,i)==0,:,i);
    nu_z_spectral{i} = nu_z(sep(:,i)==0,:,i);
end

% Do "nearby" evaluation with spheroidal harmonics
[LP_spectral_cell_x,LP_spectral_cell_y,LP_spectral_cell_z] = spheroidalSP(params,X_spectral,nu_x_spectral,nu_y_spectral,nu_z_spectral);

% Do "far" evaluation with smooth quadrature

SPx=zeros(nt,nf); SPy=SPx; SPz=SPx;

% Calculate effect from each particle on each target point
for i=1:ns

    % Add contributions from spectral method
    SPx(sep(:,i)==0,:) = SPx(sep(:,i)==0,:) + LP_spectral_cell_x{i}; 
    SPy(sep(:,i)==0,:) = SPy(sep(:,i)==0,:) + LP_spectral_cell_y{i};
    SPz(sep(:,i)==0,:) = SPz(sep(:,i)==0,:) + LP_spectral_cell_z{i}; 

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
        
        pot='dSL_L_3D';
        KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-8,2,400,1);
        KEparams.dim = 3;
        [~, gwt]=g_grid(p+1);
        wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
        wt = wt(:);
        Wns = Sns.geoProp.W; Wns= Wns.*wt;
        KEparams.X = Xself;
        KEparams.W2 = Wns.';

        Nrns=nu_x(sep(:,i)==1,:,i);
        KEparams.nor = Nrns; 
        LP_smooth = Kernel_Eval(X_smooth,Xself,KEparams)*sigma(:,:,i);
        SPx(sep(:,i)==1,:) = SPx(sep(:,i)==1,:) + LP_smooth;

        Nrns=nu_y(sep(:,i)==1,:,i);
        KEparams.nor = Nrns; 
        LP_smooth = Kernel_Eval(X_smooth,Xself,KEparams)*sigma(:,:,i);
        SPy(sep(:,i)==1,:) = SPy(sep(:,i)==1,:) + LP_smooth;

        Nrns=nu_z(sep(:,i)==1,:,i);
        KEparams.nor = Nrns; 
        LP_smooth = Kernel_Eval(X_smooth,Xself,KEparams)*sigma(:,:,i);
        SPz(sep(:,i)==1,:) = SPz(sep(:,i)==1,:) + LP_smooth;

    end 
end

if isReal
    SPx = real(SPx);
    SPy = real(SPy);
    SPz = real(SPz);
end

end

clear;
close all
parr=[8];
errs=zeros(size(parr));
obl=0;
is_flux=1;
use_kernel=0;

for l=1:length(parr)

p=parr(l);
np=2*p*(p+1);
% u0=2/sqrt(3);
u0=1.1;
if ~obl
    a=1/u0;
else
    a=1/sqrt(u0^2+1);
end
params=SpheroidalParameters;
params.p=p;
params.u0=u0; params.a=a; params.isReal=0; params.oblate=obl;
params.centers=zeros(1,3);
params.matvec_eta=10;
[theta,phi]=gl_grid(p);
v=cos(theta);

% if ~obl
%     Shape='ellipseZ'; 
% else
%     Shape='oblateZ';
% end
% Sns = SurfaceSph(shape_gallery(p,Shape)); 
if ~obl
    Xself=prolate_spheroid_shape(p,u0,a);
else
    Xself=oblate_spheroid_shape(p,u0,a);
end
Sns=SurfaceSph(Xself);

% First: Create point charge. Add up point charge strength times green
% function to get potential that we want to solve for on the surface.
% Either do this in cartesian and then convert, or do it all in spheroidal
% coordinates.

Xptch=[[0,0.02,0];[.15,.1,.1]];
% Xptch=[0,0,0];
ptch=[1.5,2]';
% ptch=[1];

Ys=[u0*ones(length(v),1) , v, phi];
Y=spheroidal2cart(Ys,a,obl);
NrY=params.get_Norm(p);
% quiver3(Xself(:,1),Xself(:,2),Xself(:,3),NrY(:,1),NrY(:,2),NrY(:,3)); hold on;
% plot3(Xptch(:,1),Xptch(:,2),Xptch(:,3),'r*'); 
% quiver3(Y(:,1),Y(:,2),Y(:,3),NrY(:,1),NrY(:,2),NrY(:,3));hold off;

% We will use this for surface boundary conditions
if ~is_flux
    if use_kernel
        truesolnSurf=KernelPot(ptch,Xptch,Y);
    else
        truesolnSurf=PtChargePotential(ptch,Xptch,Y); 
    end
else
    if use_kernel
        truesolnSurf=KernelFlux(ptch,Xptch, Y,NrY);
    else
        truesolnSurf=PtChargeFlux(ptch,Xptch,Y,NrY);
    end
end


% Construct DL/SL on-surface matrices using KernelDLap singular quadrature
%----------------------------------------------------------%
% Laplace SL, S' and DL on the spheroid surface
% [SM,Sp,DM] = kernelDLap(Sns);
Sp=spheroidalMatVecKernel(params,'SP',p);

%----------------------------------------------------------%
[~, gwt]=g_grid(p+1);
wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
wt = wt(:)';
W = Sns.geoProp.W;
RY=sqrt(sum(Y.^2,2));
if ~is_flux
    C=(1./RY)*(W'.*wt); % Completion Term 1/||x|| * integral of sigma over surface
    % Final operator: 1/2 I + D + C
    K=.5*eye(np)+DM+C;
else
    C = ones(np,1)*(W'.*wt);
    % Final operator: -1/2 I + SP + C
    K=-.5*eye(np) + Sp;

    % %% test C
    % sint=integrateOverS(Sns,sum(ptch));
    % smatvec=C*repmat(sum(ptch),np,1);
end
%----------------------------------------------------------%

sigma=K\truesolnSurf;
% sigma = reshape(sigma_vec,np,[],ns);
params.sigma=sigma;
params.get_shc;     

% Check if they match 

%Test points: scale u coordinate up a little bit so we're just off the
%surface.
if ~obl
    normal=1./sqrt((u0^2-1).*(u0^2-v.^2)).*[u0.*sqrt(u0.^2-1).*sqrt(1-v.^2).*cos(phi) u0.*sqrt(u0.^2-1).*sqrt(1-v.^2).*sin(phi) (u0.^2-1).*v];
else
    normal=1./sqrt((u0^2+1).*(u0^2+v.^2)).*[u0.*sqrt(u0.^2+1).*sqrt(1-v.^2).*cos(phi) u0.*sqrt(u0.^2+1).*sqrt(1-v.^2).*sin(phi) (u0.^2+1).*v];
end
% Xeval=spheroidal2cart(Ys,a,obl);
% Xeval=Xeval+1.1.*normal;
Xeval1=Y+0.8.*normal;
Xeval=params.set_X_targets({Xeval1});
RXeval=sqrt(sum(Xeval.^2,2));

% Solution from solving integral equation
sigmaSurfInt=integrateOverS(Sns,sigma); 

if ~is_flux
    soln=spheroidalDL(params,Xeval)+sigmaSurfInt./RXeval;
else
    soln=spheroidalSL(params,Xeval);
end

% Actual solution
truesoln=PtChargePotential(ptch,Xptch,Xeval);

errs(l)=max(abs(truesoln-soln));

fprintf('p=%d: infinity-norm err=%.6e\n',p, max(abs(truesoln-soln)))

% close all
figure;
semilogy(1:np,abs(truesoln-soln)/abs(truesoln))
xlabel('evaluation point');
ylabel('relative error')

end




function pcp = PtChargePotential(ptch,Xptch,Y)
M=length(ptch);
np=length(Y);
Rptch=zeros(np,M);
for i=1:M
    Rptch(:,i)=sqrt(sum((Xptch(i,:)-Y).^2,2));
end
pcp=1./(4*pi*Rptch)*ptch; 
end

function flux=PtChargeFlux(ptch,Xptch,Y,NrY)
    M=length(ptch);
    np=length(Y);
    Rptch=zeros(np,M);
    EdotN=zeros(np,M);
    for i=1:M
        r=Xptch(i,:)-Y;
        Rptch(:,i)=sqrt(sum(r.^2,2));
        EdotN(:,i)=dot(NrY,r./Rptch(:,i),2);
    end
    flux=EdotN./(4*pi*Rptch.^2)*ptch;
end

function flux=KernelPot(ptch, Xptch, Y)
    pot='SL_L_3D';
    KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-8,2,400,1);
    KEparams.dim = 3;
    KEparams.X = Y;
    KEparams.W2 = ones(size(ptch))';
    SLmat = Kernel_Eval(Y,Xptch,KEparams);
    flux=SLmat*ptch;
end

function flux=KernelFlux(ptch, Xptch, Y, NrY)
    pot='dSL_L_3D';
    KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-8,2,400,1);
    KEparams.dim = 3;
    KEparams.X = Y;
    KEparams.W2 = ones(size(ptch))';
    KEparams.nor = NrY; 
    SPmat = Kernel_Eval(Y,Xptch,KEparams);
    flux=SPmat*ptch;
end


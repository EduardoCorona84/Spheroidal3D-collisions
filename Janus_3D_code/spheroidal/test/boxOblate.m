% Test: loop over u-values getting closer to spheroid surface.
%         (o) Solve point charge potential BIE to get surface density
%         (o) Evaluate DL, SL of density with smooth quadrature
%         (o) Evaluate DL, SL of density with spheroidal harmonics
%         (o) Compare each against true point charge potential
% Calculate mean error, put box plots on a plot over -log10(distance to surface)
% 1,2,3, and 4. Do for multiple p values.

parr=[4,8,12,18,24];  %p-values to test
distance=10.^(-(1:6)); %Distance from spheroid surface
errs=cell(length(parr),1);

u0=2/sqrt(3);
a=1/sqrt(1+u0^2);
params=SpheroidalParameters; 
params.u0=u0; params.a=a; params.isReal=true; 
params.oblate=1;

for l=1:length(parr)

    p=parr(l);
    np=2*p*(p+1);
    [theta,phi]=gl_grid(p);
    v=cos(theta);

    this_p_errs=zeros(np,length(distance));

    sph_norm=[1,0,0];
    Xself=oblate_spheroid_shape(p,u0,a);
    Sself=cart2spheroidal(Xself,a,1);
    Sns=SurfaceSph(Xself);
    [normal,~]=spheroidalNu2cart(sph_norm,Sself,a,1);
    % could have used line below.
    normal_expl=1./sqrt((u0^2+1).*(u0^2+v.^2)).*[u0.*sqrt(u0.^2+1).*sqrt(1-v.^2).*cos(phi) u0.*sqrt(u0.^2+1).*sqrt(1-v.^2).*sin(phi) (u0.^2+1).*v];
    % display(norm(abs(normal_expl-normal)));
    
    % First: Create point charge. Add up point charge strength times green
    % function to get potential that we want to solve for on the surface.
    % Either do this in cartesian and then convert, or do it all in spheroidal
    % coordinates.
    
    Xptch=[[0,0,0];[.1,.1,.1];[.2,-.2,.1]];
    % Xptch=[0,0,0];
    ptch=[1,2,-.5]';
    
    Ys=[u0*ones(length(v),1) , v, phi];
    Y=spheroidal2cart(Ys,a,1);
    
    % We will use this for surface boundary conditions
    truesolnSurf=PtChargePotential(ptch,Xptch,Y); 
    
    % Construct DL/SL on-surface matrices using KernelDLap singular quadrature
    %----------------------------------------------------------%
    % Laplace SL, S' and DL on the spheroid surface
    [SM,Sp,DM] = kernelDLap(Sns);
    %----------------------------------------------------------%
    
    % Completion Term 1/||x|| * integral of sigma over surface
    %----------------------------------------------------------%
    [~, gwt]=g_grid(p+1);
    wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
    wt = wt(:)';
    W = Sns.geoProp.W;
    RY=sqrt(sum(Y.^2,2));
    C=(1./RY)*(W'.*wt); %Completion term
    %----------------------------------------------------------%
    
    % Final operator: -1/2 I + D + C
    K=.5*eye(np)+DM+C;
    % K=-.5*eye(np)+DM;
    sigma=K\truesolnSurf;
    params.sigma=sigma;
    
    for j=1:length(distance)

        % Check if they match 
        
        %Test points: scale u coordinate up a little bit so we're just off the
        %surface.
        Xeval=spheroidal2cart(Ys,a,1);
        Xeval=Xeval+distance(j).*normal;
        RXeval=vecnorm(Xeval,2);
        % TODO: bit different from NearSingularTest
        
        % Solution from solving integral equation
        sigmaSurfInt=integrateOverS(Sns,sigma); 
        
        soln=spheroidalDL(params,Xeval)+sigmaSurfInt./RXeval;
        % soln=spheroidalDL(params,Xeval,'cart');
        % TODO: size here is np x 3 whereas in NearSingTest is np x 1.
        
        % Actual solution
        truesoln=PtChargePotential(ptch,Xptch,Xeval);
        
        this_p_errs(:,j)=log10(vecnorm(abs(truesoln-soln)')'./max(abs(truesoln)));
        
%         close all
%         figure;
%         semilogy(1:np,abs(truesoln-soln)/abs(truesoln))
%         xlabel('evaluation point');
%         ylabel('relative error')
        
    end

    errs{l}=this_p_errs;

end


close all
figure;
hold on;
colors=['k','b','r','g','y'];
for l=1:length(parr)
    boxchart(errs{l});
end
legend("p="+string(parr));
xlabel("-log10(distance)")
ylabel("log10 relative error")
hold off; 
pause;

%%
% close all

% % testing normal vectors
% p=18;
% [theta,phi]=gl_grid(p);
% v=cos(theta);
% normal=1./sqrt((u0^2+1).*(u0^2+v.^2)).*[u0.*sqrt(u0.^2+1).*sqrt(1-v.^2).*cos(phi) u0.*sqrt(u0.^2+1).*sqrt(1-v.^2).*sin(phi) (u0.^2+1).*v];
% % Plot x-z slice of spheroid normals for theta=0 and theta=pi
% ind=[1:p+1 p*(p+1)+1:(p+1)^2]';
% x1=normal(ind,1); 
% z1=normal(ind,3);
% [X1,Z1]=meshgrid(x1,z1);
% Ss=[u0*ones(length(v),1) , v, phi];
% S=spheroidal2cart(Ss,a);
% x=S(ind,1);
% z=S(ind,3);
% [X,Z]=meshgrid(x,z);
% 
% figure;
% hold on
% axis equal
% scatter(x,z)
% quiver(x,z,x1,z1)
% hold off

function pcp = PtChargePotential(ptch,Xptch,Y)
M=length(ptch);
np=length(Y);
Rptch=zeros(np,M);
for i=1:M
    Rptch(:,i)=sqrt(sum((Xptch(i,:)-Y).^2,2));
end
pcp=1./(4*pi*Rptch)*ptch; 
end


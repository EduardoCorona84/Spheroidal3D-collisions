




% Coming up with heuristic for spheroid distance.

pars=SpheroidalParameters;
pars.sigma=[func(1,1,16),func(1,0,16)];
pars.u0=[1.1 1.2];
pars.a= 1./pars.u0;
pars.centers = [0 0 0; 1.5 0 0];
pars.thetas=[0 pi*4/5];
pars.phis=[0 pi*11/7];

[Xt,Xs]=pars.get_X(1);
X=[Xt; Xs];

np=size(Xt,1);

distances=[];
for i=1:np
    Xi=repmat(Xs(i,:),np,1);
    distances = [distances; vecnorm(Xi - Xt, 2,2)];
end

closest= min(distances);
closest

%  Now heuristic

vd1= Xs(1,:) - Xs(pars.p+1,:);
vd2= Xt(1,:) - Xt(pars.p+1,:);

c2=pars.centers(2,:);

alignment1 = abs((vd1/norm(vd1))*(c2/norm(c2))');
alignment2 = abs((vd2/norm(vd2))*(c2/norm(c2))');
sphere_est1 = pars.a(1) * sqrt(pars.u0(1)^2-1+alignment1);
sphere_est2 = pars.a(2) * sqrt(pars.u0(2)^2-1+alignment2);

sphere_est1A = pars.a(1) * sqrt(pars.u0(1)^2-1) + (pars.a(1)*pars.u0(1) - pars.a(1) * sqrt(pars.u0(1)^2-1))*alignment1;
sphere_est2A = pars.a(2) * sqrt(pars.u0(2)^2-1) + (pars.a(2)*pars.u0(2) - pars.a(2) * sqrt(pars.u0(2)^2-1))*alignment2;

closest_est = norm(c2) - sphere_est1-sphere_est2;
closest_est

closest_estA = norm(c2) - sphere_est1A-sphere_est2A;
closest_estA %This one seems to be working better


% Both of these heuristics err on the side of caution, will categorize
% spheroids as too close even when they aren't.

theta=pars.thetas(2); phi=pars.phis(2);
vd=[sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)];

close all
figure;
scatter3(X(:,1),X(:,2),X(:,3));
xlabel('x')
ylabel('y')
zlabel('z')
hold on
quiver3(c2(1),c2(2),c2(3),vd(1),vd(2),vd(3))
legend('spheroids','vd')
hold off
axis equal


% -----------------------------
% Procedure to separate near and far given pars, spheroid number i

pars=SpheroidalParameters;
pars.sigma=[func(1,1,16),func(1,0,16),func(1,-1,16)];
pars.u0=[1.1 1.2 1.3];
pars.a= 1./pars.u0;
pars.centers = [0 0 0; 1.5 0 0; 2 2 2];
pars.thetas=[0 0 pi/3];
pars.phis=[0 0 pi/10];
i=1;
[Xt,Xs]=pars.get_X(i);
X=[Xt; Xs];

[np,ns]=size(pars.sigma);
ind=[1:i-1 i+1:ns];
at = pars.a(ind); %target a's
u0t = pars.u0(ind); %target u0's
ct=pars.centers(ind,:);
thetat=pars.thetas(ind);
phit=pars.phis(ind);
ai=pars.a(i); %self a
u0i=pars.u0(i); %self a
ci=pars.centers(i,:);
thetai=pars.thetas(i);
phii=pars.phis(i);

vdt=[sin(thetat').*cos(phit') sin(thetat').*sin(phit') cos(thetat')];
vdt=vdt./repmat(vecnorm(vdt,2,2),1,3);
% cdt=ct./repmat(vecnorm(cdt,2,2),1,3);
vdi=[sin(thetai).*cos(phii) sin(thetai).*sin(phii) cos(thetai)];
vdi=vdi./norm(vdi);
% cdi=ci/norm(ci);

displacements = ct-repmat(ci,ns-1,1);
d=vecnorm(displacements,2,2);
cd = displacements./repmat(d,1,3); 

alignment = abs(sum(vdt.*cd,2))';
self_alignment = abs(vdi*cd');
target_radii = at.*sqrt(u0t.^2-1) + (at.*u0t - at.*sqrt(u0t.^2-1)).*alignment;
% Radii of approximate spheres around centered ("self") spheroid
self_radii = ai.*sqrt(u0i.^2-1) + (ai.*u0i - ai.*sqrt(u0i.^2-1)).*self_alignment;

closest_est = d' - self_radii - target_radii;

% separation = -eye(ns); %-1=self, 0 =close, 1= far
% separation(i:end,:) = 

% True distance
distances=zeros(1,ns-1);
xs=Xs(:,1); ys=Xs(:,2); zs=Xs(:,3);
for k=1:ns-1
    xt=Xt((k-1)*np+1:k*np,1); 
    yt=Xt((k-1)*np+1:k*np,2); 
    zt=Xt((k-1)*np+1:k*np,3);

    D=sqrt( (xs-xt').^2+(ys-yt').^2+(zs-zt').^2) ;
    distances(k)=min(D(:));
end

closest_est
distances

% Plot 1: particle i vs first target
[thplt,phplt]=gl_grid(16);
self_sphere_points1 = [self_radii(1).*sin(thplt).*cos(phplt) self_radii(1).*sin(thplt).*sin(phplt) self_radii(1).*cos(thplt)];
self_sphere_points2 = [self_radii(2).*sin(thplt).*cos(phplt) self_radii(2).*sin(thplt).*sin(phplt) self_radii(2).*cos(thplt)];
target_sphere_points1 = repmat(ct(1,:),16*2*17,1) + [target_radii(1).*sin(thplt).*cos(phplt) target_radii(1).*sin(thplt).*sin(phplt) target_radii(1).*cos(thplt)];
target_sphere_points2 = repmat(ct(2,:),16*2*17,1) + [target_radii(2).*sin(thplt).*cos(phplt) target_radii(2).*sin(thplt).*sin(phplt) target_radii(2).*cos(thplt)];

X1=[X; Xt(1:2*pars.p*(pars.p+1),:)];
spheres1=[self_sphere_points1; target_sphere_points1];
np_spheroids=size(X1,1);
np_spheres = size(spheres1,1);
X2=[X; Xt(2*pars.p*(pars.p+1)+1:end,:)];
spheres2=[self_sphere_points2; target_sphere_points2];

points1=[X1;spheres1];
points2=[X2;spheres2];
colors=[repmat([1 0 0],np_spheroids,1); repmat([0 0 1],np_spheres,1)];

figure;
scatter3(points1(:,1),points1(:,2),points1(:,3),36,colors);
xlabel('x')
ylabel('y')
zlabel('z')
axis equal

figure;
scatter3(points2(:,1),points2(:,2),points2(:,3),36,colors);
xlabel('x')
ylabel('y')
zlabel('z')
axis equal

% -----------------------------

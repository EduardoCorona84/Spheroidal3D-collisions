u0=2/sqrt(3); a=1/u0; p=8;

pars=SpheroidalParameters;
pars.sigma=[func(1,1,p) func(1,0,p)];
pars.u0 = [u0 u0];
pars.a = [a a];
pars.centers = [0 -1 0 ; 0 3 0];
pars.thetas = [0  pi/2];
pars.phis = [0 0];

pars


[Xt,Xs]=pars.get_X(1);
X=[Xt; Xs];

close all
figure;
scatter3(X(:,1),X(:,2),X(:,3));
axis equal

[Xt,Xs]=pars.get_X(2);
X=[Xt; Xs];

figure;
scatter3(X(:,1),X(:,2),X(:,3));
axis equal


pars.sigma=[func(1,1,p)];
pars.u0 = u0;
pars.a = a ;
pars.centers = [0 0 0];
pars.thetas = 0;
pars.phis = 0;

X = pars.get_X;
[x,y]=meshgrid(X(:,1),X(:,2));
z=repmat(X(:,3),1,length(X(:,3)));
figure;
surf(x,y,z)
axis equal

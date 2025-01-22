clear;clc;
hold off

addpath ./support;
C=[-1.5,1,0; -3,0,0; -1,0,2];
r=[0.5; 1; 0.7];
lambda=0.1;
compute = 0;

if (compute ==1)
[X,L,x1]=BIEtest_func(C,r,lambda);
Size=size(L,1);
Pvals=[4 8 12 16 24];

save 'BIE_TEST_DATA'
cd ..
else
     load('BIE_TEST_DATA.mat');
end

% plot spheres 
hold on
sphere_nodes = size(L,1) / size(C,1);
L=smooth3(L);
errormap = L(:,:,4);
dotsize=[100,125,100];
for k=1:3
    sz=dotsize(k);
    radius=r(k);
    Cen = C(k,:);
    Interval = sphere_nodes*(k-1) +[1:sphere_nodes];
    error=errormap(Interval,3);
    Colormap=meshgrid(error,error);
    [X, Y, Z]=sphere(sphere_nodes - 1);
    X=radius*X + Cen(1); Y=radius*Y + Cen(2); Z=radius*Z + Cen(3);
    scatter3(x1(Interval,1),x1(Interval,2),x1(Interval,3),repmat(sz,size(Interval),1),error,'filled');
    
end
hcb=colorbar;
hcb.Label.String = 'evaluation error (log10) ';
hcb.FontSize = 15;
title('Error in BIE Test','FontSize',15);
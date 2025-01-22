clear;clc;
load('BIE_TEST_DATA');


figure
hold on
for i=1:size(r)
    
    [Sx, Sy, Sz]=sphere(40);
    Sx=r(i)*Sx;
    Sy=r(i)*Sy;
    Sz=r(i)*Sz;
    Sc=C(i,:);
    mesh(Sx+Sc(1),Sy+Sc(2),Sz+Sc(3))%,'edgecolor', 'k');
    
end
view(50,10);
ax=gca; ax.FontSize = 16;

% clear;clc;
% load('BIE_TEST_DATA');
set(0,'DefaultFigureWindowStyle','docked')
figure
ax = axes;
ax.FontSize = 16;
grid on
view(50,1)
axis vis3d
% axis tight manual

grp = hgtransform(Parent=ax);
n3 = size(Ct{1},1);
nt = numel(Ct);
if ~exist('r','var')
    r = ones(n3,1);
end

for i = 1:nt
    C = Ct{i};
    for k=1:n3
        [Sx, Sy, Sz]=sphere(40);
        Sx=r(k)*Sx;
        Sy=r(k)*Sy;
        Sz=r(k)*Sz;
        Sc=C(k,:);
        surf(Sx+Sc(1),Sy+Sc(2),Sz+Sc(3),Parent=grp,EdgeColor='k');
    end
    drawnow
end


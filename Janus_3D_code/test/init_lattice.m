function [C, init_dir] = init_lattice(n, Cdst) 
lx=0:Cdst:Cdst*(n-1);
lx = lx - mean(lx); 
[xx,yy,zz] = meshgrid(lx);
C = [xx(:) yy(:) zz(:)]; 
% randomize centers and/or radii
try %#ok<TRYNC>
    rng('default'); 
end
meanSep = Cdst - 2;
while true
    C = C + meanSep/4*rand(size(C));
    minDist = computePairwiseDistance(C);
    if all(minDist > meanSep/4 )
        break
    end
end
% initial particle orientations
%% NIC: set the initial direction to be towards the center 
init_dir= -C;
for i = 1:size(C,1)
    init_dir(i,:) = init_dir(i,:) / norm(init_dir(i,:));
end

function minDist = computePairwiseDistance(C)
n = size(C,1);
[Y_g1,  X_g1  ] = meshgrid(C(:,1), C(:,1));
[Y_g2,  X_g2  ] = meshgrid(C(:,2), C(:,2));
[Y_g3,  X_g3  ] = meshgrid(C(:,3), C(:,3));
d1 = (X_g1 - Y_g1); d2 = (X_g2 - Y_g2); d3 = (X_g3 - Y_g3); 
pairwiseDistance = sqrt(d1.^2 + d2.^2 + d3.^2);
minDist = zeros(n-1,1);
for i =1:n-1
    minDist(i) = min(pairwiseDistance(i,i+1:end));
end
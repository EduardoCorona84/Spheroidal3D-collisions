function [C, rd, init_dir] = init_lattice(n, Cdst, rd, polydisperseRatio) 
lx=0:Cdst:Cdst*(n-1);
lx = lx - mean(lx); 
[xx,yy,zz] = meshgrid(lx);
C = [xx(:) yy(:) zz(:)]; 
n3 = size(C,1);
rd=rd*ones(n3,1)+polydisperseRatio*randn(n3,1);
% randomize centers and/or radii
try %#ok<TRYNC>
    rng('default'); 
end
eps = .1;
while true
    C = C + eps*randn(size(C));
    minDist = computePairwiseDistance(C,rd);
    if all(minDist > .1 )
        break
    end
end
% initial particle orientations
%% NIC: set the initial direction to be towards the center 
init_dir = -C;
for i = 1:size(C,1)
    init_dir(i,:) = init_dir(i,:) / norm(init_dir(i,:));
end

function minDist = computePairwiseDistance(C,rd)
n = size(C,1);
minDist = zeros(n-1,1);
for i =1:n-1
    pairwiseDistance = arrayfun(@(j) norm(C(i,:) - C(j,:)) - rd(i) - rd(j), i+1:n);
    minDist(i) = min(pairwiseDistance);
end
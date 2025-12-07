function [C, rd, init_dir] = init_lattice(n, Cdst, meanRadius, polydisperseRatio) 
while true
    lx=0:Cdst:Cdst*(n-1);
    lx = lx - mean(lx); 
    [xx,yy,zz] = meshgrid(lx);
    C = [xx(:) yy(:) zz(:)]; 
    n3 = size(C,1);
    rd=meanRadius*(1+polydisperseRatio*rand(n3,1));
    delta = meanRadius*2*eps*(.5 - rand(n3,1));
    C = C + delta;
    minDist = computePairwiseDistance(C,rd);
    if all(minDist > 0.1 )
        break
    end
    Cdst = Cdst*1.01;
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
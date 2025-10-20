function [C, init_dir] = init_lattice(n, Cdst) 
lx=0:Cdst:Cdst*(n-1);
lx = lx - mean(lx); 
[xx,yy,zz] = meshgrid(lx);
C = [xx(:) yy(:) zz(:)]; 
% randomize centers and/or radii
try %#ok<TRYNC>
    rng('default'); 
end
meanSep = cDist - 2;
while true
    C = C + meanSep/4*rand(size(C));
    [~,~,mindst,~] = LOCAL_check_collision_sph(C,Fparams);
    if all(mindst > ep/10 )
        break
    end
end
% initial particle orientations
%% NIC: set the initial direction to be towards the center 
init_dir= -C;
for i = 1:size(C,1)
    init_dir(i,:) = init_dir(i,:) / norm(init_dir(i,:));
end

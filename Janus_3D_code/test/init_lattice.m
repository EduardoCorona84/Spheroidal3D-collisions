function [C, rd, init_dir] = init_lattice(n, Cdst, meanRadius, ...
    polydisperseRatio, colThresh, desiredNumCol, desiredTol, debug) 
if ~exist('colThresh','var') || isempty(colThresh)
    colThresh = 0.3;
end
if ~exist('desiredNumCol','var') || isempty(desiredNumCol)
    desiredNumCol = floor(n^3 / 2);
end
if ~exist('desiredTol','var') || isempty(desiredTol)
    desiredTol = floor(n^3 / 10);
end
if ~exist('debug','var') || isempty(debug)
    debug = false;
end
lx=0:Cdst:Cdst*(n-1);
lx = lx - mean(lx); 
[xx,yy,zz] = meshgrid(lx);
C0 = [xx(:) yy(:) zz(:)]; 
n3 = size(C0,1);
rd=meanRadius*(1+polydisperseRatio*rand(n3,1));
delta = (meanRadius/10)*(.5 - rand(n3,1));
C0 = C0 + delta;
%% Obtain the desired number (potential) of initial collisions
% this is done by performing an binary search over the multiplicative
% scaling on the initial configuration. Because the config is centered at 
% the origin, as long as the initial configuration is valid, the the
% desired scaling should be possible. 
lb = 0;
cur = 1;
ub = cur;
C = C0;
numCol = computePairwiseDistance(C,rd, colThresh);
while numCol > 0
    ub = ub*2;
    numCol = computePairwiseDistance(C,rd, colThresh);
    C = C0*ub;
end

if debug 
    figure;
    hold on 
    Gamma = lb:.01:ub;
    plot(Gamma, arrayfun(@(gamma) computePairwiseDistance(C0*gamma, rd, colThresh), Gamma), "Color",'blue','LineWidth',5)
    ylabel('Number of Collisions')
    xlabel('\gamma')
end

while true
    numCol = computePairwiseDistance(C,rd, colThresh);
    if desiredNumCol - desiredTol <= numCol && numCol <= desiredNumCol + desiredTol
        break
    elseif desiredNumCol - desiredTol < numCol
        lb = cur;
        cur = cur + (ub - cur)/2;
    else 
        ub = cur;
        cur = lb + (cur - lb)/2;
    end
    C = C0*cur;
end
fprintf('Initializing Configuration\n')
fprintf('- Desired number of (potential) collisions %d\n', desiredNumCol)
fprintf('- Tolerance on the desired number of collisions %d\n', desiredTol)
fprintf('- Actual numbrer of (potential) Collision %d\n', numCol)
%% Initial Particle Orientations
% set the initial direction to be towards the center 
% so that the particles want to pack together
init_dir = -C;
for i = 1:size(C,1)
    init_dir(i,:) = init_dir(i,:) / norm(init_dir(i,:));
end

if debug 
    figure
    ax = axes;
    ax.FontSize = 16;
    grid on
    view(0,90)
    R = max(arrayfun(@(i) norm(C(i,:)), 1:n));
    lims = [-2*R 2*R];
    ax = gca;
    xlim(ax,lims);
    ylim(ax, lims);
    zlim(ax, lims);
    grp = hgtransform(Parent=ax);
    hndls = cell(n3, 1);
    for k=1:n3
        [Sx, Sy, Sz]=sphere(40);
        Sx=rd(k)*Sx;
        Sy=rd(k)*Sy;
        Sz=rd(k)*Sz;
        Sc=C(k,:);
        hndls{k} = surf(Sx+Sc(1),Sy+Sc(2),Sz+Sc(3),Parent=grp,EdgeColor='k');
    end
    drawnow

    [~, pairwiseDistance] = computePairwiseDistance(C,rd, colThresh);
    figure;
    hold on
    imagesc(pairwiseDistance)
    colormap(hsv(512))
    colorbar
end

function [numCol, pairwiseDistance] = computePairwiseDistance(C,rd, colThresh)
n = size(C,1);
pairwiseDistance = Inf*ones(n,n);
for i = 2:n
    for j = i+1:n
        pairwiseDistance(i,j) = norm(C(i,:) - C(j,:)) - rd(i) - rd(j);
    end
end
% matches logic in LOCAL_check_collision_sph
numCol = sum(pairwiseDistance(:) < 1.1*colThresh); 
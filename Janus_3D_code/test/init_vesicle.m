function [C, init_dir] = init_vesicle(n, Cdst)
if ~exist('n','var') || isempty(n)
    n=100; 
end
if ~exist('Cdst','var') || isempty(Cdst)
    Cdst=2.5;
end

plotFlag = true;
ratioOut = 2/3;
nOut = round(n * ratioOut);
nIn = n - nOut;
Cin = generateEvenlySpacedPointsOnSpere(nIn, [], Cdst, plotFlag);
Cout = generateEvenlySpacedPointsOnSpere(nOut, norm(Cin(1,:)), Cdst, plotFlag);
C = [Cout ; Cin];
if plotFlag
    % Plot Spheres
    R = norm(C(1,:));
    plot_spheres({C}, [-2*R 2*R -2*R 2*R -2*R 2*R])
    minDist = computePairwiseDistance(C);
    figure
    histogram(minDist)
    title('MinDist of the final result')
end
% Point the amphiphilic tails of each particle on the outer shell towards the origin
init_dir_out= -Cout;
for i = 1:size(Cout,1)
    init_dir_out(i,:) = init_dir_out(i,:) / norm(init_dir_out(i,:));
end
% Point the amphiphilic tails of each particle on the inner shell away from the origin
init_dir_in= Cin;
for i = 1:size(Cin,1)
    init_dir_in(i,:) = init_dir_in(i,:) / norm(init_dir_in(i,:));
end
init_dir = [init_dir_out ; init_dir_in];
function C = generateEvenlySpacedPointsOnSpere(n, innerRadius, Cdst, plotFlag)
%% Generate points on the surface of a unit sphere
% using spherical fibonacci grid points
C = zeros(n,3);
offset = 0.5;
ix = 0:n-1 + offset;
phi = acos(1 - 2*ix/n);
theta = pi * (1 + sqrt(5)) * ix;
C(:,1) = cos(theta) .* sin(phi);
C(:,2) = sin(theta) .* sin(phi);
C(:,3) = cos(phi);
%% Compute pairwise distance
minDist = computePairwiseDistance(C);
% Plot minimum distance
if plotFlag
    figure()
    subplot(1,2,1)
    histogram(minDist)
    title('MinDist before removal outliers')
end
%% Remove trouble makers
xbar = mean(minDist);
sigma = std(minDist);
ToRemove = false(n,1);
for i = 1:n-1
    if abs(minDist(i) - xbar) > 2*sigma
        ToRemove(i) = true;
    end
end
C(ToRemove,:) = [];
minDist = computePairwiseDistance(C);
% Plot Final result minimum distance
if plotFlag
    subplot(1,2,2)
    histogram(minDist)
    title('MinDist after removal of outliers')
end
d = min(minDist);
R = Cdst / d;
if ~isempty(innerRadius)
    if R - innerRadius < Cdst
        R = R - (R - innerRadius) + Cdst;
    end
end
C = C*R;

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
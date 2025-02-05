function[x1, x2, distances, neighbors] = allWithPairs(par, tol, maxIter, eta, alg)

if strcmp(alg, 'mb')
    pairDistanceFunction = str2func('movingBallsPair');
elseif strcmp(alg,'GJK')
    pairDistanceFunction = str2func('GJKPair');
else
    disp('No algorithm specified, terminating.');
    return;
end

%number of spheroids
n = length(par); 

% Initialize arrays 
centers = zeros(n,3); 
longestSemis = zeros(n,1);
 
%Find the centers and largest semi axis of each ellipsoid
for i=1:n
    centers(i,:) = par(i).C; 
    longestSemis(i) = max([par(i).a par(i).b par(i).c]); 
end

%Treat as spheres and compute spherical distances, if distances are close enough use the moving balls algorithm and find a list of neighbors and corresponding contact points. For efficiency, only the contact points of particles with higher indices are computed at each time.
distances = zeros(n,n);
x1 = cell(n,1);
x2 = cell(n,2);
neighbors = cell(n,1);
for i = 1:n 
    x1Local = [];
    x2Local = [];
    neighborsLocal = [];
    for j = i + 1:n
        centerDistance = norm(centers(i) - centers(j));
        semiDistance = longestSemis(i) + longestSemis(j);
        if centerDistance >= longestSemis(i) + longestSemis(j)
            distances(i,j) = centerDistance - semiDistance;
        end
        longestSemi = max([longestSemis(i) longestSemis(j)]);
        if (distances(i,j) <= 2*eta*longestSemi) 
            [x1New, x2New, distance] = pairDistanceFunction(par(i), par(j), tol, maxIter);
            neighborsLocal =[neighborsLocal; j];
            x1Local = [x1Local; x1New];
            x2Local = [x2Local; x2New];
            distances(i,j) = distance;
        end
    end
    x1{i} = x1Local;
    x2{i} = x2Local;
    neighbors{i} = neighborsLocal;
end

end

    
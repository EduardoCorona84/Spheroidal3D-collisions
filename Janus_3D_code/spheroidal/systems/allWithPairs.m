function[x1, x2, distances, neighbors] = allWithPairs(par, tol, maxIter, eta, alg)

if strcmp(alg, 'mb')
    pairDistanceFunction = str2func('movingBallsPair');
elseif strcmp(alg,'GJKJohn')
    pairDistanceFunction = str2func('GJKJohnsonPair');
elseif strcmp(alg, 'GJKSigned')
    pairDistanceFunction = str2func('GJKSignedVolumesPair');
elseif strcmp(alg,'GJKJohnN')
    pairDistanceFunction = str2func('GJKJohnsonNestPair');
elseif strcmp(alg, 'GJKSignedN')
    pairDistanceFunction = str2func('GJKSignedVolumesNestPair');
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

%Treat as spheres and compute spherical distances, if distances are close enough use a distance algorithm and find a list of neighbors and corresponding contact points. 
distances = zeros(n,n);
x1 = cell(n,n);
x2 = cell(n,n);
neighbors = logical(false(n,n));
for i = 1:n 
    for j = i + 1:n
        sphericalDiameter = longestSemis(i) + longestSemis(j);
        sphericalDistance = norm(centers(i) - centers(j)) - sphericalDiameter;
        if sphericalDistance <= eta*sphericalDiameter
            [x1New, x2New, distance] = pairDistanceFunction(par(i), par(j), tol, maxIter, false);
            distances(i,j) = distance;
            x1{i,j} = x1New;
            x2{j,i} = x1New;
            x1{j,i} = x2New;
            x2{i,j} = x2New;
            neighbors(i,j) = 1;
            neighbors(j,i) = 1;
        else
            distances(i,j) = sphericalDistance;
        end
        distances(j,i) = distances(i,j);
    end

end

end

    
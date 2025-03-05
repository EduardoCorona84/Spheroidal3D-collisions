function [x1, x2, distances, neighbors] = allWithWall(par, eta)
    n = length(par); 

    % Initialize arrays 
    centers = zeros(n,3); 
    longestSemis = zeros(n,1);
    
    %Find the centers and largest semi axis of each ellipsoid
    for i=1:n
        centers(i,:) = par(i).C; 
        longestSemis(i) = max([par(i).a par(i).b par(i).c]); 
    end

    distances = zeros(n,1);
    x1 = cell(n,1);
    x2 = cell(n,1);
    neighbors = logical(false(n,1));

    for i=1:n
        [x1New, x2New, d] = planeEllipsoidIntersection(par(i), [0 0 1].', [0 0 0].');
        distances(i) = d;
        if d<= eta*longestSemis(i)
            x1{i} = x1New;
            x2{i} = x2New;
            neighbors(i) = true;
        end
    end
end
function [weights, newBitString, v, vNorm] = signedVolumes(Y, card)

    %This implements the signed volumes distance subroutine for the GJK algorithm. This is from montanari et al, Improving the GJK Algorithm for Faster and More Reliable Distance Queries Between Convex Objects.

    %Checks the cardinality of the active subset
    r = card - 1;
    if r == 3
        [weights, newBitString, v, vNorm] = S3D(Y);
    elseif r == 2
        currBitString = logical([1 1 1 0]);
        [weights, newBitString, v, vNorm] = S2D(Y, currBitString);
    elseif r == 1
        currBitString = logical([1 1 0 0]);
        [weights, newBitString, v, vNorm] = S1D(Y, currBitString);
    else 
        weights = 1;
        newBitString = logical([1 0 0 0]);
        v = Y(:,1);
        vNorm = norm(v);
    end

end


function [weights, newBitString, v, vNorm] = S3D(Y)

    M = [Y ; 1 1 1 1];
    C = zeros(4,1);
    detM = 0;
    for j = 1:4
        C(j) = (-1)^(j + 4) * det3by3(M(1:3, 1:4 ~= j));
        detM = detM + C(j);
    end
    if all(detM.*C > 0)
        weights = C./detM;
        newBitString = logical([1 1 1 1]);
        v = Y*weights;
        vNorm = norm(v);
        return;
    else
        vNorm = Inf;
        for j = 1:3
            if -detM*C(j) > 0 
                tempBitString = logical([1 1 1 1]);
                tempBitString(j) = 0;
                [weightsCandidate, newBitStringCandidate, vCandidate, vNormCandidate] = S2D(Y, tempBitString);
                if vNormCandidate < vNorm 
                    weights = weightsCandidate;
                    newBitString = newBitStringCandidate;
                    v = vCandidate;
                    vNorm = vNormCandidate;
                end
            end
        end
    end
end

function [weights, newBitString, v, vNorm] = S2D(Y, currBitString)

    activeIndices = find(currBitString);
    n = cross(Y(:,activeIndices(2)) - Y(:,activeIndices(3)), Y(:, activeIndices(1)) - Y(:, activeIndices(3)));
    p0 = dot(Y(:,activeIndices(3)), n).*n./dot(n,n);
    muMax = 0;
    for i = 1:3
        mu = (-1)^i*n(i);
        if abs(mu) > abs(muMax)
            muMax = mu;
            J = i;
        end
        
    end
    keptCoordinates = find(1:3 ~= J);
    k = 2;
    l = 3;
    C = zeros(3,1);
    %Calculate the determinants of the matrices with p0 substituted in required for cramers rule.
    for i = 1:3
        C(i) = p0(keptCoordinates(1))*Y(keptCoordinates(2), activeIndices(k)) + p0(keptCoordinates(2))*Y(keptCoordinates(1),activeIndices(l)) + Y(keptCoordinates(1),activeIndices(k))*Y(keptCoordinates(2), activeIndices(l)) - p0(keptCoordinates(1))*Y(keptCoordinates(2), activeIndices(l)) - p0(keptCoordinates(2))*Y(keptCoordinates(1), activeIndices(k)) - Y(keptCoordinates(2), activeIndices(k))*Y(keptCoordinates(1), activeIndices(l));
        k = l;
        l = i;  
    end
    if all(muMax.*C > 0) 
        weights = C./muMax;
        newBitString = currBitString;
        v = Y(:,newBitString)*weights;
        vNorm = norm(v);
    else
        vNorm = Inf;
        for j = 1:2
            if -muMax*C(j) > 0
                tempBitString = currBitString;
                tempBitString(activeIndices(j)) = 0;
                [weightsCandidate, newBitStringCandidate, vCandidate, vNormCandidate] = S1D(Y, tempBitString);
                if vNormCandidate < vNorm
                    weights = weightsCandidate;
                    newBitString = newBitStringCandidate;
                    v = vCandidate;
                    vNorm = vNormCandidate;
                end
            end
        end
    end
end

function [weights, newBitString, v, vNorm] = S1D(Y, currBitString)
    activeIndices = find(currBitString);
    t = Y(:,activeIndices(1)) - Y(:,activeIndices(2));
    p0 = -(dot(Y(:,activeIndices(2)), t)/dot(t,t)).*t + Y(:, activeIndices(2));
    muMax = 0;
    for i = 1:3
        mu = Y(i, activeIndices(1)) - Y(i, activeIndices(2));
        if abs(mu) > abs(muMax)
            muMax = mu;
            I = i;
        end
    end
    k = 2;
    C = zeros(2,1);
    for j=1:2
        C(j) = (-1)^j*(Y(I,activeIndices(k)) - p0(I));
        k = j;
    end
    if all(muMax.*C > 0)
        weights = C./muMax;
        newBitString = currBitString;
        v = Y(:, newBitString)*weights;
        vNorm = norm(v);
    else
        weights = 1;
        newBitString = currBitString;
        newBitString(activeIndices(1)) = 0;
        v = Y(:,newBitString);
        vNorm = norm(v);
    end
end


function det = det3by3(M)
    det = M(1,1)*(M(2,2)*M(3,3) - M(2,3)*M(3,2)) - M(1,2)*(M(2,1)*M(3,3) - M(2,3)*M(3,1)) + M(1,3)*(M(2,1)*M(3,2) - M(2,2)*M(3,1));
    
end

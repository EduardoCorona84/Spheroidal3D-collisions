function [x1, x2, distance, distanceIters] = GJKSignedVolumesPair(par1, par2, tol, maxIter, val)

    %Get parameters of two ellipsoids.
    C1 = par1.C; R1 = par1.R; a1 = par1.a; b1=par1.b; c1=par1.c; 
    C2 = par2.C; R2 = par2.R; a2 = par2.a; b2=par2.b; c2=par2.c;

    %Check if row vectors, if so make column vectors.
    if isequal(size(C1),[1 3])
    C1 = C1.';
    C2 = C2.';
    end

    %Allocate Space for needed arrays
    distanceIters = zeros(maxIter);
    P = zeros(3,4);
    Q = zeros(3,4);
    Y = zeros(3,4);

    if val 
        v = [2*C1(1) 2*b1 2*b1 ].';
    else
        v = C1 - C2;
    end

    %Find norm of the initial point and initialize needed values
    vNorm = norm(v);
    mu0 = 0;
    card = 0;
    weights = 1;


    for iter = 0:maxIter 
        
        %record distances for order calculation
        distanceIters(iter + 1) = vNorm;

        %Find open index and populate needed arrays
        newIndex = card + 1;
        P(:, newIndex) = ellipsoidSupportMapping(-v, a1, b1, c1, C1, R1);
        Q(:, newIndex) = ellipsoidSupportMapping(v, a2, b2, c2, C2, R2);
        Y(:, newIndex) = P(:, newIndex) - Q(:, newIndex);

        %Calculate error
        delta = dot(v, Y(:, newIndex))/vNorm;
        mu1 = max([mu0, delta]);
        mu0 = mu1;

        %If within tolerance, find values and return
        if vNorm - mu0 < tol
            if iter == 0
                x1 = P(:,1);
                x2 = Q(:,1);
            else
                x1 = P(:, 1:card)*weights;
                x2 = Q(:, 1:card)*weights;
            end
            distance = vNorm;
            return;
        end

        %Compute weights, active subset of vertices, and the new distance using the signed volume algorithm.
        [weights, activeBitString, v, vNorm] = signedVolumes(Y, newIndex);

        card = sum(activeBitString);
        P(:, 1:card) = P(:, activeBitString);
        Q(:, 1:card) = Q(:, activeBitString);
        Y(:, 1:card) = Y(:, activeBitString);
        

            %Check if overlapping
        if vNorm < tol
            distanceIters(iter + 1) = vNorm; 
            x1 = P(:,1:card)*weights;
            x2 = Q(:,1:card)*weights;
            distance = 0;
            return
        end

    end 
    
end
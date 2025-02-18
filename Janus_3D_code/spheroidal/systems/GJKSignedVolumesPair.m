function [x1, x2, d] = GJKSignedVolumesPair(par1, par2, tol, maxIter, val)

    %Get parameters of two ellipsoids.
    C1 = par1.C; R1 = par1.R; a1 = par1.a; b1=par1.b; c1=par1.c; 
    C2 = par2.C; R2 = par2.R; a2 = par2.a; b2=par2.b; c2=par2.c;

    %Check if row vectors, if so make column vectors.
    if isequal(size(C1),[1 3])
    C1 = C1.';
    C2 = C2.';
    end

    %Allocate Space for needed arrays
    P = zeros(3,4);
    Q = zeros(3,4);
    Y = zeros(3,4);

    %{
    Perform part of the 0th iteration
    --------------------------------------------
    %}

    if val 
        v = [2*C1(1) 2*b1 2*b1 ].';
    else
        v = C1 - C2;
    end

    vNorm = norm(v);
    p = ellipsoidSupportMapping(-v, a1, b1, c1, C1, R1);
    q = ellipsoidSupportMapping(v, a2, b2, c2, C2, R2);
    activeBitString = logical([0 0 0 0]);
    mu0 = 0;

    %Check if centroids are overlapping, if so, return.
    if vNorm < tol
        x1 = p;
        x2 = q;
        d = 0;
        return
    end

    %{
    --------------------------------------------
    End 0th Iteration
    %}

    %Begin while loop iterations
    iter = 0;
    while vNorm - mu0 > tol && iter < maxIter
        if iter == 0
            P(:,1) = p;
            Q(:,1) = q;
            v = p - q;
            Y(:,1) = v;
            vNorm = norm(v);
            weights = 1;
            card = 1;
            activeBitString(1) = 1;
        else
            P(:, card + 1) = p;
            Q(:, card + 1) = q;
            Y(:, card + 1) = p - q;
            card = card + 1;
            %Compute weights, active subset of vertices, and the new distance using the signed volume algorithm.
            [weights, activeBitString, v, vNorm] = signedVolumes(Y, card);
            card = sum(activeBitString);
            P(:, 1:card) = P(:, activeBitString);
            Q(:, 1:card) = Q(:, activeBitString);
            Y(:, 1:card) = Y(:, activeBitString);
        
        end

            %Check if overlapping
        if vNorm < tol
            x1 = P(:,1:card)*weights;
            x2 = Q(:,1:card)*weights;
            d = 0;
            return
        end

        
        %Find next support point and compute the error.
        p = ellipsoidSupportMapping(-v, a1, b1, c1, C1, R1);
        q = ellipsoidSupportMapping(v, a2, b2, c2, C2, R2);
        w = p - q;
        delta = dot(v, w)/vNorm;
        mu1 = max([mu0, delta]);
        iter = iter + 1;
        mu0 = mu1;
    end 
    %Compute final quantities.
    x1 = P(:,1:card)*weights;
    x2 = Q(:,1:card)*weights;
    d = vNorm;
    
end
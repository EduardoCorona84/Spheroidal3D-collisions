function [x1, x2, d] = GJKSignedVolumesNestPair(par1, par2, tol, maxIter, val)

    %Get parameters
    C1 = par1.C; R1 = par1.R; a1 = par1.a; b1=par1.b; c1=par1.c; 
    C2 = par2.C; R2 = par2.R; a2 = par2.a; b2=par2.b; c2=par2.c;

    %check if row vectors
    if isequal(size(C1),[1 3])
        C1 = C1.';
        C2 = C2.';
    end

    %Allocate Space for needed arrays
    P = zeros(3,4);
    Q = zeros(3,4);
    Y = zeros(3,4);


    if val
        xCurr = [2*C1(1) 2*b1 2*b1 ].';
    else
        xCurr = C1 - C2;
    end

    sOld = xCurr;
    dOld = xCurr;
    card = 0;
    criterion = false;

    for iter = 0:maxIter
        alpha = (iter + 1)/(iter + 3);
        y = alpha.*xCurr + (1 - alpha).*sOld;

        if criterion
            dCurr = xCurr;
        else
            dCurr = alpha.*dOld + 2.*(1 - alpha).*y;
        end

        P(:, card + 1) = ellipsoidSupportMapping(-dCurr, a1, b1, c1, C1, R1);
        Q(:, card + 1) = ellipsoidSupportMapping(dCurr, a2, b2, c2, C2, R2);

        sCurr = P(:, card + 1) - Q(:, card + 1);
        if 2*dot(xCurr, xCurr - sCurr) <= tol
            if norm(dCurr - xCurr) <= tol
                x1 = P(:, 1:card)*weights;
                x2 = Q(:, 1:card)*weights;
                d = xNewNorm;
                return
            else
                P(:, card + 1) = ellipsoidSupportMapping(-2.*xCurr, a1, b1, c1, C1, R1);
                Q(:, card + 1) = ellipsoidSupportMapping(2.*xCurr, a2, b2, c2, C2, R2);
                sCurr = P(:, card + 1) - Q(:, card + 1);
                criterion = true;
            end
        end 

        Y(:, card + 1) = sCurr;
        card = card + 1;
        if iter == 0
            xNew = sCurr;
            xNewNorm = norm(xNew);
        else
            %Compute weights, active subset of vertices, and the new distance using the signed volume algorithm.
            [weights, activeBitString, xNew, xNewNorm] = signedVolumes(Y, card);
            card = sum(activeBitString);
            P(:, 1:card) = P(:, activeBitString);
            Q(:, 1:card) = Q(:, activeBitString);
            Y(:, 1:card) = Y(:, activeBitString);

        end
            %Check if overlapping
        if xNewNorm < tol
            x1 = P(:,1:card)*weights;
            x2 = Q(:,1:card)*weights;
            d = 0;
            return
        end

        sOld = sCurr;
        dOld = dCurr;
        xCurr = xNew;
    end 
end
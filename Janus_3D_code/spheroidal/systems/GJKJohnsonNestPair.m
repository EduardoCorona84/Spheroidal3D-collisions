function [x1, x2, distance, distanceIters] = GJKJohnsonNestPair(par1, par2, tol, maxIter, val)

    %Get parameters
    C1 = par1.C; R1 = par1.R; a1 = par1.a; b1=par1.b; c1=par1.c; 
    C2 = par2.C; R2 = par2.R; a2 = par2.a; b2=par2.b; c2=par2.c;

    %check if row vectors
    if isequal(size(C1),[1 3])
        C1 = C1.';
        C2 = C2.';
    end

    %Allocate Space for needed arrays
    distanceIters = zeros(maxIter);
    P = zeros(3,4);
    Q = zeros(3,4);
    Y = zeros(3,4);
    dotProducts = zeros(4,4);

    DeltaX = zeros(16,4);
    DeltaX(2,1) = 1;
    DeltaX(3,2) = 1;
    DeltaX(5,3) = 1;
    DeltaX(9,4) = 1;

    %Arrays used for indexing
    bitStrings = dec2bin(0:2^4 - 1) - '0';
    bitStrings = logical(fliplr(bitStrings));
    cardIndices = [2 3 5 9 4 6 7 10 11 13 8 12 14 15 16];
    binaryWeights = [1 2 4 8];
    activeBitString = logical([0 0 0 0]);

    %Check if validation case
    if val
        v = [2*C1(1) 2*b1 2*b1 ].';
    else
        v = C1 - C2;
    end

    %Find norm of the initial point and initialize needed values
    vNorm = norm(v);
    s = v;
    d = v;
    criterion = true;

    for iter = 0:maxIter

        distanceIters(iter + 1) = vNorm;

        %Check if using nesterov accelerated or reverting back to standard GJK.
        if criterion
            alpha = (iter + 1)/(iter + 3);
            y = alpha.*v + (1 - alpha).*s;
            d = alpha.*d + 2.*(1 - alpha).*y;
        else
            d = v;
        end

        %Find open index and populate needed arrays
        newIndex = find(~activeBitString, 1, 'first');
        P(:, newIndex) = ellipsoidSupportMapping(-d, a1, b1, c1, C1, R1);
        Q(:, newIndex) = ellipsoidSupportMapping(d, a2, b2, c2, C2, R2);
        activeBitString(newIndex) = 1;

        s = P(:, newIndex) - Q(:, newIndex);
        %Check convergence criterion (duality gap).
        if 2*dot(v, v - s) <= tol
            %Check fixed point condition
            if norm(d - v) <= tol
                decimal = dot(activeBitString, binaryWeights) + 1;
                x1 = P(:,activeBitString)*DeltaX(decimal, activeBitString).'./sum(DeltaX(decimal, activeBitString));
                x2 = Q(:,activeBitString)*DeltaX(decimal, activeBitString).'./sum(DeltaX(decimal, activeBitString));
                distance = vNorm;
                return
            else
                P(:, newIndex) = ellipsoidSupportMapping(-2.*v, a1, b1, c1, C1, R1);
                Q(:, newIndex) = ellipsoidSupportMapping(2.*v, a2, b2, c2, C2, R2);
                s = P(:, newIndex) - Q(:, newIndex);
                criterion = false;
            end
        end 

        Y(:, newIndex) = s;

        %Compute dot products
        dotProducts(activeBitString, newIndex) = (Y(:,activeBitString).'*Y(:, newIndex)).';
        dotProducts(newIndex, activeBitString) = dotProducts(activeBitString, newIndex).';
        
        %Iterate up cardinality
        for j = cardIndices
            %Now check its a subset of the active set.
            ifSubset = activeBitString | bitStrings(j,:);
            %Check that the subset is a subset of the active set.
            if all(ifSubset == activeBitString)
                minIndex = find(bitStrings(j,:), 1, 'first');
                inactiveSubsetIndices = find(bitStrings(j,:) ~= activeBitString);
                %If the subset contains the new Index, calculate the Delta X values
                if bitStrings(j, newIndex) == 1
                    nonPositivityCondition = true;
                    for k = 1:length(inactiveSubsetIndices)
                        newBitString = bitStrings(j,:);
                        newBitString(inactiveSubsetIndices(k)) = 1;
                        newDecimal = dot(newBitString, binaryWeights) + 1;
                        DeltaX(newDecimal, inactiveSubsetIndices(k)) = dot(DeltaX(j, bitStrings(j,:)), dotProducts(minIndex, bitStrings(j,:)) - dotProducts(inactiveSubsetIndices(k), bitStrings(j,:)));

                        if DeltaX(newDecimal, inactiveSubsetIndices(k)) <= 0
                            newNonPositivityCondition = true;
                        else
                            newNonPositivityCondition = false;
                        end
                        nonPositivityCondition = nonPositivityCondition && newNonPositivityCondition;
                    end

                    positivityCondition =  all(DeltaX(j, bitStrings(j,:)) > 0);
                    
                    %Check if the solution conditions are satisfied
                    if positivityCondition && nonPositivityCondition
                        activeBitString = bitStrings(j,:);
                        v = Y(:,activeBitString)*DeltaX(j, activeBitString).'./sum(DeltaX(j,activeBitString));
                        vNorm = norm(v);
                        break;
                    end
                %If newIndex not a part of the subset, compute DeltaX for subsets containing the current subset (bitString) and newIndex.    
                else
                    newBitString = bitStrings(j,:);
                    newBitString(newIndex) = 1;
                    newDecimal = dot(newBitString, binaryWeights) + 1;
                    DeltaX(newDecimal, newIndex) = dot(DeltaX(j, bitStrings(j,:)), dotProducts(minIndex, bitStrings(j,:)) - dotProducts(newIndex, bitStrings(j,:)));
                end
            end
        end
    
        %Check if overlapping
        if vNorm < tol
            distanceIters(iter + 1) = vNorm;
            decimal = dot(activeBitString, binaryWeights) + 1;
            x1 = P(:,activeBitString)*DeltaX(decimal, activeBitString).'./sum(DeltaX(decimal, activeBitString));
            x2 = Q(:,activeBitString)*DeltaX(decimal, activeBitString).'./sum(DeltaX(decimal, activeBitString));
            distance = 0;
            return
        end
    end 
end

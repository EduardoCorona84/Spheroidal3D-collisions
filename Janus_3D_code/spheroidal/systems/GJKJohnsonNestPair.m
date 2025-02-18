function [x1, x2, d] = GJKJohnsonNestPair(par1, par2, tol, maxIter, val)

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

    bitStrings = dec2bin(0:2^4 - 1) - '0';
    bitStrings = logical(fliplr(bitStrings));
    cardIndices = [2 3 5 9 4 6 7 10 11 13 8 12 14 15 16];
    binaryWeights = [1 2 4 8];
    activeBitString = logical([0 0 0 0]);


    DeltaX = zeros(16,4);
    DeltaX(2,1) = 1;
    DeltaX(3,2) = 1;
    DeltaX(5,3) = 1;
    DeltaX(9,4) = 1;

    dotProducts = zeros(4,4);

    if val
        xCurr = [2*C1(1) 2*b1 2*b1 ].';
    else
        xCurr = C1 - C2;
    end

    sOld = xCurr;
    dOld = xCurr;
    criterion = false;

    for iter = 0:maxIter
        alpha = (iter + 1)/(iter + 3);
        y = alpha.*xCurr + (1 - alpha).*sOld;

        if criterion
            dCurr = xCurr;
        else
            dCurr = alpha.*dOld + 2.*(1 - alpha).*y;
        end

        newIndex = find(~activeBitString, 1, 'first');
        P(:, newIndex) = ellipsoidSupportMapping(-dCurr, a1, b1, c1, C1, R1);
        Q(:, newIndex) = ellipsoidSupportMapping(dCurr, a2, b2, c2, C2, R2);
        activeBitString(newIndex) = 1;

        sCurr = P(:, newIndex) - Q(:, newIndex);
        if 2*dot(xCurr, xCurr - sCurr) <= tol
            if norm(dCurr - xCurr) <= tol
                decimal = dot(activeBitString, binaryWeights) + 1;
                x1 = P(:,activeBitString)*DeltaX(decimal, activeBitString).'./sum(DeltaX(decimal, activeBitString));
                x2 = Q(:,activeBitString)*DeltaX(decimal, activeBitString).'./sum(DeltaX(decimal, activeBitString));
                d = xNewNorm;
                return
            else
                P(:, newIndex) = ellipsoidSupportMapping(-2.*xCurr, a1, b1, c1, C1, R1);
                Q(:, newIndex) = ellipsoidSupportMapping(2.*xCurr, a2, b2, c2, C2, R2);
                sCurr = P(:, newIndex) - Q(:, newIndex);
                criterion = true;
            end
        end 

        Y(:, newIndex) = sCurr;

        if iter == 0
            dotProducts(1,1) = dot(sCurr,sCurr);
            xNew = sCurr;
            xNewNorm = norm(xNew);
        else

            %Now begin the Johnson distance algorithm.
            %update bit string, compute dot products, and find cardinality
            dotProducts(activeBitString, newIndex) = (Y(:,activeBitString).'*Y(:, newIndex)).';
            dotProducts(newIndex, activeBitString) = dotProducts(activeBitString, newIndex).';
            
            %Iterate up cardinality
            for j = cardIndices
                %Now check each subset
                ifSubset = activeBitString | bitStrings(j,:);
                %Check the size of the subset and that the subset is a subset of the active set. 
                if all(ifSubset == activeBitString)
                    minIndex = find(bitStrings(j,:), 1, 'first');
                    inactiveSubsetBitString = bitStrings(j,:) ~= activeBitString;
                    inactiveSubsetIndices = find(inactiveSubsetBitString);
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
                        
                        if positivityCondition && nonPositivityCondition
                            activeBitString = bitStrings(j,:);
                            xNew = Y(:,activeBitString)*DeltaX(j, activeBitString).'./sum(DeltaX(j,activeBitString));
                            xNewNorm = norm(xNew);
                            break;
                        end
                    else
                        newBitString = bitStrings(j,:);
                        newBitString(newIndex) = 1;
                        newDecimal = dot(newBitString, binaryWeights) + 1;
                        DeltaX(newDecimal, newIndex) = dot(DeltaX(j, bitStrings(j,:)), dotProducts(minIndex, bitStrings(j,:)) - dotProducts(newIndex, bitStrings(j,:)));
                    end
                end
            end
        
        end
            %Check if overlapping
        if xNewNorm < tol
            decimal = dot(activeBitString, binaryWeights) + 1;
            x1 = P(:,activeBitString)*DeltaX(decimal, activeBitString).'./sum(DeltaX(decimal, activeBitString));
            x2 = Q(:,activeBitString)*DeltaX(decimal, activeBitString).'./sum(DeltaX(decimal, activeBitString));
            d = 0;
            return
        end
        sOld = sCurr;
        dOld = dCurr;
        xCurr = xNew;
    end 
end

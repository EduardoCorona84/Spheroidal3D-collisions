function [x1, x2, d] = GJKJohnsonPair(par1, par2, tol, maxIter, val)

    %Get parameters
    C1 = par1.C; R1 = par1.R; a1 = par1.a; b1=par1.b; c1=par1.c; 
    C2 = par2.C; R2 = par2.R; a2 = par2.a; b2=par2.b; c2=par2.c;

    %Check if row vectors
    if isequal(size(C1),[1 3])
        C1 = C1.';
        C2 = C2.';
    end

    %Allocate space for needed arrays
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

    %{
    Perform the first iteration
    ---------------------------------
    %}

    %Check if validation case
    if val 
        v = [2*C1(1) 2*b1 2*b1 ].';
    else
        v = C1 - C2;
    end

    %Find initial point and its norm
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
    End some parts of 0th Iteration
    %}

    iter = 0;

    while vNorm - mu0 > tol && iter < maxIter

        if iter == 0
            P(:,1) = p;
            Q(:,1) = q;
            v = p - q;
            Y(:,1) = v;
            dotProducts(1,1) = dot(v,v);
            vNorm = norm(v);
            activeBitString(1) = 1;
            
        else
            newIndex = find(~activeBitString, 1, 'first');
            P(:, newIndex) = p;
            Q(:, newIndex) = q;
            Y(:, newIndex) = P(:, newIndex) - Q(:, newIndex);

            %Now begin the Johnson distance algorithm.
            %update bit string, compute dot products, and find cardinality
            activeBitString(newIndex) = 1;
            dotProducts(activeBitString, newIndex) = (Y(:,activeBitString).'*Y(:, newIndex)).';
            dotProducts(newIndex, activeBitString) = dotProducts(activeBitString, newIndex).';
            
            %Iterate up cardinality
            for j = cardIndices
                %Now check each subset
                ifSubset = activeBitString | bitStrings(j,:);
                %Check the size of the subset and that the subset is a subset of the active set. 
                if all(ifSubset == activeBitString)
                    minIndex = min(activeSubsetIndices);
                    inactiveSubsetIndices = find(bitStrings(j,:) ~= activeBitString);
                    if bitStrings(j, newIndex) == 1
                        nonPositivityCondition = true;
                        for k = 1:length(inactiveSubsetIndices)
                            newBitString = bitStrings(j,:);
                            newBitString(inactiveSubsetIndices(k)) = 1;
                            newDecimal = dot(newBitString, binaryWeights) + 1;
                            DeltaX(newDecimal,inactiveSubsetIndices(k)) = dot(DeltaX(j, bitStrings(j,:)), dotProducts(minIndex, bitStrings(j,:)) - dotProducts(inactiveSubsetIndices(k), bitStrings(j,:)));

                            if DeltaX(newDecimal,inactiveSubsetIndices(k)) <= 0
                                newNonPositivityCondition = true;
                            else
                                newNonPositivityCondition = false;
                            end
                            nonPositivityCondition = nonPositivityCondition && newNonPositivityCondition;
                        end

                        positivityCondition =  all(DeltaX(j, bitStrings(j,:)) > 0);
                        
                        if positivityCondition & nonPositivityCondition
                            activeBitString = bitStrings(j,:);
                            v = Y(:,activeBitString)*DeltaX(j, activeBitString).'./sum(DeltaX(j,activeBitString));
                            vNorm = norm(v);
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
        if vNorm < tol
            decimal = dot(activeBitString, binaryWeights) + 1;
            weights = DeltaX(decimal, activeBitString).'./sum(DeltaX(decimal, activeBitString));
            x1 = P(:,activeBitString)*weights;
            x2 = Q(:,activeBitString)*weights;
            d = 0;
            return
        end
        
        %Find next support point and compute the error.
        p = ellipsoidSupportMapping(-v, a1, b1, c1, C1, R1);
        q = ellipsoidSupportMapping(v, a1, b1, c1, C2, R2);
        w = p - q;
        delta = dot(v, w)/vNorm;
        mu1 = max([mu0, delta]);
        iter = iter + 1;
        mu0 = mu1;
    end 
    %Compute final quantities.
    decimal = dot(activeBitString, binaryWeights) + 1;
    weights = DeltaX(decimal, activeBitString).'./sum(DeltaX(decimal, activeBitString));
    x1 = P(:,activeBitString)*weights;
    x2 = Q(:,activeBitString)*weights;
    d = vNorm;
end
 


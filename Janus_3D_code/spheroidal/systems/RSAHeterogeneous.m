function spheroids = RSAHeterogeneous(cubeLength, Ar, Deq, parameter)

    %Randomized Sequential Adsorption: This algorithm generates spheroids within some sort of container. This is done by randomly generating one, then checking if it intersects another spheroid or the container.  
    
    volumeFraction = parameter*(4/3)*pi*(Deq/2)^3/(cubeLength^3);
    fprintf('Desired Volume Fraction:%f\n', volumeFraction);
    
    a = cell(parameter, 1);
    b = cell(parameter, 1);
    c = cell(parameter, 1);
    C = cell(parameter, 1);
    R = cell(parameter, 1);
    
    total = 0;
    ratios = size(Ar);
    
    while total < parameter
        ArCurr = Ar(randi(ratios));
        a1 = ArCurr^(2/3)*(Deq/2);
        b1 = ArCurr^(-1/3)*(Deq/2);
        c1 = b1; 
        C1 = rand(3,1).*cubeLength - (cubeLength/2);
        [R1,~] = qr(randn(3,3));
    
        overlap = boxOverlap(cubeLength, a1, b1, c1, C1, R1);
        if overlap == true
            continue;
        end
    
        if total == 0
            total = total + 1;
            a{total} = a1;
            b{total} = b1;
            c{total} = c1;
            C{total} = C1;
            R{total} = R1;
        else
            d = 1;
            for i = 1:total
                par1 = struct('a', a1, 'b', b1, 'c', c1, 'C', C1, 'R', R1);
                par2 = struct('a', a{i}, 'b', b{i}, 'c', c{i}, 'C', C{i}, 'R', R{i});
                [~, ~, dNew] = movingBallsPair(par1, par2, 1e-5, 500);
                d = d*dNew;
                if d == 0
                    break;
                end
            end
    
            if d == 0
                continue;
            else
                total = total + 1;
                a{total} = a1;
                b{total} = b1;
                c{total} = c1;
                C{total} = C1;
                R{total} = R1;
            end
        end
    end
    
    spheroids = struct('a', a, 'b', b, 'c', c, 'C', C, 'R', R);
    
    end
    
    function overlap = boxOverlap(L, a, b, c, C, R)
        %L is the cube side length.
        %Determines if an ellipsoid overlaps with a box. 
        %a, b, c, are the semi axes of the ellipsoid
        %C and R are the center and rotation matrix 
        overlap = false;
    
        nx = [L/2 0 0].';
        suppMap = ellipsoidSupportMapping(nx, a, b, c , C, R);
        if suppMap(1) - L/2 > 0
            overlap = true;
            return;
        end
        suppMap = ellipsoidSupportMapping(-nx, a, b, c , C, R);
        if suppMap(1) + L/2 < 0
            overlap = true;
            return;
        end
    
        ny = [0 L/2 0].';
        suppMap = ellipsoidSupportMapping(ny, a, b, c , C, R);
        if suppMap(2) - L/2 > 0
            overlap = true;
            return;
        end
    
        suppMap = ellipsoidSupportMapping(-ny, a, b, c , C, R);
        if suppMap(2) + L/2 < 0
            overlap = true;
            return;
        end
    
        nz = [0 0 L/2].';
        suppMap = ellipsoidSupportMapping(nz, a, b, c , C, R);
        if suppMap(3) - L/2 > 0
            overlap = true;
            return;
        end
    
        suppMap = ellipsoidSupportMapping(-nz, a, b, c , C, R);
        if suppMap(3) + L/2 < 0
            overlap = true;
            return;
        end
    
    end
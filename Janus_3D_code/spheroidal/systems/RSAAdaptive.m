function spheroids = RSAHeterogeneous(cubeLength, Ar, Deq, parameter)

    %Randomized Sequential Adsorption: This algorithm generates spheroids within some sort of container. This is done by randomly generating one, then checking if it intersects another spheroid or the container. An array of aspects ratio can be passed to the function, and the function will randomly select one of these aspect ratios for each spheroid or the desired volume fraction is reached.
    
    if parameter < 1
        total = ceil(parameter*(4/3*pi*(Deq/2)^3/(cubeLength^3))^(-1));
        volumeFraction = parameter;
    else
        volumeFraction = parameter*(4/3)*pi*(Deq/2)^3/(cubeLength^3);
        total = parameter;
    end
    fprintf('Desired Volume Fraction:%f\n', volumeFraction);
    fprintf('Desired Number of Spheroids:%f\n', total);
    
    a = cell(total, 1);
    b = cell(total, 1);
    c = cell(total, 1);
    C = cell(total, 1);
    R = cell(total, 1);
    
    currTotal = 0;
    ratios = size(Ar);
    
    while currTotal < total
        DeqCurr = (10/9)*Deq;
        C1 = rand(3,1).*cubeLength - (cubeLength/2);
        [R1,~] = qr(randn(3,3));


        %keep trying for this specific aspect ratio until the spheroid is within the box and not overlapping with any other spheroid
        accepted = false;
        while accepted == false
            DeqCurr = DeqCurr*0.9;
            if DeqCurr < (1/2)*Deq
                break;
            end
            a1 = Ar^(2/3)*(DeqCurr/2);
            b1 = Ar^(-1/3)*(DeqCurr/2);
            c1 = b1; 
            cubeOverlap = boxOverlap(cubeLength, a1, b1, c1, C1, R1);
            if cubeOverlap == true
                continue;
            end
            if currTotal == 0
                currTotal = currTotal + 1;
                a{currTotal} = a1;
                b{currTotal} = b1;
                c{currTotal} = c1;
                C{currTotal} = C1;
                R{currTotal} = R1;
                accepted = true;
            else
                spheroidOverlap = false;
                for i = 1:currTotal
                    par1 = struct('a', a1, 'b', b1, 'c', c1, 'C', C1, 'R', R1);
                    par2 = struct('a', a{i}, 'b', b{i}, 'c', c{i}, 'C', C{i}, 'R', R{i});
                    [~, ~, d] = movingBallsPair(par1, par2, 1e-5, 500);
                    if d == 0
                        spheroidOverlap = true;
                        break;
                    end
                end
                if spheroidOverlap == true
                    continue;
                else
                    currTotal = currTotal + 1;
                    a{currTotal} = a1;
                    b{currTotal} = b1;
                    c{currTotal} = c1;
                    C{currTotal} = C1;
                    R{currTotal} = R1;
                    accepted = true;
                end
            end
        end
    end
    
    %Make the final array of structs
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
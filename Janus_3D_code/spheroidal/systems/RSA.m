function spheroids = RSA(cubeLength, Ar, Deq, parameter)

%Randomized Sequential Adsorption: This algorithm generates spheroids within some sort of container. This is done by randomly generating one, then checking if it intersects another spheroid or the container.  

volumeFraction = parameter*(4/3)*pi*(Deq/2)^3/(cubeLength^3);
fprintf('Desired Volume Fraction:%f\n', volumeFraction);

a = cell(parameter, 1);
b = cell(parameter, 1);
c = cell(parameter, 1);
C = cell(parameter, 1);
R = cell(parameter, 1);

total = 0;
a1 = Ar^(2/3)*(Deq/2);
b1 = Ar^(-1/3)*(Deq/2);
c1 = b1; 


while total < parameter
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

%This uses the Jia et all modified overlap detection algorithm, but I am currently unable to get it to work so we will instead use the moving balls algorithm.
%{
function overlap = ellipsoidOverlap(a1, b1, c1, C1, R1, a2, b2, c2, C2, R2)

    overlap = true;
    %Get New Coordinates
    S1 = diag([a1^-2 b1^-2 c1^-2 -1]);
    M1 = [R1.' -R1.'*C1; zeros(1, size(R1,1)) 1];
    A1 = M1.'*S1*M1;

    S2 = diag([a2^-2 b2^-2 c2^-2 -1]);
    M2 = [R2.' -R2.'*C2; zeros(1, size(R2,1)) 1];
    A2 = M2.'*S2*M2;

    A1inv = A1.^-1; 
    p = poly(-A1inv*A2);
    
    var = 0;
    for i = 1:3
        if p(i)*p(i + 1) < 0
            var = var + 1;
        end
    end

    if var == 2
        bbar = -p(2)/4;
        cbar = p(3)/6;
        dbar = -p(4)/4;
        ebar = p(5);

        Delta2 = bbar^2 - cbar;
        Delta3 = cbar^2- bbar*dbar;
        W1 = dbar - dbar*cbar;
        W2 = bbar*ebar - cbar*dbar;
        W3 = ebar - bbar*dbar;

        T = -9*W1^2 + 27*Delta2*Delta3 - 3*W3*Delta2;
        A = W3 +3*Delta3;
        B = -dbar*W1 - ebar*Delta2 - cbar*Delta3;
        T2 = A*W1 - 3*bbar*B;
        Delta1 = A^3 - 27*B^2;

        sr22 = Delta2;
        sr20 = -W3;
        sr11 = T;
        sr10 = T2;
        sr0 = Delta1;

        if (sr22 > 0 && sr11 > 0 && sr0) || (sr22 > 0 && sr11 > 0 && sr10 > 0 && sr == 0)
            overlap = false;
        end
    end

end
%}
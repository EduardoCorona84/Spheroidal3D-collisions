function spheroidalPairValidation(Ar, Deq, xc, alg, tol, maxIter)
    %This function takes in spheroid parameters, a point xc, an algorithm and then constructs a distance problem between a pair of spheroids with a known solution. The computed solution from the given algorithm is compared to the analytic solution.

    %Set up structs with params for distance algorithm.
    a = Ar^(2/3)*(Deq/2);
    b = Ar^(-1/3)*(Deq/2);
    C = [xc 0 0].';
    r = sqrt(2)/2';
    R1 = [r -r 0; r r 0; 0 0 1];
    R2 = [r r 0; -r r 0; 0 0 1];
    par1 = struct('a', a, 'b', b, 'c', b, 'C', C, 'R', R1);
    par2 = struct('a', a, 'b', b, 'c', b, 'C', -C, 'R', R2);
    spheroids = [par1 par2];
    plotSpheroids(spheroids);

    %Find algorithm values
    distanceAlgo = str2func(alg);
    [x1 x2 d] = distanceAlgo(par1, par2, tol, maxIter);

    %Find true values using the support mapping
    n = [1 0 0].';
    x1true = ellipsoidSupportMapping(-n, a, b, b, C, R1);
    x2true = ellipsoidSupportMapping(n, a, b, b, -C, R2);
    dtrue = norm(x1true - x2true);

    %Compare values 
    fprintf('X1 Error:%f\n', norm(x1 - x1true));
    fprintf('X2 Error:%f\n', norm(x2 - x2true));
    fprintf('Distance Error:%f\n', norm(d - dtrue));




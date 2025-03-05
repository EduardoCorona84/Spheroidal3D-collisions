cubeLength = 4;
eta = 1;
maxIter = 500;

%Test for large aspect ratios 

Ar = 8;
Deq = 1;
total = 2;
tol = 1e-6;

rateLargeMB = zeros(10,1);
rateLargeGJKJohn = zeros(10,1);
rateLargeGJKJohnN = zeros(10,1);
rateLargeGJKSigned = zeros(10,1);
rateLargeGJKSignedN = zeros(10,1);

for i=1

    spheroids = RSA(cubeLength, Ar, Deq, total);

    [x1, x2, d, dIters] = movingBallsPair(spheroids(1), spheroids(2), tol, maxIter, false);
    iters = find(~dIters, 1, 'first') - 1; 
    disp(dIters(1:iters));
    xAxis = log(abs((dIters(6:iters - 4) - dIters(5:iters - 5))./(dIters(5:iters - 5) - dIters(4:iters - 6))));
    yAxis = log(abs((dIters(5:iters - 5) - dIters(4:iters - 6))./(dIters(4:iters - 6) - dIters(3:iters - 7))));
    p = polyfit(xAxis, yAxis, 1);
    rateLargeMB(i) = p(1);

    [x1, x2, d, dIters] = GJKJohnsonPair(spheroids(1), spheroids(2), tol, maxIter, false);
    iters = find(~dIters, 1, 'first') - 1; 
    disp(dIters(1:iters));
    xAxis = log(abs((dIters(6:iters - 4) - dIters(5:iters - 5))./(dIters(5:iters - 5) - dIters(4:iters - 6))));
    yAxis = log(abs((dIters(5:iters - 5) - dIters(4:iters - 6))./(dIters(4:iters - 6) - dIters(3:iters - 7))));
    p = polyfit(xAxis, yAxis, 1);
    rateLargeGJKJohn(i) = p(1);

    [x1, x2, d, dIters] = GJKJohnsonNestPair(spheroids(1), spheroids(2), tol, maxIter, false);
    iters = find(~dIters, 1, 'first') - 1; 
    disp(dIters(1:iters));
    xAxis = log(abs((dIters(6:iters - 4) - dIters(5:iters - 5))./(dIters(5:iters - 5) - dIters(4:iters - 6))));
    yAxis = log(abs((dIters(5:iters - 5) - dIters(4:iters - 6))./(dIters(4:iters - 6) - dIters(3:iters - 7))));
    p = polyfit(xAxis, yAxis, 1);
    rateLargeGJKJohnN(i) = p(1);

    [x1, x2, d, dIters] = GJKSignedVolumesPair(spheroids(1), spheroids(2), tol, maxIter, false);
    iters = find(~dIters, 1, 'first') - 1; 
    disp(dIters(1:iters));
    xAxis = log(abs((dIters(6:iters - 4) - dIters(5:iters - 5))./(dIters(5:iters - 5) - dIters(4:iters - 6))));
    yAxis = log(abs((dIters(5:iters - 5) - dIters(4:iters - 6))./(dIters(4:iters - 6) - dIters(3:iters - 7))));
    p = polyfit(xAxis, yAxis, 1);
    rateLargeGJKSigned(i) = p(1);

    [x1, x2, d, dIters] = GJKSignedVolumesNestPair(spheroids(1), spheroids(2), tol, maxIter, false);
    iters = find(~dIters, 1, 'first') - 1; 
    disp(dIters(1:iters));
    xAxis = log(abs((dIters(6:iters - 4) - dIters(5:iters - 5))./(dIters(5:iters - 5) - dIters(4:iters - 6))));
    yAxis = log(abs((dIters(5:iters - 5) - dIters(4:iters - 6))./(dIters(4:iters - 6) - dIters(3:iters - 7))));
    p = polyfit(xAxis, yAxis, 1);
    rateLargeGJKSignedN(i) = p(1);

end

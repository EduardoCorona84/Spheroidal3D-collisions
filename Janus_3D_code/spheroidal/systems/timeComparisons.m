cubeLength = 4;
eta = 0.1;
maxIter = 1000;

%Test for large aspect ratios 

Ar = 8;
Deq = 1;
total = 15;

timeLargeMBLow = zeros(10,1);
timeLargeGJKJohnLow = zeros(10,1);
timeLargeGJKJohnNLow = zeros(10,1);

timeLargeMBHigh = zeros(10,1);
timeLargeGJKJohnHigh = zeros(10,1);
timeLargeGJKJohnNHigh = zeros(10,1);

for i=1:10
    spheroids = RSA(cubeLength, Ar, Deq, total);
    
    tol = 1e-4;
    f = @() allWithPairs(spheroids, tol, maxIter, eta, 'mb');
    time = timeit(f);
    timeLargeMBLow(i) = time;
    f = @() allWithPairs(spheroids, tol, maxIter, eta, 'GJKJohn');
    time = timeit(f);
    timeLargeGJKJohnLow(i) = time;
    f = @() allWithPairs(spheroids, tol, maxIter, eta, 'GJKJohnN');
    time = timeit(f);
    timeLargeGJKJohnNLow(i) = time;

    tol = 1e-8;
    f = @() allWithPairs(spheroids, tol, maxIter, eta, 'mb');
    time = timeit(f);
    timeLargeMBHigh(i) = time;
    f = @() allWithPairs(spheroids, tol, maxIter, eta, 'GJKJohn');
    time = timeit(f);
    timeLargeGJKJohnHigh(i) = time;
    f = @() allWithPairs(spheroids, tol, maxIter, eta, 'GJKJohnN');
    time = timeit(f);
    timeLargeGJKJohnNHigh(i) = time;
end

boxchart([timeLargeMBLow timeLargeGJKJohnLow timeLargeGJKJohnNLow]);
xlabel('algorithm');
ylabel('time (s)');
saveas(gcf, 'LargeAspectRatioLowAccuracy.pdf');


boxchart([timeLargeMBHigh timeLargeGJKJohnHigh timeLargeGJKJohnNHigh]);
xlabel('algorithm');
ylabel('time (s)');
saveas(gcf, 'LargeAspectRatioHighAccuracy.pdf');


%Test for medium aspect rations

Ar = 2;
timeMediumMBLow = zeros(10,1);
timeMediumGJKJohnLow = zeros(10,1);
timeMediumGJKJohnNLow = zeros(10,1);

timeMediumMBHigh = zeros(10,1);
timeMediumGJKJohnHigh = zeros(10,1);
timeMediumGJKJohnNHigh = zeros(10,1);

for i=1:10
    spheroids = RSA(cubeLength, Ar, Deq, total);
    
    tol = 1e-4;
    f = @() allWithPairs(spheroids, tol, maxIter, eta, 'mb');
    time = timeit(f);
    timeMediumMBLow(i) = time;
    f = @() allWithPairs(spheroids, tol, maxIter, eta, 'GJKJohn');
    time = timeit(f);
    timeMediumGJKJohnLow(i) = time;
    f = @() allWithPairs(spheroids, tol, maxIter, eta, 'GJKJohnN');
    time = timeit(f);
    timeMediumGJKJohnNLow(i) = time;

    tol = 1e-8;
    f = @() allWithPairs(spheroids, tol, maxIter, eta, 'mb');
    time = timeit(f);
    timeMediumMBHigh(i) = time;
    f = @() allWithPairs(spheroids, tol, maxIter, eta, 'GJKJohn');
    time = timeit(f);
    timeMediumGJKJohnHigh(i) = time;
    f = @() allWithPairs(spheroids, tol, maxIter, eta, 'GJKJohnN');
    time = timeit(f);
    timeMediumGJKJohnNHigh(i) = time;
end

boxchart([timeMediumMBLow timeMediumGJKJohnLow timeMediumGJKJohnNLow]);
xlabel('algorithm');
ylabel('time (s)');
saveas(gcf, 'MediumAspectRatioLowAccuracy.pdf');


boxchart([timeMediumMBHigh timeMediumGJKJohnHigh timeMediumGJKJohnNHigh]);
xlabel('algorithm');
ylabel('time (s)');
saveas(gcf, 'MediumAspectRatioHighAccuracy.pdf');

%Test for small aspect ratios
Ar = 1/8;
timeSmallMBLow = zeros(10,1);
timeSmallGJKJohnLow = zeros(10,1);
timeSmallGJKJohnNLow = zeros(10,1);

timeSmallMBHigh = zeros(10,1);
timeSmallGJKJohnHigh = zeros(10,1);
timeSmallGJKJohnNHigh = zeros(10,1);

for i=1:10
    spheroids = RSA(cubeLength, Ar, Deq, total);
    
    tol = 1e-4;
    f = @() allWithPairs(spheroids, tol, maxIter, eta, 'mb');
    time = timeit(f);
    timeSmallMBLow(i) = time;
    f = @() allWithPairs(spheroids, tol, maxIter, eta, 'GJKJohn');
    time = timeit(f);
    timeSmallGJKJohnLow(i) = time;
    f = @() allWithPairs(spheroids, tol, maxIter, eta, 'GJKJohnN');
    time = timeit(f);
    timeSmallJKJohnNLow(i) = time;

    tol = 1e-8;
    f = @() allWithPairs(spheroids, tol, maxIter, eta, 'mb');
    time = timeit(f);
    timeSmallMBHigh(i) = time;
    f = @() allWithPairs(spheroids, tol, maxIter, eta, 'GJKJohn');
    time = timeit(f);
    timeSmallGJKJohnHigh(i) = time;
    f = @() allWithPairs(spheroids, tol, maxIter, eta, 'GJKJohnN');
    time = timeit(f);
    timeSmallGJKJohnNHigh(i) = time;
end

boxchart([timeSmallMBLow timeSmallGJKJohnLow timeSmallGJKJohnNLow]);
xlabel('algorithm');
ylabel('time (s)');
saveas(gcf, 'SmallAspectRatioLowAccuracy.pdf');


boxchart([timeSmallMBHigh timeSmallGJKJohnHigh timeSmallGJKJohnNHigh]);
xlabel('algorithm');
ylabel('time (s)');
saveas(gcf, 'SmallAspectRatioHighAccuracy.pdf');

%Test for heterogenous mixture of spheroid aspect ratios

Ar = [1/8 2 8];
timeHeteroMBLow = zeros(10,1);
timeHeteroGJKJohnLow = zeros(10,1);
timeHeteroGJKJohnNLow = zeros(10,1);

timeHeteroMBHigh = zeros(10,1);
timeHeteroGJKJohnHigh = zeros(10,1);
timeHeteroGJKJohnNHigh = zeros(10,1);

for i=1:10
    spheroids = RSAHeterogeneous(cubeLength, Ar, Deq, total);
    
    tol = 1e-4;
    f = @() allWithPairs(spheroids, tol, maxIter, eta, 'mb');
    time = timeit(f);
    timeHeteroMBLow(i) = time;
    f = @() allWithPairs(spheroids, tol, maxIter, eta, 'GJKJohn');
    time = timeit(f);
    timeHeteroGJKJohnLow(i) = time;
    f = @() allWithPairs(spheroids, tol, maxIter, eta, 'GJKJohnN');
    time = timeit(f);
    timeHeteroGJKJohnNLow(i) = time;

    tol = 1e-8;
    f = @() allWithPairs(spheroids, tol, maxIter, eta, 'mb');
    time = timeit(f);
    timeHeteroMBHigh(i) = time;
    f = @() allWithPairs(spheroids, tol, maxIter, eta, 'GJKJohn');
    time = timeit(f);
    timeHeteroGJKJohnHigh(i) = time;
    f = @() allWithPairs(spheroids, tol, maxIter, eta, 'GJKJohnN');
    time = timeit(f);
    timeHeteroGJKJohnNHigh(i) = time;
end

boxchart([timeHeteroMBLow timeHeteroGJKJohnLow timeHeteroGJKJohnNLow]);
xlabel('algorithm');
ylabel('time (s)');
saveas(gcf, 'HeteroAspectRatioLowAccuracy.pdf');


boxchart([timeHeteroMBHigh timeHeteroGJKJohnHigh timeHeteroGJKJohnNHigh]);
xlabel('algorithm');
ylabel('time (s)');
saveas(gcf, 'HeteroAspectRatioHighAccuracy.pdf');
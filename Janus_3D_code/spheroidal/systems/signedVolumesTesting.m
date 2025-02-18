Ar = 6;
Deq = 1;
tol = 1e-9;
maxIter = 1000;

a = Ar^(2/3)*(Deq/2);
b = Ar^(-1/3)*(Deq/2);
r = sqrt(2)/2';
R1 = [r -r 0; r r 0; 0 0 1];
R2 = [r r 0; -r r 0; 0 0 1];

%Overlapping Case, make sure all algorithms return 0 distance
xc = 0.5;
C = [xc 0 0].';
par1 = struct('a', a, 'b', b, 'c', b, 'C', C, 'R', R1);
par2 = struct('a', a, 'b', b, 'c', b, 'C', -C, 'R', R2);

%Overlapping Case, make sure all algorithms return 0 distance

[x1 x2 dMB] = movingBallsPair(par1, par2, tol, maxIter, true);
fprintf('d MB:%f\n', dMB);
[x1 x2 dSV] = GJKJohnsonNestPair(par1, par2, tol, maxIter, true);
fprintf('d SV:%f\n', dSV);



%Nonoverlaping case for various distances, compare to the actual solution and find an average computation time

Ar = 2;
xc = 3;
C = [xc 0 0].';
par1 = struct('a', a, 'b', b, 'c', b, 'C', C, 'R', R1);
par2 = struct('a', a, 'b', b, 'c', b, 'C', -C, 'R', R2);

n = [1 0 0].';
x1true = ellipsoidSupportMapping(-n, a, b, b, C, R1);
x2true = ellipsoidSupportMapping(n, a, b, b, -C, R2);
dtrue = norm(x1true - x2true);

[x1 x2 d] = movingBallsPair(par1, par2, tol, maxIter, true);
fprintf('MB X1 Error:%f\n', norm(x1 - x1true));
fprintf('MB X2 Error:%f\n', norm(x2 - x2true));

[x1 x2 d] = GJKJohnsonNestPair(par1, par2, tol, maxIter, true);
fprintf('GJK John N X1 Error:%f\n', norm(x1 - x1true));
fprintf('GJK John N X2 Error:%f\n', norm(x2 - x2true));
fprintf('GJK John N Distance Error:%f\n', norm(d - dtrue));

spheroids = RSA(10, 2, 1, 10);

[x1 x2 d n] = allWithPairs(spheroids, tol, 1000, 1, 'mb');

disp(d);

[x1 x2 d n] = allWithPairs(spheroids, tol, 1000, 1, 'GJKSignedN');

disp(d);




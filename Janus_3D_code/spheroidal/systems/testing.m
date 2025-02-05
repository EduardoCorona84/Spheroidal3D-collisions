Ar = 4;
Deq = 1;
xc = 0.5;

a = Ar^(2/3)*(Deq/2);
b = Ar^(-1/3)*(Deq/2);
C = [xc 0 0].';
r = sqrt(2)/2';
R1 = [r -r 0; r r 0; 0 0 1];
R2 = [r r 0; -r r 0; 0 0 1];
par1 = struct('a', a, 'b', b, 'c', b, 'C', C, 'R', R1);    par2 = struct('a', a, 'b', b, 'c', b, 'C', -C, 'R', R2);

spheroids = [par1 par2];
plotSpheroids(spheroids);

%Find algorithm values
[x1 x2 d] = GJKPair(par1, par2, 1e-5, 1000);
disp(d);
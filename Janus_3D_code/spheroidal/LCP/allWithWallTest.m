eta = 1;
%None of these should be detected as neighbors
spheroids = RSA(10, 2, 2, 4);

for i = 1:length(spheroids)
    spheroids(i).C = spheroids(i).C + [0 0 10].';
end
[x1, x2, distances, neighbors] = allWithWall(spheroids, eta);
disp(nnz(neighbors));
disp(distances);


%Only the third should be detected as a neighbor
a = cell(3,1);
b = cell(3,1);
c = cell(3,1);
C = cell(3,1);
R = cell(3,1);

for i=1:3
    a{i} = 1;
    b{i} = 4;
    c{i} = 1;
end 

C{1} = [0 0 10].';
C{2} = [1 3 5.5].';
C{3} = [-1 4 3.5].';

[R{1}, ~] = qr(randn(3,3));
R{2} = eye(3);
R{3} = eye(3);

spheroids = struct('a', a, 'b', b, 'c', c, 'C', C, 'R', R);
[x1, x2, distances, neighbors] = allWithWall(spheroids, eta);
disp(nnz(neighbors));
disp(distances);

%The third and the second should be detected as neighbors
r = sqrt(2)/2';
R{2} = [1 0 0; 0 r -r ;0 r r];
spheroids = struct('a', a, 'b', b, 'c', c, 'C', C, 'R', R);
[x1, x2, distances, neighbors] = allWithWall(spheroids, eta);
disp(nnz(neighbors));
disp(distances);
%Basic case, test for just a sphere and its distance to the xy plane.
a = cell(1,1);
b = cell(1,1);
c = cell(1,1);
a{1} = 1;
b{1} = 1;
c{1} = 1;

C = [0 0 10].';
[R, ~] = qr(randn(3,3));

rotConfig = struct('a', a, 'b', b, 'c', c,'C', C, 'R', R);

planeNormal = [0 0 1].';
planePoint = [0 0 0].';

[x1, x2, d] = planeEllipsoidIntersection(rotConfig, planeNormal, planePoint);
fprintf('Sphere with center [0 0 10] with radius 1\n');
fprintf('Distance to plane: %f\n', d);
fprintf('x1 value: %f\n', x1);
fprintf('x2 value: %f\n', x2);

C = [1 3 10].';
rotConfig = struct('a', a, 'b', b, 'c', c,'C', C, 'R', R);
[x1, x2, d] = planeEllipsoidIntersection(rotConfig, planeNormal, planePoint);
fprintf('Sphere with center [1 3 10] with radius 1\n');
fprintf('Distance to plane: %f\n', d);
fprintf('x1 value: %f\n', x1);
fprintf('x2 value: %f\n', x2);

a{1} = 2;
b{1} = 2;
c{1} = 2;

C = [-1 7 8].';
rotConfig = struct('a', a, 'b', b, 'c', c,'C', C, 'R', R);
[x1, x2, d] = planeEllipsoidIntersection(rotConfig, planeNormal, planePoint);
fprintf('Sphere with center [-1 7 8] with radius 2\n');
fprintf('Distance to plane: %f\n', d);
fprintf('x1 value: %f\n', x1);
fprintf('x2 value: %f\n', x2);

%Test for an ellipsoid and its distance to the xy plane.
a{1} = 1;
b{1} = 1;
c{1} = 1;

%ellipsoid rotated  45 degrees
a{1} = 1;
b{1} = 2;
c{1} = 1;
C = [1 5 10].';
r = sqrt(2)/2';
R = [1 0 0; 0 r -r ;0 r r];
rotConfig = struct('a', a, 'b', b, 'c', c,'C', C, 'R', R);
plotSpheroids(rotConfig, 12);
[x1, x2, d] = planeEllipsoidIntersection(rotConfig, planeNormal, planePoint);
fprintf('Sphere with center [-1 7 8] with radius 2\n');
fprintf('Distance to plane: %f\n', d);
fprintf('x1 value: %f\n', x1);
fprintf('x2 value: %f\n', x2);


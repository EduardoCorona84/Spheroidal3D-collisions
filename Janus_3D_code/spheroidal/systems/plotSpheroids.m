function plotSpheroids(spheroids, upper)

    n = length(spheroids);
    for i = 1:n
        [X, Y, Z] = ellipsoid(0, 0, 0, spheroids(i).a, spheroids(i).b, spheroids(i).c);
        s = surf(X,Y,Z);
        h = hgtransform;
        set(s, "Parent", h);
        T = makehgtform("translate", spheroids(i).C);
        R = zeros(4,4);
        R(1:3, 1:3) = spheroids(i).R;
        R(4,4) = 1;
        A = T*R;
        h.Matrix = A;
        axis([-5 5 -5 5 0 upper]);
        hold on;
    end
end
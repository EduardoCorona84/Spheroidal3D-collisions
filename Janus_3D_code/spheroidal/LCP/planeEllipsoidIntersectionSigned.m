function [x1, x2, d] = planeEllipsoidIntersectionSigned(par, planeNormal, planePoint)
    
    %Find the rotation matrix for the ellipsoid
    R = par.R;
    
    %Find the center of the ellipsoid
    C = par.C;
    
    %Find the semi-axes of the ellipsoid
    a = par.a;
    b = par.b;
    c = par.c;
    
    %This closed solution is given by lagrange multipliers, see my notes.
    temp = R*diag([a b c].^2)*R.'*planeNormal;
    x1 = C - sqrt(dot(planeNormal, temp))^(-1).*temp;
    x2 = x1;
    x2(3) = 0;
    %find intersection point of the plane and the line
    %t = dot(planeNormal, planePoint - x1)/dot(planeNormal, planeNormal);
    %x2 = x1 +t*planeNormal;

    d = dot(planeNormal, x1 - planePoint);
    %d = max([d, 0]);


end
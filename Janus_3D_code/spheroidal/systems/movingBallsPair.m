function [x1, x2, d] = movingBallsPair(par1, par2, tol, maxIter, val)

%get parameters of the ellipsoids.
%for notation:
%capital C refers to the center of the ellipsoid
%little a,b,c are the semi axes of the ellipsoid
%capital R denotes rotation matrices
%S denotes center of spheres in the moving balls algorithm
C1 = par1.C; R1 = par1.R; a1 = par1.a; b1=par1.b; c1=par1.c; 
C2 = par2.C; R2 = par2.R; a2 = par2.a; b2=par2.b; c2=par2.c;

%check if row vectors
if isequal(size(C1),[1 3])
    C1 = C1.';
    C2 = C2.';
end

%make relevant SPD matrices for ellipsoids
D1 = diag([a1 b1 c1].^-2);
D2 = diag([a2 b2 c2].^-2);

A1 = R1*D1*R1.';
A2 = R2*D2*R2.';

%find gammas and set tolerance
smallestSemi = min([a1 b1 c1]);
largestSemi = max([a1 b1 c1]);
maxCurv1 = max([smallestSemi^2/largestSemi largestSemi^2/smallestSemi]);
gamma1 = smallestSemi^2;

smallestSemi = min([a2 b2 c2]);
largestSemi = max([a2 b2 c2]);
maxCurv2 = max([smallestSemi^2/largestSemi largestSemi^2/smallestSemi]);
gamma2 = smallestSemi^2;

eps_theta = sqrt((2*tol)/(maxCurv1 + maxCurv2));

%Find initial points using the centers of the ellipsoids and compute distances
iter = 1;
[x1, x2, d] = findIntersection(C1, C2, A1, A2, C1, C2);
if d == 0
    return;
end
bothTol = withinTol(x1 - x2, A1*(x1 - C1), eps_theta) || withinTol(x2 - x1, A1*(x2 - C2), eps_theta);
iter = iter + 1;

%Begin iterations >= 2, this time with centers of spheres determined by the update rule
while (~bothTol && iter < maxIter)
    S1 = x1 - gamma1.*(A1*(x1 - C1));
    S2 = x2 - gamma2.*(A2*(x2 - C2));
    
    [x1New, x2New, dNew] = findIntersection(S1, S2, A1, A2, C1, C2);
    if (dNew == 0)
        x1 = x1New;
        x2 = x2New;
        d = dNew;
        return;
    else
        x1 = x1New;
        x2 = x2New;
        d = dNew;
        bothTol = withinTol(x1 - x2, A1*(x1 - C1), eps_theta) || withinTol(x2 - x1, A1*(x2 - C2), eps_theta);
        iter = iter + 1;
    end

end

end

%Find points of intersection of line segment connecting the sphere (S1, S2) centers and determine if the ellipsoids overlap.
function [x1, x2, d] = findIntersection(S1, S2, A1, A2, C1, C2)
v = S2 - S1;

a = v.'*A1*v;
b = 2.*(S1 - C1).'*A1*v;
c = (S1 - C1).'*A1*(S1 - C1) - 1;

t1 = (-b + sqrt(b^2 - 4*a*c))/(2*a);

a = v.'*A2*v;
b = 2.*(S1 - C2).'*A2*v;
c = (S1 - C2).'*A2*(S1 - C2) - 1;

t2 = (-b - sqrt(b^2 - 4*a*c))/(2*a);

x1 = S1 + t1.*v;
x2 = S1 + t2.*v;

if t1 >= t2
    d = 0;
else
    d = norm(x1 - x2);

end

end

%Use the alternative check from Girault 2022 which avoids square root and arccos computation.
function withinTol = withinTol(x1, x2, eps)

testQuantity1 = dot(x1,x2)^2;
testQuantity2 = (1 - (eps^2)/2)^2*dot(x1,x1)*dot(x2,x2);
if testQuantity1 >= testQuantity2 
    withinTol = true;
else
    withinTol = false;
end

end
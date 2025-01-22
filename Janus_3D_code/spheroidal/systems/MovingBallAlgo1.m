function [x1,x2,dst,x1l,x2l,dstl] = MovingBallAlgo1(par1,par2,tol,maxit)

% Get spheroid params
C1 = par1.C; R1 = par1.R; ax1 = par1.a; by1=par1.b; cz1=par1.c; 
C2 = par2.C; R2 = par2.R; ax2 = par2.a; by2=par2.b; cz2=par2.c;

% Diagonal matrices, Am matrices and qm. 
% Ellipsoids are defined by qm(x) <= 1, and surface is qm(x)=1. 
D1 = diag([ax1 by1 cz1].^(-2)); A1 = R1*D1*R1.'; %q1 = @(x) (x-C1)*(A1*(x-C1)');
D2 = diag([ax2 by2 cz2].^(-2)); A2 = R2*D2*R2.'; %q2 = @(x) (x-C2)*(A2*(x-C2)');

%Set eps_th based on max curvature radii and target tol
rmx1 = max([ax1 by1 cz1]); rmn1 = min([ax1 by1 cz1]);
Rmx1 = max(rmx1^2/rmn1, rmn1^2/rmx1);
gamma1 = rmn1^2;

rmx2 = max([ax2 by2 cz2]); rmn2 = min([ax2 by2 cz2]);
Rmx2 = max(rmx2^2/rmn2, rmn2^2/rmx2);
gamma2 = rmn2^2;

eps_th = sqrt((2*tol)/(Rmx1 + Rmx2)); 

% Pick initial guess based on centers 
[x1,x2,dst]=LOCAL_FindInters(C1,C2,C1,C2,A1,A2);
x1l = x1; x2l = x2; dstl = dst; 

iter=1; 

% Start iteration 
while ((costest(x2-x1,A1*x1(:) - A1'*(C1'),eps_th) || costest(x1-x2,A2*x2(:) - A2'*(C2'),eps_th)) && iter<maxit)

    ck1 = x1 - gamma1*(A1*(x1') - A1'*(C1'))';
    ck2 = x2 - gamma2*(A2*(x2') - A2'*(C2'))';
    [x1,x2,dst]=LOCAL_FindInters(ck1,ck2,C1,C2,A1,A2);
    
    display([iter,dst])

    if dst==0
        break; 
    else
        x1l(iter+1,:) = x1; 
        x2l(iter+1,:) = x2; 
        dstl(iter+1) = dst; 
    end

    iter=iter+1; 

end

end

function [x1n,x2n,dst]=LOCAL_FindInters(x1,x2,C1,C2,A1,A2)

d = x2-x1; 
a1 = 0.5*d*(A1*d');
b1 = d*(A1*(x1-C1)');
c1 = 0.5*((x1-C1)*(A1*(x1-C1)') - 1);
t1 = (-b1+sqrt(b1^2-4*a1*c1))/(2*a1);

a2 = 0.5*d*(A2*d');
b2 = d*(A2*(x1-C2)');
c2 = 0.5*((x1-C2)*(A2*(x1-C2)') - 1);
t2 = (-b2-sqrt(b2^2-4*a2*c2))/(2*a2);

x1n = x1 + t1*d;
x2n = x1 + t2*d; 

if t2<t1
    dst=0; 
else
    dst=norm(x2n-x1n); 
end 

end

function tst = costest(v1,v2,eps_th)
cs12 = v1(:)'*v2(:); 
tst = cs12 < (1 - 0.5*eps_th^2)*sqrt((v1(:)'*v1(:))*(v2(:)'*v2(:))); 

end
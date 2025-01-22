function [x1,x2,dst] = MovingBallNeighbors(par1,par2,tol,maxit)
% Given a parameter struct for spheroid 1 and struct array par2 for neighboring 
% spheroids, this code returns pairs of contact points x1, x2 and distances

nn = length(par2); %number of neighbors 
R2 = cell(nn,1); D2=R2; A2=D2;
C2 = zeros(nn,3); 
ax2=zeros(nn,1); by2=ax2; cz2=ax2; Rmx2=ax2; gamma2=ax2; 
x1 = zeros(nn,3); x2=zeros(nn,3); dst=zeros(nn,1); g2=x2; 

% Get spheroid params and define objects for source spheroid 1 
C1 = par1.C; R1 = par1.R; ax1 = par1.a; by1=par1.b; cz1=par1.c; 
% Diagonal matrices, Am matrices and qm. 
D1 = diag([ax1 by1 cz1].^(-2)); A1 = R1*D1*R1.'; %q1 = @(x) (x-C1)*(A1*(x-C1)');
%Set eps_th based on max curvature radii and target tol
rmx1 = max([ax1 by1 cz1]); rmn1 = min([ax1 by1 cz1]);
Rmx1 = max(rmx1^2/rmn1, rmn1^2/rmx1);
gamma1 = rmn1^2;

for i=1:nn
    C2(i,:) = par2(i).C; R2{i} = par2(i).R; 
    ax2(i) = par2(i).a; by2(i)=par2(i).b; cz2(i)=par2(i).c;
    % Diagonal matrices, Am matrices and qm. 
    D2{i} = diag([ax2(i) by2(i) cz2(i)].^(-2)); 
    A2{i} = R2{i}*D2{i}*R2{i}.'; 

    rmx2 = max([ax2(i) by2(i) cz2(i)]); rmn2 = min([ax2(i) by2(i) cz2(i)]);
    Rmx2(i) = max(rmx2^2/rmn2, rmn2^2/rmx2);
    gamma2(i) = rmn2^2;

    % Pick initial guess based on centers 
    [x1(i,:),x2(i,:),dst(i)]=LOCAL_FindInters(C1,C2(i,:),C1,C2(i,:),A1,A2{i});

    g2(i,:) = (A2{i}*((x2(i,:)-C2(i,:))'))';
end

eps_th = sqrt((2*tol)./(Rmx1 + Rmx2)); 
v1 = x2-x1; g1 = (A1*(x1-repmat(C1,nn,1))')';
cstest1 = costest(v1,g1,eps_th);

v2 = x1-x2; 
cstest2 = costest(v2,g2,eps_th);

criteria = cstest1 | cstest2; 
scrit = sum(criteria)>0; 

iter=1; 

% Start iteration 
while (scrit && iter<maxit)

    ck1 = x1 - gamma1*(A1*(x1') - repmat(A1'*(C1'),1,nn))';
    ck2=x2; 

    for i=1:nn
        if criteria(i)
            ck2(i,:) = x2(i,:) - gamma2(i)*(A2{i}*(x2(i,:)') - A2{i}'*(C2(i,:)'))';
            [x1(i,:),x2(i,:),dst(i)]=LOCAL_FindInters(ck1(i,:),ck2(i,:),C1,C2(i,:),A1,A2{i});
            g2(i,:) = (A2{i}*((x2(i,:)-C2(i,:))'))';
        end
    end
    
    %display([iter sum(criteria)/nn min(dst) mean(dst) max(dst)])

    v1 = x2-x1; g1 = (A1*(x1-repmat(C1,nn,1))')';
    cstest1 = costest(v1,g1,eps_th);

    v2 = x1-x2; %g2 = (A2*((x2-C2)'))';
    cstest2 = costest(v2,g2,eps_th);

    criteria = (cstest1 | cstest2) & dst>0; 
    scrit = sum(criteria)>0; 

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
cs12 = sum(v1.*v2,2); 
n1 = sum(v1.*v1,2); n2 = sum(v2.*v2,2);
tst = cs12 < (1 - 0.5*eps_th.^2).*sqrt(n1.*n2); 

end
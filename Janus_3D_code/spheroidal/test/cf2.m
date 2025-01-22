function [H,iters]=cf2(n,m,x)
%Compute the continued fraction that relates Q_n^m with Q_{n-1}^m, to
%desired tolerance, using Modified Lentz's Method.

tol=1e-12;
max_iters=100000;

tiny=1e-300;
f0=tiny;
c0=f0;
d0=0;

%sprintf('Initial values:')
%sprintf("f_0=%.6e",f0)
%sprintf("c_0=%.6e",c0)
%sprintf("d_0=%.6e",d0)

for k=1:max_iters
    %sprintf("iteration %d",k-n+1)

    d1=b(n,k,m,x)+a(n,k,m).*d0;
    if d1==0
        d1=tiny;
    end

    c1=b(n,k,m,x)+a(n,k,m)./c0;
    if c1==0
        c1=tiny;
    end

%   sprintf("c_%d=%.6e",k-n+1,c1)
%   sprintf("d_%d=%.6e",k-n+1,d1)

    d1=1./d1;
    delta=c1.*d1;
    f1=delta.*f0;

%   sprintf("f_%d=%.6e",k-n+1,f1)
%   sprintf("abs(1-delta)=%.6e",abs(1-delta))

    iters=k;

    if abs(delta-1)<tol
        break
    end

    if k==max_iters
        error("Continued fraction algorithm did not converge.")
    end

    %update values for recursion
    c0=c1;
    d0=d1;
    f0=f1;
   
end
H=f1;
end

function bj=b(n,j,m,x)
%b_j(n,m,x) coefficient in recursion relation
% bj= (2*(n+j-1)+1).*x ./(n+j-1+m);

bj=-(2*(n+j-1)+1).*x ./(n+j-m);
end

function aj=a(n,j,m)
%a_j coefficient in recursion relation
% if j==1
%     aj=1;
% else
%     aj=-(n+j-1-m)./(n+j-2+m);
% end

aj=-(n+j-1+m)./(n+j-m);
end


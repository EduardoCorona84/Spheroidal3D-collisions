function [y,z] = ALF(n,x,debug)
%computes the associated Legendre functions 
%of degree N and order M = 0, 1, ..., N, evaluated for each element
%of X.  N must be a scalar integer and X must contain real values greater
%than 1.

if nargin==2
    %optional input, if deubg=1 it will print intermediate results.
    debug=0;
end

x=reshape(x,1,length(x));
y=zeros([n+1,length(x)]);
z=zeros([n+1,length(x)]);

%Use forward recursion to compute P_n^m for m=0,...,n
for m=0:n

    %calculate initial legendre functions P_m^m and P_{m+1}^m
    P00=1;
    Pmm=P00;  %P_m^m=1 if m=0
    if m>0
      Pmm=prod(1:2:2*m-1).*(x.^2-1).^(m/2);   %P_m^m(x)=(2m-1)!! (x^2-1)^(m/2) if m>0
    elseif m<0
      error("Not valid m. Function only accepts nonnegative values of m.")
    end

    %Identity: P_{m+1}^m = (2m+1)x P_m^m(x) 
    Pmp1m = (2*m+1).* x.* Pmm;
    
    %testing/deubg
    if debug
    sprintf("this m=%d",m)
    sprintf("Computed P_%d^%d(x)",[m,m])
    if length(x)<5
    sprintf("P_%d^%d(x)= ",[m,m])
    disp(Pmm)
    end
    end    

    %testing/deubg
    if debug
    sprintf("computed P_%d^%d(x)",[m+1,m])
    if length(x)<5       
    sprintf("P_%d^%d(x)= ",[m+1,m])
    disp(Pmp1m)
    end
    end

    %no recurrence is needed for these cases (n=m or n=m+1) so we treat
    %them separately.
    if n <= m+1    
        H=cf(m+1,m,x);
        Qmm=factorial(2*m).* (-1)^m ./(Pmp1m-H.*Pmm);
        Qmp1m=H.* Qmm;
        if n==m
            Pnm=Pmm;
            Qnm=Qmm;
        elseif n==m+1
            Pnm=Pmp1m;
            Qnm=Qmp1m;
        end

        
    else

      Plm=Pmp1m;       %P_l^m
      Plminus1m=Pmm;    %P_{l-1}^m
      for l=m+1:n-1

         %Given P_l^m and P_{l-1}^m, use recuscursion to calculate
         %P_{l+1}^m up until P_n^m.
         Plp1m= (2*l+1)/(l-m+1) .* x .* Plm - (l+m)/(l-m+1) .* Plminus1m;

         %update for the next iteration
         Plminus1m=Plm;
         Plm=Plp1m;

      end
      Pnm=Plm;             %actual P_n^m we want
      Pnminus1m=Plminus1m;  %We will need this to compute Q_n^m
      
      %calculate Q_n^m and Q_{n-1}^m
      H=cf(n,m,x);
      Qnminus1m=factorial(n+m-1)/factorial(n-m) .* (-1)^m ./(Pnm-H.*Pnminus1m);
      Qnm=H.* Qnminus1m;
    
    end
    
    %store into y matrix
    y(m+1,:)=Pnm;
    z(m+1,:)=Qnm;
    
    if debug
    sprintf("P_%d^%d has been found:",[n,m])
    if length(x)<5
    disp(Pnm)
    end
    sprintf("Q_%d^%d has been found:",[n,m])
    if length(x)<5
    disp(Qnm)
    end
    end

    %Use backward recusrion to compute Q_n^m for m=n,n-1,...,0
    %To do later. Only needed if we want to provide additional Q's beyond
    %just Qnm.
    
    end
end

function [H,iters]=cf(n,m,x)
%Compute the continued fraction that relates Q_n^m with Q_{n-1}^m, to
%desired tolerance, using Modified Lentz's Method.

tol=1e-10;
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


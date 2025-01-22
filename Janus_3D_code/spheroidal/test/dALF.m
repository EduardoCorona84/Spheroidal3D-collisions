function [dy,dz]=dALF(n,m,x,mode,test)
%compute the derivative of P_n^m(x) and Q_n^m(x) with respect to x, x>1

if nargin==3
    test=false;
    mode=0;
elseif nargin==4;
    test=false;
end

[pnp1,qnp1]=ALF(n+1,x);
[pn,qn]=ALF(n,x);

%coefficient to multiply for negative m values
mc0=1;
mc1=1;
mcm1=1;
if m<0
    mc0=factorial(n+m)./factorial(n-m);
    mc1=factorial(n+1+m)./factorial(n+1-m);
    if m>-n
        mcm1=factorial(n-1+m)./factorial(n-1-m);
    end
end
am=abs(m);

% %Strategy (A)
% %----------------------------------------------------
% % recursion formula: (1-x^2)dP_n^m/dx = (m-n-1)P_{n+1}^m + (n+1)xP_n^m
% dy=((m-n-1).*mc1.*pnp1(am+1,:)+(n+1).*x.*mc0.*pn(am+1,:))./(1-x.^2);    %dP_n^m/dx
% %Plug dP_n^m/dx into Wronskian, then solve for dQ_n^m/dx
% dz=( factorial(n+m)./factorial(n-m).*(-1).^m./(1-x.^2) + dy.*mc0.*qn(am+1,:) )./(mc0.*pn(am+1,:));      %dQ_n^m/dx
% %----------------------------------------------------

%Strategy (B)
%----------------------------------------------------
% % recursion formula: (1-x^2)dP_n^m/dx = (m+n)P_{n-1}^m - nxP_n^m
% dy_numerator=-n.*x.*mc0.*pn(am+1,:);
% if n ~= abs(m)
%    [pnm1,qnm1]=ALF(n-1,x);
%    dy_numerator= dy_numerator + (m+n).*mcm1.*pnm1(am+1,:) ;
% end
% dy=dy_numerator./(1-x.^2);
% %Plug dP_n^m/dx into Wronskian, then solve for dQ_n^m/dx
% dz=( factorial(n+m)./factorial(n-m).*(-1).^m./(1-x.^2) + dy.*mc0.*qn(am+1,:) )./(mc0.*pn(am+1,:));      %dQ_n^m/dx
%----------------------------------------------------



%recursion formula: (1-x^2)dP_n^m/dx = (m-n-1)P_{n+1}^m + (n+1)xP_n^m
dy=((m-n-1).*mc1.*pnp1(am+1,:)+(n+1).*x.*mc0.*pn(am+1,:))./(1-x.^2);    %dP_n^m/dx
dz=((m-n-1).*mc1.*qnp1(am+1,:)+(n+1).*x.*mc0.*qn(am+1,:))./(1-x.^2);    %dQ_n^m/dx

%scaled by (n-m)!/(n+m)!
if mode==1
if m<0
    dy=( -(n+m+1).*pnp1(am+1,:) +(n+1).*x.*pn(am+1,:) )./(1-x.^2);
    dz=( -(n+m+1).*qnp1(am+1,:) +(n+1).*x.*qn(am+1,:) )./(1-x.^2); 
else
    dy=factorial(n-m)./factorial(n+m).*((m-n-1).*pnp1(am+1,:)+(n+1).*x.*pn(am+1,:))./(1-x.^2);
    dz=factorial(n-m)./factorial(n+m).*((m-n-1).*qnp1(am+1,:)+(n+1).*x.*qn(am+1,:))./(1-x.^2); 
end
end

%recursion formula: (1-x^2)dP_n^m/dx = (m+n)P_{n-1}^m - nxP_n^m
% if n==abs(m)
%    dy_numerator=-n.*x.*mc0.*pn(am+1,:);
%    dz_numerator=-n.*x.*mc0.*qn(am+1,:) + (n==0); %equals 1 if n=0, or first term if n>0
% else
%    [pnm1,qnm1]=ALF(n-1,x);
%    dy_numerator= (m+n).*mcm1.*pnm1(am+1,:) -n.*x.*mc0.*pn(am+1,:) ;
%    dz_numerator= (m+n).*mcm1.*qnm1(am+1,:) -n.*x.*mc0.*qn(am+1,:) ;
% end
% dy=dy_numerator./(1-x.^2);
% dz=dz_numerator./(1-x.^2);



if test

%test recursion formula: P_n^m dQ_n^m/dx - dP_n^m/dx Q_n^m = (n+m)!/(n-m)! * (-1)^m / (1-x^2)
Err=abs( mc0.*pn(am+1,:).*dz - mc0.*qn(am+1,:).*dy - factorial(n+m).*(-1).^m./(factorial(n-m)*(1-x.^2)) );

disp(sum(Err==0))

figure;
plot(x,log10(Err))
title("Wronskian error")


if m<0
    dytestfac1=( -(n+m+1).*pnp1(am+1,:) +(n+1).*x.*pn(am+1,:) )./(1-x.^2);
    dytestfac2= factorial(n-m)./factorial(n+m).*dy;
    dztestfac1=( -(n+m+1).*qnp1(am+1,:) +(n+1).*x.*qn(am+1,:) )./(1-x.^2);
    dztestfac2= factorial(n-m)./factorial(n+m).*dz;

    figure;
    plot(x,log10(abs(dytestfac2-dytestfac1)))
    title("dP test: error")

    figure;
    plot(x,log10(abs(dztestfac2-dztestfac1)))
    title("dQ test: error")

    figure;
    plot(x,log10(abs(dytestfac1+dztestfac1)))
    title("dP test: scale of fac*(dP+dQ)")

end

end
end
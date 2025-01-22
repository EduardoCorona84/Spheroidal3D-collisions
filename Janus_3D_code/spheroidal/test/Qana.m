function Qnm=Qana(n,m,x)
%Qnm in terms of hypergeometric functions
F=hypergeom([(n+m+2)/2,(n+m+1)/2],n+3/2,1./x.^2);
Qnm=2^(-n-1).*(-1).^m .* sqrt(pi).*factorial(n+m)./gamma(n+3/2).*x.^(-n-m-1) .*(x.^2-1).^(m/2).*F;
end
p=16; mu=2/sqrt(3);

f=func(0,0,p);
SLf=spheroidalSL(f,mu);
SLfsh=spheroidalAna(SLf);
[d1,d2]=size(SLfsh);

[u,v]=gl_grid(p);
fac=(1./sqrt(mu^2-cos(u).^2)); 
Y2=fac.*f;

G=Gmatrix(p,mu);
Gf=G*spheroidalAna(f);

% ss=sum((Y2-Gf).^2);
% disp(ss)

eigval=zeros(d1,1);
for i=1:d1
    n=floor(sqrt(i-1));
    m=i-1-n*(n+1);
    [P,Q]=ALF(n,mu);      
    Pnm=(factorial(n+m)./factorial(n-m)*(m<0)+(m>=0)).*P(abs(m)+1,:);
    Qnm=(factorial(n+m)./factorial(n-m)*(m<0)+(m>=0)).*Q(abs(m)+1,:);
    eigval(i)=factorial(n-m)./factorial(n+m).*(-1).^(m+1) .*sqrt(mu.^2-1).*Pnm.*Qnm;
end

err=abs(eigval-SLfsh);
figure;
plot(1:d1,log10(err))

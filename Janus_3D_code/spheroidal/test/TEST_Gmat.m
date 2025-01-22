clear;
u0=2/sqrt(3); a=1/sqrt(1+u0^2); p=8;
sp=(p+1)^2;
% ii = (1:sp)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;
% 
[u,v]=gl_grid(p);
Y=zeros(size(u,1),sp);
Y2=zeros(size(u,1),sp);
for n=0:p
    Yn=Ynm(n,[],u,v);
    fac=1./sqrt(u0^2+cos(u).^2); Ynfac=Yn.*fac;
    Y(:,n^2+1:(n+1)^2)=Yn;
    Y2(:,n^2+1:(n+1)^2)=Ynfac;
end
this_G = Gmatrix(p,u0,0,1);
shcfac=shAna(Y2);
Gshcfac=this_G\shcfac;
display(Gshcfac);
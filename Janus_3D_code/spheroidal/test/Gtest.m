p=4;
mu=2/sqrt(3);

nsh=(p+1).^2; %number of spherical harmonics in expansion for p
ii = (1:nsh)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;

G=zeros(nsh);

for i=1:nsh
    ni= nn(i); mi = mm(i); 
    sprintf("Projection: Y00 -> Y_%d^%d: %.5f", ni,mi, project(0,0,ni,mi,nump,mu))
    for j=1:nsh
        nj= nn(j); mj = mm(j); 
        G(i,j)=project(ni,mi,nj,mj,50,mu);
    end
end

close all
figure;
imagesc(G)
colorbar

%%

n=2;
m=1;
u=pi/3;
v=pi/5;

Y=Ynm(n,m,u,v);

Pn=legendre2(n,cos(u));
Pnm=Pn(abs(m)+1,:);
Ytest=Pnm.*exp(1i.*m.*v)*(-1).^n.*sqrt((2*n+1)/(4*pi)*factorial(n-m)/factorial(n+m));

sprintf("Y=%.5f + (%.5f) i",real(Y),imag(Y))
sprintf("Ytest=%.5f + (%.5f) i",real(Ytest),imag(Ytest))

%%

for nump=40:50
    err=log10(abs(project(2,0,2,0,nump,mu)-project(2,0,2,0,nump-1,mu)));
    sprintf("num points = %d, err=%.4f", nump,err)
end

%%



function proj=project(ni,mi,nj,mj,unp,mu)
% Compute projection of Y_{ni}^{mi}(u,v) onto new basis element
% Y_{nj}^{mj}(u,v)/sqrt(mu^2-u^2):
% 
%    <Y_{ni}^{mi}(u,v),Y_{nj}^{mj}(u,v)/sqrt(mu^2-u^2)> 
%    --------------------------------------------------    on the spheroid. 
%         ||Y_{nj}^{mj}(u,v)/sqrt(mu^2-u^2)||^2
% 
% *Notation: mu= 1/eccentricity of generating ellipse, 
% Ynm(u,v)=Cnm*P_n^m(u)*e^{imv}, -1<=u<=1, and 0<=v<=2pi, 
% Cnm normalization constant.
% unp = number of u points for GL quadrature

if mi ~= mj
    proj=0;
    return
else
    m=mi;
    mci=((-1).^m.*factorial(ni+m)./factorial(ni-m)*(m<0)+(m>=0)); %coefficient for negative m
    mcj=((-1).^m.*factorial(nj+m)./factorial(nj-m)*(m<0)+(m>=0)); %coefficient for negative m

    [up,uw]=g_grid(unp); %Gauss-Legendre Nodes and weights
    
    Pni=legendre2(ni,up);
    Pnmi=mci.*Pni(abs(m)+1,:);
    Pnj=legendre2(nj,up);
    Pnmj=mcj.*Pnj(abs(m)+1,:);
    
    %inner product
    integrand1=Pnmi.*Pnmj ./ ((mu.^2-up.^2).^(3/2));
    integral1=sum(uw.*integrand1);

    %L2 norm in Ynm/sqrt(mu^2-u^2) basis
    integrand2=(Pnmj./(mu.^2-up.^2)).^2;
    integral2=sum(uw.*integrand2);

    nci=(-1).^ni.*sqrt((2*ni+1) .* factorial(ni-m)./factorial(ni+m) );
    ncj=(-1).^nj.*sqrt((2*nj+1) .* factorial(nj-m)./factorial(nj+m) );
    proj=nci./ncj .* integral1 ./ integral2;
end
end



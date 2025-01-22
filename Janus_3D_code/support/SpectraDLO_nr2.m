function [ DLO ] =SpectraDLO_nr2( n,lambda,r )

%{
This function computes the solid harmonic coefficients for the Modified
(Screened) Laplace equation Delta*u - lambda*u = 0, lambda>0.
 
Inputs
n      (int) spherical harmonics degree
lambda (double) positive real value for lambda
r      (double) target radius

%}
mxlr=4;
lr=lambda*r;
dm=n;
if nargin<4
    out=1; 
end

p = (0:n)+0.5;
ilambda = (1i)*lambda;
J_pos = besselj(p,ilambda);
Jr_pos = besselj(p,ilambda*r);
J_neg = besselj(-p,ilambda);
Jr_neg = besselj(-p,ilambda*r);
MJ_first_pos = (1i.^(-p)).*J_pos;
MJr_first_pos = (1i.^(-p)).*Jr_pos;
MJ_first_neg = (1i.^p).*J_neg;
MJr_first_neg =(1i.^p).*Jr_neg; 


MSJ_first = sqrt(pi/(2*lambda))*MJ_first_pos;
der_MSJ_first = der_MSJ_f(MSJ_first, lambda);

if lambda*r<mxlr
    MJr_second = (pi./(2*sin(p*pi))).*(MJr_first_neg - MJr_first_pos);
    MSJr_second = sqrt(pi/(2*lambda*r))*MJr_second;
    DLO = -real((-(2*lambda^2/pi)*der_MSJ_first .*MSJr_second));
else
   
    LR = repmat((1/lr).^(0:dm),n+1,1);
    ank = zeros(n+1,dm+1); 
    ank(1,:) = [1 zeros(1,dm)]; 
        
    for j=1:n
        nck = zeros(dm,1); 
        for i=0:min(j,dm)
            nck(i+1) = nchoosek(j,i);
            ank(j+1,i+1)=nck(i+1)*prod(j+1:j+i)/(2^i);
        end   
    end
        
    SLRank = sum(LR.*ank,2); 
    MSJr_second = (exp(-lambda*r)/r)*SLRank.'; 

DLO = -real((-(lambda)*der_MSJ_first .*MSJr_second));
end
end





function [SL] = Spectra_nr2(n,lambda,r,out)
%{
This function computes the solid harmonic coefficients for the Modified
(Screened) Laplace equation Delta*u - lambda*u = 0, lambda>0.
 
Inputs
n      (int) spherical harmonics degree
lambda (double) positive real value for lambda
r      (double) target radius
out    (bool)   exterior (true) or interior (false)

Output
SL     (double) n+1 array of coefficients dependent on r and lambda. 
%}

if nargin<4
    out=1; 
end

p = (0:n)+0.5;
ilambda = (1i)*lambda;

if out
    % modified Bessel functions of first and second kind
    %J_n(lambda)
    J_pos = besselj(p,ilambda);
    MJ_first_pos = (1i.^(-p)).*J_pos;
    
    mxlr=4; 
    
    if lambda*r < mxlr
        %K_n(r*lambda)
        Jr_pos = besselj(p,ilambda*r);
        Jr_neg = besselj(-p,ilambda*r);
    
        MJr_first_pos = (1i.^(-p)).*Jr_pos;
        MJr_first_neg =(1i.^p).*Jr_neg; 
        MJr_second = (pi./(2*sin(p*pi))).*(MJr_first_neg - MJr_first_pos);

        % modified spherical Bessel functions of first and second kind
        MSJ_first = MJ_first_pos; %sqrt(pi/(2*lambda))*MJ_first_pos;
        MSJr_second = sqrt(1/r)*MJr_second; %sqrt(pi/(2*lambda*r))*MJr_second;

        %spectra
        SL = real(MSJ_first.*MSJr_second);
    else
        %approximation k_n(z) ~ (pi/2z)*exp(-z) as z->infty
        MSJ_first = sqrt(pi/(2*lambda))*MJ_first_pos;
        %MSJr_second = (exp(-lambda*r)/r)*ones(size(MSJ_first));
        lr = lambda*r; 
        
        dm = n; 
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
        
        SL = real(MSJr_second.*MSJ_first); 
    end
else

    % modified Bessel functions of first and second kind
    %J_n(r*lambda)
    Jr_pos = besselj(p,ilambda*r);
    
    MJr_first_pos = (1i.^(-p)).*Jr_pos;
    
    %K_n(lambda)
    J_pos = besselj(p,ilambda);
    J_neg = besselj(-p,ilambda);
    
    MJ_first_pos = (1i.^(-p)).*J_pos;
    MJ_first_neg = (1i.^p).*J_neg;
    MJ_second = (pi./(2*sin(p*pi))).*(MJ_first_neg - MJ_first_pos);

    % modified spherical Bessel functions of first and second kind
    MSJr_first = MJr_first_pos; %sqrt(pi/(2*lambda))*MJr_first_pos;
    MSJ_second = sqrt(1/r)*MJ_second; %sqrt(pi/(2*lambda*r))*MJ_second;

    %spectra
    SL = real(MSJr_first.*MSJ_second);
    
end


end


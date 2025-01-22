function [ SL ] = Spectra_n(n, lambda)
%n-fixed and lambda fixed value
%   Detailed explanation goes here

p = [0:n]+0.5;
ilambda = (1i)*lambda;

% modified Bessel functions of first and second kind
J_pos = besselj(p,ilambda);
J_neg = besselj(-p,ilambda);
MJ_first_pos = (1i.^(-p)).*J_pos;
MJ_first_neg = (1i.^p).*J_neg;
MJ_second = (pi./(2*sin(p*pi))).*(MJ_first_neg - MJ_first_pos);

% modified spherical Bessel functions of first and second kind
MSJ_first = sqrt(pi/(2*lambda))*MJ_first_pos;
MSJ_second = sqrt(pi/(2*lambda))*MJ_second;

%spectra
SL = real((2*lambda*MSJ_first.*MSJ_second)/pi);
end


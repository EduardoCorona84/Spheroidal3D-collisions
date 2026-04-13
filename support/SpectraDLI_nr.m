function [DL_I,der_MSJ_second] = SpectraDLI_nr(n,lambda,r )
% UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
p = [0:n]+0.5;
ilambda = (1i)*lambda;

% modified Bessel functions of first and second kind
Jr_pos = besselj(p,ilambda*r);
MJr_first_pos = (1i.^(-p)).*Jr_pos;

J_pos = besselj(p,ilambda);
J_neg = besselj(-p,ilambda);
MJ_first_pos = (1i.^(-p)).*J_pos;
MJ_first_neg = (1i.^p).*J_neg;
MJ_second = (pi./(2*sin(p*pi))).*(MJ_first_neg - MJ_first_pos);

% modified spherical Bessel functions of first and second kind
% MSJ_first = sqrt(pi/(2*lambda))*MJ_first _pos;
MSJr_first = sqrt(pi/(2*lambda*r))*MJr_first_pos;
MSJ_second = sqrt(pi/(2*lambda))*MJ_second;
% MSJr_second = sqrt(pi/(2*lambda*r))*MJr_second;
% MSJ_second(1:2)

% calculate the derivative of modified spherical Bessel function of first kind 
der_MSJ_second = der_MSJ_s(MSJ_second, lambda);

% spectra of double layer (inner) 
DL_I = -real((-2*lambda^2*der_MSJ_second.*MSJr_first)/pi);


end


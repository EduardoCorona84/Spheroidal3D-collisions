function [der_MSJ_second] = der_MSJ_s(MSJ_second,lambda)
% UNTITLED2 Summary of this function goes here
lambda=lambda';
q = size(MSJ_second,2);
s=size(lambda,1);
der_MSJ_second = zeros(s,q);
der_MSJ_second(:,1) = -((pi/2)*exp(-lambda).*(lambda+1))./lambda.^2;

for k = 2:q
der_MSJ_second(:,k) = (-1)^(k-2)*MSJ_second(:,k-1) - (k*(-1)^(k-1))*MSJ_second(:,k)./lambda;
der_MSJ_second(:,k) = (-1)^(k-1)*der_MSJ_second(:,k);
end

der_MSJ_second; 


end



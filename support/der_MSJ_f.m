
function [der_MSJ_first] = der_MSJ_f(MSJ_first,lambda)
% UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
lambda=lambda';
s=size(lambda,1);
q = size(MSJ_first,2);

der_MSJ_first = zeros(s,q);
%display(size((lambda.*cosh(lambda)-sinh(lambda))./(lambda.^2)));
der_MSJ_first(:,1) = (lambda.*cosh(lambda)-sinh(lambda))./(lambda.^2);


%MSJ_first=repmat(MSJ_first,s,1);
Lambda=repmat(lambda,1,q-1);

der_MSJ_first(:,2:q) = MSJ_first(:,1:q-1)-repmat((2:q),s,1).*MSJ_first(:,2:q)./Lambda;




end



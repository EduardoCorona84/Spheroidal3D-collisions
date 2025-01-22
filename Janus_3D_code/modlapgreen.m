function [G] = modlapgreen(x,y,lambda)
%y is center of the green's function (1 by 3)
%x is list of points to evaluate (#np by 3)
   y=repmat(y,size(x,1),1);
   del=x-y;
   r=sqrt(del(:,1).^2+del(:,2).^2+del(:,3).^2);
   G=exp(-lambda*r)./r;
    G=G/(8*pi);
end


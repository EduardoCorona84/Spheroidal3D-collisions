function S=cart2spheroidal(X,a,oblate)

if nargin<3
    oblate=false;
end

x2=size(X,2);
if fix(x2/3)~=x2/3
    error('Given points are not correct dimensions. #columns must be a multiple of 3.')
end

x3=size(X,3);
if length(a)==1
    a = a*ones(1,x3);
end

if length(a)~=x3
    error('length of a does not match with size of X')
end

S=zeros(size(X));

for i=1:length(a)
    ai=a(i);
    for j=1:x2/3
        x=X(:,3*j-2,i); 
        y=X(:,3*j-1,i); 
        z=X(:,3*j  ,i);
        if ~oblate
            u=1./(2*ai).*(sqrt(x.^2+y.^2+(z+ai).^2)+sqrt(x.^2+y.^2+(z-ai).^2));
            v=1./(2*ai).*(sqrt(x.^2+y.^2+(z+ai).^2)-sqrt(x.^2+y.^2+(z-ai).^2));
        else
            uhat=1./(2*ai).*(sqrt((sqrt(x.^2+y.^2)+ai).^2+z.^2)+sqrt((sqrt(x.^2+y.^2)-ai).^2+z.^2));
            vhat=1./(2*ai).*(sqrt((sqrt(x.^2+y.^2)+ai).^2+z.^2)-sqrt((sqrt(x.^2+y.^2)-ai).^2+z.^2));
            u=sqrt(uhat.^2-1);
            v=sqrt(-vhat.^2+1);
            v(z<0)=-1.*sqrt(-vhat(z<0).^2+1);
        end
        phi=atan2(y,x);
    
        S(:,3*j-2:3*j,i)=[u,v,phi];
    end
end
end

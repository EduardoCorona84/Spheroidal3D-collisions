function X=spheroidal2cart(S,a,oblate)

if nargin<3
    oblate=false;
end

s2=size(S,2);
if fix(s2/3)~=s2/3
    error('Given points are not correct dimensions. #columns must be a multiple of 3.')
end

s3=size(S,3);
if length(a)==1
    a = a*ones(1,s3);
end

if length(a)~=s3
    error('length of a does not match with size of X')
end

X=zeros(size(S));

for i=1:length(a)
    ai=a(i);
    for j=1:s2/3
        u=S(:,3*j-2,i); 
        v=S(:,3*j-1,i); 
        phi=S(:,3*j,i);
    
        if ~oblate
            x=ai.*sqrt(u.^2-1).*sqrt(1-v.^2).*cos(phi);
            y=ai.*sqrt(u.^2-1).*sqrt(1-v.^2).*sin(phi);
        else
            x=ai.*sqrt(u.^2+1).*sqrt(1-v.^2).*cos(phi);
            y=ai.*sqrt(u.^2+1).*sqrt(1-v.^2).*sin(phi);
        end
        z=ai.*u.*v;
    
        X(:,3*j-2:3*j,i)=[x,y,z];
    end
end



end
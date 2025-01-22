function [x1,x2,dst,nghlist]=MovingBallAlltoAll(par,tol,maxit,eta)

if nargin<4
    eta=1; 
end

nn = length(par); %number of spheroids

% Initialize arrays and cells 
C = zeros(nn,3); 
ax=zeros(nn,1); by=ax; cz=ax; rmx=ax; 
x1 = cell(nn,1); x2=x1; 

for i=1:nn
    C(i,:) = par(i).C; 
    ax(i) = par(i).a; by(i)=par(i).b; cz(i)=par(i).c;
    rmx(i) = max([ax(i) by(i) cz(i)]); 
end

% Initialize distance matrix using spheres with r=rmx. 
dst0 = dstSph(C,rmx); 
dst=dst0; 
idx = 1:nn; 
nghlist = cell(nn,1); 

% For neighbors, defined as pairs (i,j) for which dst(i,j)<= 2*eta*rmx(i), 
% explicitly calculate distances using moving ball algorithm

for i=1:nn
    neigh = (dst0(:,i) <= 2*eta*rmx(i) & dst0(:,i)~=0); 
    parneigh = par(neigh);

    nghlist{i} = idx(neigh); 

    [x1ngh,x2ngh,dstngh] = MovingBallNeighbors(par(i),parneigh,tol,maxit); 

    % update distances 
    dst(i,neigh) = dstngh; 
    dst(neigh,i) = dst(i,neigh)'; 

    x1{i} = x1ngh; 
    x2{i} = x2ngh; 
end

end

function dst = dstSph(C,R)

nC = size(C,1); 

[Y_g1,  X_g1  ] = meshgrid(C(:,1), C(:,1));
[Y_g2,  X_g2  ] = meshgrid(C(:,2), C(:,2));
[Y_g3,  X_g3  ] = meshgrid(C(:,3), C(:,3));
d1 = (X_g1 - Y_g1); d2 = (X_g2 - Y_g2); d3 = (X_g3 - Y_g3); 
den = sqrt(d1.^2 + d2.^2 + d3.^2);

Rmat = repmat(R(:),1,nC) + repmat(R(:)',nC,1) - diag(2*R(:)); 
dst = den-Rmat; 

end
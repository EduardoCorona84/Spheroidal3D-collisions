function [x1,x2,dst,nghlist,pars,counts,dst2,difnear,diffar] = test_MovingBallAlgo(nn,tol,p,comp)

if nargin<4
    comp=false; 
end

rng("default"); 
dbox = 8; 
CMat = dbox*2*(rand(nn,3)-0.5); %pick centers randomly from [-dbox,dbox]^3
Ccell = cell(nn,1); R=Ccell; a=R; b=R; c=R; shp=R; 
intr = zeros(nn,1); 
shapes = {'ellipseZ','oblateZ'};

for i=1:nn
    Ccell{i}=CMat(i,:);

    intr(i) = ceil(2*rand()); %random unif 1 or 2
    shp{i}=shapes{intr(i)}; 
    if intr(i)==1
        % prolate spheroid 
        a{i}=0.5; b{i}=0.5; c{i}=1; 
    elseif intr(i)==2
        % oblate spheroid 
        a{i}=1; b{i}=1; c{i}=2/sqrt(7);
    else
        fprintf("Error with intr");
        pause; 
    end

    [R{i},~]=qr(rand(3,3));
end

pars = struct('C',Ccell,'R',R,'a',a,'b',b,'c',c,'shape',shp);
maxit=100; neareta = 1; 
tic; 
[x1,x2,dst,nghlist]=MovingBallAlltoAll(pars,tol,maxit,neareta);
timealgo = toc; 
fprintf("\n Moving algo time for %d spheroids was %1.2e seconds",nn,timealgo); 

% Test distance against minimum distance between spheroidal meshes for mid
% to high p
if nargin<3
    p = 16;
end

np=2*p*(p+1); czero=0; cnear=0; cfar=0; 
Surf=cell(length(shapes)); Xs = Surf; 
dst2 = zeros(nn,nn); 

if comp
for k=1:length(shapes)
    Surf{k} = SurfaceSph(shape_gallery(p,shapes{k}));
    Xs{k} = reshape(Surf{k}.cart.to_array,[],3);
end
end

for i=1:nn-1
    for j=i+1:nn
        D = dst(i,j); 

        if D==0
            czero=czero+1; 
            fprintf("%d: dst(%d,%d) = 0 \n",czero,i,j);
        else
            if comp
            Xi = (R{i}*Xs{intr(i)}')' + repmat(CMat(i,:),np,1);
            Xj = (R{j}*Xs{intr(j)}')' + repmat(CMat(j,:),np,1);

            den = LOCAL_dstpts(Xi,Xj); den=den(:); 
            dst2(i,j) = min(den); 
            end

            if D <= 2*neareta
                cnear=cnear+1;
                difnear(cnear) = dst2(i,j)-dst(i,j);
                %fprintf("dst(%d,%d) <= 2*eta, dif = %1.2e \n",i,j,difnear(cnear));
            else
                cfar=cfar+1; 
                diffar(cfar) = dst2(i,j)-dst(i,j); 
                %fprintf("dst(%d,%d) > 2*eta, dif = %1.2e \n",i,j,diffar(cfar));
            end
        end

    end
end

for j=1:nn-1
    for i=j+1:nn
        dst2(i,j) = dst2(j,i); 
    end
end

counts = [czero cnear cfar];

end

function den = LOCAL_dstpts(X1,X2)

[Y_g1,  X_g1  ] = meshgrid(X2(:,1), X1(:,1));
[Y_g2,  X_g2  ] = meshgrid(X2(:,2), X1(:,2));
[Y_g3,  X_g3  ] = meshgrid(X2(:,3), X1(:,3));
d1 = (X_g1 - Y_g1); d2 = (X_g2 - Y_g2); d3 = (X_g3 - Y_g3); 
den = sqrt(d1.^2 + d2.^2 + d3.^2);

end
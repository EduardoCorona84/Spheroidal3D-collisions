function [x1,x2,dst,nghlist,pars,counts,dst2,difnear,diffar] = in_sph_MovingBallAlgo(nn,tol,p,comp,dbox)

if nargin<4
    comp=false; 
end

rng("default"); 
% dbox = 8; 
% % Default: square boxes
% CMat = dbox*2*(rand(nn,3)-0.5); %pick centers randomly from [-dbox,dbox]^3

% Special: Oblate spheroidal boundary with A=dbox and C=2/sqrt(7)*dbox
% CMat_box = dbox*2*(rand(5*nn,3)-0.5);
CMat_x=(dbox-1)*2.*(rand(5*nn,1)-0.5);
CMat_yz=3/5*(dbox-1)*2.*(rand(5*nn,2)-0.5);
CMat_box=cat(2,CMat_x,CMat_yz);
sph_domain=@(x,y,z) 1-((x+(-1).^(x<0)).^2./dbox^2+((z+(-1).^(z<0)).^2+(y+(-1).^(y<0)).^2)./(3*dbox/5)^2);
% sph_domain=@(x,y,z) 1-(((x+(-1).^(x<0)).^2+(y+(-1).^(y<0)).^2)./dbox^2+(z+(-1).^(z<0)).^2./(2/sqrt(7)*dbox)^2);

if_in_domain=sph_domain(CMat_box(:,1),CMat_box(:,2),CMat_box(:,3))>1e-3;
CMat=CMat_box(if_in_domain==1,:);
nn=size(CMat,1);

% scatter3(CMat(:,1),CMat(:,2),CMat(:,3),'r*'); hold on;
% Yout=prolate_spheroid_shape(8*2,5/4,dbox*4/5);
% Yout=([0 0 1;0 1 0;-1 0 0]*Yout')';
% scatter3(Yout(:,1),Yout(:,2),Yout(:,3),'bo'); hold off;
% axis equal


Ccell = cell(nn,1); R=Ccell; a=R; b=R; c=R; shp=R; 
intr = zeros(nn,1); 
shapes = {'ellipseZ','oblateZ','ellipseZ3','oblate85','oblate75'};

for i=1:nn
    Ccell{i}=CMat(i,:);

    % intr(i) = ceil(2*rand()); %random unif 0 or 1
    intr(i) = randi(5); % one random sample in {1,2,3,4,5}
    shp{i}=shapes{intr(i)}; 
    if intr(i)==1
        % prolate spheroid 'ellipseZ'
        a{i}=0.5; b{i}=0.5; c{i}=1; 
    elseif intr(i)==2
        % oblate spheroid  'oblateZ'
        a{i}=1; b{i}=1; c{i}=2/sqrt(7);
    elseif intr(i)==3
        % prolate spheroid 'ellipseZ3'
        a{i}=1/3; b{i}=1/3; c{i}=1;
    elseif intr(i)==4
        % oblate spheroid 'oblate85'
        a{i}=1; b{i}=1; c{i}=.47;
    elseif intr(i)==5
        % oblate spheroid 'oblate75'
        a{i}=1; b{i}=1; c{i}=.366;
    % elseif intr(i)==6
    %     % oblate spheroid 'oblate65'
    %     a{i}=1; b{i}=1; c{i}=.29;
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

difnear=[];diffar=[];

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
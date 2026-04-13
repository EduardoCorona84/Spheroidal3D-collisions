function [colevent,collist,mindst,mindstsh] = LOCAL_check_collision_sph(C,Fparams)

n3 = size(C,1); 

rd = Fparams.parbd.rd; 
mxrd = Fparams.parbd.mxrd; 
diam = Fparams.parbd.diam; 
eps = Fparams.parbd.eps;  

if n3>1
    % Compute center distances
    distC = LOCAL_CenterDistance(C);

    % Find pairs for which (C_i-C-j) <= (r_i+r_j)+1.1*eps*max(r_i,r_j)
    [ii,jj]=meshgrid(1:n3); 
    id = distC<=diam+1.1*eps*mxrd & ii<jj; 
    ip = ii(id); 
    jp = jj(id); 

    % Compute minimum relative distance between spheres
    mindst=min(reshape((distC-diam)./mxrd,[],1));
else
    ip=[]; jp=[]; 
    mindst=Inf; 
end

% If there is a spherical shell, compute signed distance to boundary
if isfield(Fparams,'parsh')
   rdsh = Fparams.parsh.rd; 
   sheps = Fparams.parsh.eps;
     
   NC = sqrt(sum(C.*C,2)); 
   distSh = rdsh - rd - NC; %(R-r) - ||C_i||
   mindstsh = min(distSh)/Fparams.parsh.rd; 
   iish = find(distSh <= 1.1*sheps*rdsh); 
   jjsh = (n3+1)*ones(size(iish)); %j=n3+1 -> collision with boundary
   
   colevent = mindst < 1.1*eps || mindstsh < 1.1*sheps;
   collist = [[ip ; iish] [jp ; jjsh]];
else
   colevent = mindst < 1.1*eps; 
   mindstsh = Inf;
   collist = [ip jp]; 
end

end

function den = LOCAL_CenterDistance(C)

[Y_g1,  X_g1  ] = meshgrid(C(:,1), C(:,1));
[Y_g2,  X_g2  ] = meshgrid(C(:,2), C(:,2));
[Y_g3,  X_g3  ] = meshgrid(C(:,3), C(:,3));
d1 = (X_g1 - Y_g1); d2 = (X_g2 - Y_g2); d3 = (X_g3 - Y_g3); 
den = sqrt(d1.^2 + d2.^2 + d3.^2);
den = 10000*(den==0)+den; 

end
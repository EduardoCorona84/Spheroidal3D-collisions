function Vn = Vnm(type,Sc,n,m,u,v,Nr)

ln = length(n); 

if ln==1
   Y = Ynm(n,m,u,v);
else
   % May edit in the future to allow more flexibility 
   n = max(n); 
   m = []; 
   sp=(n+1)^2; 
   Y = zeros(size(u,1),sp);
   for p=0:n
      idx = p^2+1:(p+1)^2;
      Y(:,idx) = Ynm(p,[],u,v); 
   end
end

[np,nV] = size(Y);

if nargin==6
    Nr = reshape(Sc.geoProp.nor.to_array,[],3);
    Gfun = @(mm,f) LOCAL_geoProp_GradY(Sc,f); 
else
    if isempty(Nr)
        Nr = reshape(Sc.geoProp.nor.to_array,[],3);
    end
    
    if ln>1
       Yp = [Y Ynm(n+1,[],u,v)];
       Gfun = @(mm,f) LOCAL_GradY(n,mm,u,v,Yp,ln); 
    else
       Gfun = @(mm,f) LOCAL_GradY(n,mm,u,v,Y,ln); 
    end
end

if isempty(m)
    if ~strcmp(type,'YNr')
        GYL = Gfun([],[]);
    end
else
    GYL = Gfun(m,Y); 
end

if ln>1
   nlist = floor(sqrt(0:((n+1)^2-1)));
elseif isempty(m)
   nlist = repmat(n,1,2*n+1);  
else
   nlist=n; 
end

szn = size(nlist,2); 

if ~strcmp(type,'YNr')
    GY = reshape(permute(reshape(GYL,[np 3 szn]),[1 3 2]),[],3); 
end

if ~strcmp(type,'GY')
    Nrp = repmat(Nr,szn,1); 
end

switch type
   case 'Vnm'
       Y = -repmat(nlist+1,np,1).*Y; 
       Vtmp = GY + repmat(Y(:),1,3).*Nrp; 
   case 'Wnm'
       Y = repmat(nlist,np,1).*Y;  
       Vtmp = GY + repmat(Y(:),1,3).*Nrp; 
   case 'Xnm' 
       Vtmp = cross(Nrp,GY); 
   case 'VWX'
       YV = -repmat(nlist+1,np,1).*Y;
       YW = repmat(nlist,np,1).*Y; 
       
       Vtmp = {GY + repmat(YV(:),1,3).*Nrp,GY + repmat(YW(:),1,3).*Nrp,cross(Nrp,GY)};
   case 'GY'
       Vtmp = GY;  
   case 'YNr'
       Vtmp = repmat(Y(:),1,3).*Nrp;
   otherwise
       Vtmp = repmat(Y(:),1,3);  
end

if ~strcmp(type,'VWX') && ~strcmp(type,'GY')
    Vn = reshape(permute(reshape(Vtmp,[np szn 3]),[3 1 2]),[],szn);
elseif strcmp(type,'VWX')
    Vn = {reshape(permute(reshape(Vtmp{1},[np szn 3]),[3 1 2]),[],szn),...
          reshape(permute(reshape(Vtmp{2},[np szn 3]),[3 1 2]),[],szn),...
          reshape(permute(reshape(Vtmp{3},[np szn 3]),[3 1 2]),[],szn)};
else
    Vn = reshape(permute(reshape(Vtmp,[np szn 3]),[1 3 2]),[],szn);
end

end

function GY=LOCAL_GradY(n,m,u,v,Y,ln)
nu = size(u,1); 

if ln>1
    ev = [-sin(v)./sin(u) cos(v)./sin(u) zeros(nu,1)]; 
    eu = [cot(u).*cos(v) cot(u).*sin(v) -1*ones(nu,1)];
    
    GY = zeros(3*nu,(n+1)^2); 
    %cpnm = @(n,mm) sqrt((n.^2).*((n+1).^2-mm.^2)./((2*n+1).*(2*n+3)));
    %cmnm = @(n,mm) -sqrt(((n+1).^2).*(n.^2-mm.^2)./((2*n+1).*(2*n-1))); 
    %idY = @(p) (p^2+1):(p+1)^2; 
    idY = cell(n+1,1); 
    idY{1} = 1; 
    
    for p=0:n
        mp = -p:p; 
        mp2 = (-p:p).^2; 
        mmi = repmat(1i*mp,nu,1); 
        idY{p+2} = ((p+1)^2+1):(p+2)^2; 
        %idYp = (p^2+1):(p+1)^2;
        
        Yv = mmi.*Y(:,idY{p+1});
        Yp = Y(:,idY{p+2}(2:end-1)); %Yp = Yp(:,2:end-1);
        
        cp = (p/sqrt((2*p+1)*(2*p+3)))*sqrt((p+1)^2-mp2); %cpnm(p,mm); 
        Yu = repmat(cp,nu,1).*Yp; 
        if p>0
            cm = (-(p+1)./sqrt((2*p+1)*(2*p-1)))*sqrt(p^2-mp2);
            cm = repmat(cm,nu,1);  
            Yu(:,2:end-1) = Yu(:,2:end-1)+cm(:,2:end-1).*Y(:,idY{p}); 
        end
        
        %su = repmat(sin(u),1,2*p+1); 
        %GY(:,idY(p)) = repmat(Yu,3,1).*repmat(eu(:),1,2*p+1) + repmat(Yv,3,1).*repmat(ev(:),1,2*p+1);
        idd = idY{p+1}; 
        for k=1:2*p+1
           GY(:,idd(k)) = [Yu(:,k);Yu(:,k);Yu(:,k)].*eu(:) + [Yv(:,k);Yv(:,k);Yv(:,k)].*ev(:); 
        end
    end
elseif isempty(m)
    ev = [-sin(v) cos(v) zeros(nu,1)]; 
    eu = [cos(u).*cos(v) cos(u).*sin(v) -sin(u)];
    
    mm = repmat(-n:n,nu,1); 
    %Y = Ynm(n,[],u,v); 
    Yv = 1i*mm.*Y; 
    
    cpnm = sqrt((n.^2).*((n+1).^2-mm.^2)./((2*n+1).*(2*n+3)));
    cmnm = -sqrt(((n+1).^2).*(n.^2-mm.^2)./((2*n+1).*(2*n-1))); 
    
    Yp = Ynm(n+1,[],u,v); Yp = Yp(:,2:end-1); 
    
    Yu = cpnm.*Yp; 
    if n>0
    Yu(:,2:end-1) = Yu(:,2:end-1)+cmnm(:,2:end-1).*Ynm(n-1,[],u,v); 
    end
    
    Yu = Yu./repmat(sin(u),1,2*n+1); Yv=Yv./repmat(sin(u),1,2*n+1); 
     
    GY = repmat(Yu,3,1).*repmat(eu(:),1,2*n+1) + repmat(Yv,3,1).*repmat(ev(:),1,2*n+1);
    
else
    ev = [-sin(v) cos(v) zeros(nu,1)]; 
    eu = [cos(u).*cos(v) cos(u).*sin(v) -sin(u)];
    
    Yv = 1i*m*Y; 
     
    cpnm = sqrt((n^2)*((n+1)^2-m^2)/((2*n+1)*(2*n+3)));
    cmnm = -sqrt(((n+1)^2)*(n^2-m^2)/((2*n+1)*(2*n-1)));
    
    Yu = cpnm*Ynm(n+1,m,u,v); 
    
    if abs(m)<n && n>0
       Yu = Yu+cmnm*Ynm(n-1,m,u,v); 
    end
    
    Yu = Yu./sin(u); Yv=Yv./sin(u); 
    
    GY = repmat(Yu,1,3).*eu + repmat(Yv,1,3).*ev; 
end

end

function GY = LOCAL_geoProp_GradY(Sc,f)
   GY = Sc.geoProp.Grad(f);
   GY = reshape(GY.to_array,[],3);
end


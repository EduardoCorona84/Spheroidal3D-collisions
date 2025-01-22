function U = Vsh_Kernel_Eval(F,type)

SVeg  = @(n) n./((2*n+1).*(2*n+3)); 
SWeg  = @(n) (n+1)./((2*n+1).*(2*n-1));
SpVeg = @(n) (3/2)./((2*n+1).*(2*n+3));
SpWeg = @(n) -(3/2)./((2*n+1).*(2*n-1));
SXeg  = @(n) 1./(2*n+1); 
SpXeg = @(n) -(3/2)./(2*n+1); 

[d1,d2] = size(F);
np = d1/3;
p = (sqrt(2*np+1)-1)/2;
sp=(p+1)^2;
nn = floor(sqrt((1:sp)'-1));

if nargin==1
    type='SMat'; 
end 

% phi = Sum(shV_nm*Vnm(u,v) + shW_nm*Wnm(u,v) + shX_nm*Xnm(u,v))
shc = VshAna(phi,'VW'); 

if strcmp(type,'SMat')
    Lam = repmat([SVeg(nn);SWeg(nn);SXeg(nn)],1,d2); 
elseif strcmp(type,'SpMat')
    Lam = repmat([SpVeg(nn);SpWeg(nn);SpXeg(nn)],1,d2); 
end

% u = Sum(lV_nm*shV_nm*Vnm(u,v)+lW_nm*shW_nm*Wnm(u,v)+lX_nm*shX_nm*Xnm(u,v))
U = VshSyn(Lam.*shc,'VW'); 

end


function u = Sh_Kernel_Eval(phi,type)

LSeg  = @(n) (1./(2*n+1)); 
LSpeg = @(n) -(1./(4*n+2));  

[np,d2] = size(phi); 
p = (sqrt(2*np+1)-1)/2;
sp= (p+1)^2; 
nn = floor(sqrt((1:sp)'-1)); 

if nargin==1
    type='SMat'; 
end 

% phi = Sum(sh_nm*Ynm(u,v))
shc = shAna(phi); 

if strcmp(type,'SMat')
    Lam = repmat(LSeg(nn),1,d2); 
elseif strcmp(type,'SpMat')
    Lam = repmat(LSpeg(nn),1,d2);
end

% u = Sum( l_nm*f_nm*Ynm(u,v))
u = shSyn(Lam.*shc); 

end
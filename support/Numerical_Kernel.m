function [ K ] = Numerical_Kernel( Xtrg,Sc,p,lambda,type )
%creates matrix kernel ,K, for smooth numerical quadrature
Xsource = reshape(Sc.cart.to_array,[],3); 
[u,v]=gl_grid(p); 



sp = @(p) (p+1)^2;  %(number of harmonic coefficients);
np = @(p) 2*(p+1)*p; %(number of quadrature points);
ii=(1:sp(p))';
nn=floor(sqrt(ii-1)); 

if strcmp(type,'SMat')
    partype='SL_LMOD_3D';
elseif strcmp(type,'DMat')
    partype='DL_LMOD_3D';
end
 par = Kernel_Eval_parameters(partype,0,1,1,1,1e-8,2,400,1);
    par.dim=3;
    par.lambda=lambda;
    [~, gwt]=g_grid(p+1);
    wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
    wt = wt(:);
    W = Sc.geoProp.W; W= W.*wt; 
    Nr = reshape(Sc.geoProp.nor.to_array,[],3); 
    par.X = Xsource; par.nor = Nr; par.W2 = W.'; 
    Rtrg = sqrt(sum(Xtrg'.*Xtrg'))';
    K=Kernel_Eval(Xtrg,Xsource,par);
end


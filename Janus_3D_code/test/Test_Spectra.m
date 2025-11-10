function [ Minaccuracy ] = Test_Spectra( )
%Unit tests for Get_Spectra function
%generates a smooth test potential and compares eigenvalues in the 
%non-singular case for interior and exterior problems
maxp=24;
p=16;
lambda=1e-8;
r=4;
rng('default');
sp = @(p) (p+1)^2;
swgt=rand(sp(maxp),1);

for i=1:sp(maxp)        %exponentially decaying weights
    swgt(i)=swgt(i)*(1/2)^(floor(sqrt(i-1)));
end


swgt(2)=1;
display(swgt);
%Generates Spectra to be used
ExpSmoothSpec=swgt(1:sp(p));
ExpSmoothSpec(sp(p-1)+1:sp(p))=0;  %filters out inaccurate final freq
SmoothFunc=shSyn(ExpSmoothSpec);   %creates functions from Spectra
display(size(ExpSmoothSpec));

errorr(1)=test_kernel(ExpSmoothSpec,SmoothFunc,'SMat',lambda,p,r,1);
errorr(2)=test_kernel(ExpSmoothSpec,SmoothFunc,'SMat',lambda,p,1/r,0);
errorr(3)=test_kernel(ExpSmoothSpec,SmoothFunc,'DMat',lambda,p,r,1);
errorr(4)=test_kernel(ExpSmoothSpec,SmoothFunc,'DMat',lambda,p,1/r,0);
errorr
Minaccuracy=ceil(max(errorr(1),errorr(3)));
if Minaccuracy >-10
    error('less than 10 places accuracy');
end


end

function [error]=test_kernel(Spectrum,Function,type,lambda,p,r,out)

Sc = SurfaceSph(shape_gallery(p,''));
X = reshape(Sc.cart.to_array,[],3); 
[u,v]=gl_grid(p); 



sp = @(p) (p+1)^2;  %(number of harmonic coefficients);
np = @(p) 2*(p+1)*p; %(number of quadrature points);
ii=(1:sp(p))';
nn=floor(sqrt(ii-1)); 
%{
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
    par.X = X; par.nor = Nr; par.W2 = W.'; 
    Xtrg = r*X;
    Rtrg = sqrt(sum(Xtrg'.*Xtrg'))';
    K=Kernel_Eval(Xtrg,X,par);
  %}
    K=Numerical_Kernel(r*X,Sc,p,lambda,type);
    Numerical=K*Function;
    
    TrueEvals=Get_Spectra(p,lambda,r,type,out);
    for i=1:sp(p)
        n=nn(i);
        TrueSpec(i)=Spectrum(i)*TrueEvals(n+1);
    end
    Theoretical=shSyn(TrueSpec');
   error=norm(Numerical-Theoretical)/norm(Theoretical);
   error=log10(error);



end

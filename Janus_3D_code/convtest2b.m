function [ RelVals, AbsVals,Spectra] = convtest2( p,lambda,rt,type,singular,rs)
%Tests numerically computed solutions against those computed from
%eigenvalue formula computes two cases: exponentially decaying SH
%coefficients and uniformly random coefficients. If non singular,
%the target points are on a sphere of radius rt != 1
%Inputs: 
    % p: degree of Spherical Harmonics note: coefficients are computed at 
    %       start of function, in order to allow for varying p values,
    %       a max p is set at 24. This may be modified as needed.
    %lambda: the value in Lap(f)-(lambda)f=0
    %rt:target radius (if singular==1,rt is not used
    %type: 'SMat' or 'DMat': Single or Double layer
    %singular: if true, uses singular quadrature, if false or not
    %   specified, uses smooth quadrature
    %rs: radius of source sphere-defaults to 1
%Outputs:
    %RelVals: relative error between numerically computed and formula 
        %computed vectors
    %AbsVals: absolute error
    
if nargin<4
    error('too few arguments');

elseif nargin<5
    singular=0;
    rs=1;
elseif nargin<6
    rs=1;
end

if type=='SMat'
    partype='SL_LMOD_3D';
    scaling=rs;
elseif type=='DMat'
    partype='DL_LMOD_3D';
    scaling=1;
end
maxp=24;
sp = @(p) (p+1)^2;  %(number of harmonic coefficients);
np = @(p) 2*(p+1)*p; %(number of quadrature points);
ii=(1:sp(p))';
nn=floor(sqrt(ii-1)); %indices of n
%generate random seeds
rng('default');
if (p>maxp)
    error('p exceeds maximum allowed value');
end
nswgt=rand(sp(maxp),1);  %generate uniform random weights
swgt=rand(sp(maxp),1);
%{
for i=1:sp(maxp)        %exponentially decaying weights
    swgt(i)=swgt(i)*(1/2)^(floor(sqrt(i-1)));
end
%}
swgt=zeros(sp(p),1);
swgt(4)=1;
%Generates Spectra to be used
ExpSmoothSpec=swgt(1:sp(p));
ExpNonSmoothSpec=nswgt(1:sp(p));
ExpSmoothSpec(sp(p-1)+1:sp(p))=0;  %filters out inaccurate final freq
ExpNonSmoothSpec(sp(p-1)+1:sp(p))=0;

SmoothFunc=shSyn(ExpSmoothSpec);   %creates functions from Spectra
NonSmoothFunc=shSyn(ExpNonSmoothSpec);



Sc = SurfaceSph(shape_gallery(p,''));
X = rs*reshape(Sc.cart.to_array,[],3); 
[u,v]=gl_grid(p); 

%creates kernel
if singular==1
     
    [~,~,K]=kernelModifiedLap(Sc,type,lambda);
      if type == 'SMat'
          Spectra=Spectra_n(p,lambda);
      elseif type == 'DMat'
          Spectra=(SpectraDLO_nr2(p,lambda,1)+SpectraDLI_nr(p,lambda,1))/2;
          
      end
else
    
    par = Kernel_Eval_parameters(partype,0,1,1,1,1e-8,2,400,1);
    par.dim=3;
    par.lambda=lambda;
    [~, gwt]=g_grid(p+1);
    wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
    wt = wt(:);
    W = Sc.geoProp.W; W= W.*wt; 
    Nr = reshape(Sc.geoProp.nor.to_array,[],3); 
    par.X = X; par.nor = Nr; par.W2 = rs^2*W.'; 
    Xtrg = rt*X;
    Rtrg = sqrt(sum(Xtrg'.*Xtrg'))';
    K=Kernel_Eval(Xtrg,X,par);
    if type=='SMat'
        Spectra=Spectra_nr2(p,rs*lambda,rt,(rt>1));
    elseif type=='DMat'
        if rt>1
            Spectra=SpectraDLO_nr2(p,lambda*rs,rt);
        elseif rt<1
            Spectra=SpectraDLI_nr(p,lambda*rs,rt);
        
        
        end
    end
end

NumericalSmooth=K*SmoothFunc;
NumericalNonSmooth=K*NonSmoothFunc;

%D=(shAna(NumericalSmooth));
%display(D([1 4 9 16 25]));



for i=1:sp(p)
    indexn=nn(i);
    Eigenswgt(i)=ExpSmoothSpec(i)*Spectra(indexn+1);
    Eigennswgt(i)=ExpNonSmoothSpec(i)*Spectra(indexn+1);
end
%display(Eigenswgt([1 4 9 16 25]));
Eigenswgt=Eigenswgt';
Eigennswgt=Eigennswgt';
TheoreticalSmooth=shSyn(Eigenswgt)*scaling;
TheoreticalNonSmooth=shSyn(Eigennswgt)*scaling;

RelVals(1)=norm(real(TheoreticalSmooth-NumericalSmooth))/norm(TheoreticalSmooth);
RelVals(2)=norm(real(TheoreticalNonSmooth-NumericalNonSmooth))/norm(TheoreticalNonSmooth);
RelVals=log10(RelVals);

AbsVals(1)=norm(real(TheoreticalSmooth-NumericalSmooth));
AbsVals(2)=norm(real(TheoreticalNonSmooth-NumericalNonSmooth));
AbsVals=log10(AbsVals);


end


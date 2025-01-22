function [SYnssh,SYns]=SDY(p,Mat,oblate,fvals)
% Uses singular quadrature to compute double layer of Ynm 
% or single layer of Ynm/(uu^2-cos^2(u))^(1/2) on
% surface of spheroid with eccentricity 1/uu, semimajor axis 1/uu
% Parameters: p=order of spheroidal harmonic expansion, 
% Mat=layer operator ('DL', 'SL', or 'SP') 
if nargin==2
    oblate=false;
    [u,v]=gl_grid(p);
end

%ellipse with eccentricity 1/u, u=2/sqrt(3)
uu = 2/sqrt(3);
% Shape='elipseZ'; 


% Compute Laplace kernel using singular quadrature
%----------------------------------------------------------%

% Spheroid surface
if oblate
    Shape='oblateZ';
else
    Shape='ellipseZ'; 
end
Sns = SurfaceSph(shape_gallery(p,Shape)); 
Xns = reshape(Sns.cart.to_array,[],3);
% Laplace SL, S' and DL on the spheroid surface
[SMns,Spns,DMns] = kernelDLap(Sns);

%----------------------------------------------------------%


% Compute spheroidal harmonic functions and eigenvalues from formula
%----------------------------------------------------------%
np = 2*p*(p+1); sp=(p+1)^2;
ii = (1:sp)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;
[u,v]=gl_grid(p);

if nargin<4
    Y = zeros(np,sp); SYns = Y;  Y2=Y; 
    
    for k=1:sp
        n = nn(k); m = mm(k); 
        Y(:,k) = Ynm(n,m,u,v); 
        if ~oblate
            fac=(1./sqrt(uu^2-cos(u).^2)); fac=fac(:); 
        else
            fac=(1./sqrt(uu^2+cos(u).^2)); fac=fac(:);
        end
        Y2(:,k) = fac.*Y(:,k);
    end
else
    Y = fvals; Y2 = fvals;
end

% Apply Laplace kernel to spheroidal harmonic functions, and 
% extract eigenvalues from singular quadrature
%----------------------------------------------------------%
if strcmp(Mat,'SL')
    SYns = SMns*Y2;
    SYnssh = shAna(SYns);
elseif strcmp(Mat,'DL')
    SYns = DMns*Y;
    SYnssh = shAna(SYns);
elseif strcmp(Mat,'SP')
    SYns= Spns*Y2;
    SYnssh = shAna(SYns);
else
    sprintf("Invalid input for Layer Operator. Must be 'SL' for single layer, 'DL' for double layer, or 'SP' for S'.")
    return

%----------------------------------------------------------%

end
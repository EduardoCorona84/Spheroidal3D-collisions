function DL=spheroidalDL_mat(p,u0,X,target_coords)

if nargin<4
    target_coords='cart';
end
if nargin<3
    fprintf('off-surface targets not given. Evaluating nodes on-surface.\n')
    [theta_x,phi_x]=gl_grid(p);
    nt=length(theta_x);
    u_x=u0.*ones(nt,1); %if u_x not given, calculate on-surface
    X=cat(2,u_x,cos(theta_x),phi_x);
    target_coords='spheroidal';
end

[nt,d4]=size(X);
if d4~=3
    error("Dimensions of targets should be M x 3.")
end

% sort targets by increasing u so we can split up interior/surface/exterior
% Convert targets to spheroidal coords
if strcmp(target_coords, 'cart')
    S=cart2spheroidal(X,1/u0);
elseif strcmp(target_coords, 'spheroidal')
    S=X;
else 
    error("Invalid coordinates given. Use 'cart' for cartesian or 'spheroidal' for spheroidal.")
end

u_x=S(:,1);
v_x=S(:,2);
phi_x=S(:,3);

% Calculate spectra for interior, exterior, and surface
[spectra_int,spectra_surf,spectra_ext]=DLspectrum(p,u0);
if all(u_x>u0) %exterior
    region=1; 
    spectra=spectra_ext;
elseif all(u_x==u0) %coincident
    region=0; 
    spectra=spectra_surf;
elseif all(u_x<u0) %interior
    region=-1; 
    spectra=spectra_int;
else
    error("Given targets must all be in the same region: exterior,interior, or coincident")
end

F=solid_harmonic(p,u_x,region);
Y=zeros(nt,(p+1)^2);

for n=0:p  %loop over terms in spheroidal harmonic expansion
    Yn=Ynm(n,[],acos(v_x)',phi_x);
    Y(:,n^2+1:(n+1)^2)=Yn;
end

DL=F.*Y.*repmat(spectra',nt,1);


% For debugging:
%---------------------------------------------------------------
%         sprintf('size Y: ')
%         disp(size(Y))
%         sprintf('size spectra: ')
%         disp(size(spectra))
%         sprintf('size shc: ')
%         disp(size(shc))
%         sprintf('size FY: ')
%         disp(size(FY))
%---------------------------------------------------------------

end


function [lambda_int,lambda_surf,lambda_ext]=DLspectrum(p,u0)
    sp=(p+1)^2;
    ii = (1:sp)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;
    anm=factorial(nn-mm)./factorial(nn+mm).*(-1).^mm .*(u0.^2-1);
    L=legendre_otc(p,u0,1,1,1);
    P=L{1}; Q=L{2}; dP=L{3}; dQ=L{4};

    lambda_int=anm.*dQ;
    lambda_surf=anm./2.*(P.*dQ+dP.*Q);
    lambda_ext=anm.*dP;   
end

function F=solid_harmonic(p,u_x,region)
    if region==-1
        PQ=legendre_otc(p,u_x,0);
        P=PQ{1};
        F=P';
    elseif region==0
        F=1;
    elseif region==1
        PQ=legendre_otc(p,u_x,1);
        Q=PQ{2};
        F=Q';
    else
        error("invalid region given.")
    end
end
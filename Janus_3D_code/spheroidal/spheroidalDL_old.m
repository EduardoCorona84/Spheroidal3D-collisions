function DL=spheroidalDL_old(f,u,u_x,isReal)
%--------------------------------------------------------------------%
% spheroidalDL computes the laplace Double Layer potential of a density 'f'
% on the surface of a spheroid.
% Parameters:
%    (o) f = function on the spheroid surface, as a function of (theta,phi). 
%        See gl_grid. cos(theta) are gauss-Legendre nodes and v are 
%        equispaced.
%        f is size [d0,d2] and each column of f corresponds to a different 
%        spheroid surface.
%    (o) u = 1/eccentricity of the spheroid surface.
%    (o) u_x = u coordinates (prolate spheroidal coordinate system) of 
%        off-surface evaluation. Is size [d3,d2] where each column is a set
%        of u points centered at the spheroid corresponding to that column
%        of f.
%    (o) isReal = Can be set if real output is expected. then imaginary 
%        parts are removed 
%    
% Returns DL of size [d0,d3,d2]:
%    (o) DL[:,j,k] = D[f] for a single spheroid and single u-coordinate, as
%        a function of (theta,phi).
%    (o) DL[i,:,k] = D[f] for a single spheroid at fixed (theta,phi), as a
%        function of u-coordinate of target points
%    (o) DL[i,j,:] = D[f] at fixed (theta,phi) and single u-coordinate, on 
%        different spheroid surfaces.
% 
%--------------------------------------------------------------------%

[d0,~]=size(f);
shc=spheroidalAna(f);
[d1,d2]=size(shc);

if nargin<4
    isReal=false;
end
if nargin<3
    fprintf('off-surface targets not given. Evaluating only on-surface...\n')
    u_x=u.*ones([1,d2]); %if u_x not given, calculate on-surface
end

[d3,d4]=size(u_x);
if d4~=d2
    error("Dimensions of surface density and target point array do not match.")
end

DL=zeros([d0,d3,d2]);

for k=1:d2  %loop over each spheroid surface we want to evaluate
    for j=1:d3  %for a single surface, loop over target points
        
        spectra=zeros([d1,1]);
        
        for i=1:d1  %loop over terms in spheroidal harmonic expansion
            n=floor(sqrt(i-1));
            m=i-1-n*(n+1);
            spectra(i)=spectrum(n,m,u,u_x(j,k));
        end

        DLshc=spectra.*shc(:,k);
        DL_jk=spheroidalSyn(DLshc,isReal);
        DL(:,j,k)=DL_jk;

    end

end


DL_ijk=spheroidalSyn(DLshc,isReal);

end

function val=spectrum(n,m,u,u_x)
    [P,Q]=ALF(n,u_x);
    Pnm=(factorial(n+m)./factorial(n-m)*(m<0)+(m>=0)).*P(abs(m)+1,:);
    Qnm=(factorial(n+m)./factorial(n-m)*(m<0)+(m>=0)).*Q(abs(m)+1,:);
    [dP,dQ]=dALF(n,m,u);

    anm=factorial(n-m)./factorial(n+m).*(-1).^m .*(u.^2-1);
   
    if u_x>u %exterior
        val=anm.*(Qnm.*dP);
    elseif u_x==u %surface
        val=anm.*(Pnm.*dQ+Qnm.*dP)./2;
    else %interior
        val=anm.*(Pnm.*dQ);
    end

%     [dP,dQ]=dALF(n,m,u,1)    %mode 1: scale dP and dQ
%     val=(-1).^m .*(u.^2-1).*(Pnm.*dQ+Qnm.*dP)./2;
end
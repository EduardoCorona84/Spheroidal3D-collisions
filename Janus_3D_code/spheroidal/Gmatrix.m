function G=Gmatrix(p,u,fig,oblate)
% Generates the G matrix which transforms from the basis 
% Y_n^m(v,phi) /sqrt(u^2-v^2) to Y_n^m(v,phi).
% p=order of spheroidal harmonic expansion
% u= 1/eccentricity of ellipse
% fig= generates a plot of G if 1
if nargin<4
    oblate=false;
end
if nargin<3
    fig=0;
end

nsh=(p+1).^2; %number of spherical harmonics in expansion for p
ii = (1:nsh)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;

G=zeros(nsh);
Gspy=zeros(nsh);

for i=1:nsh
    ni= nn(i); mi = mm(i); 
    for j=1:nsh
        nj= nn(j); mj = mm(j); 
        G(i,j)=project(ni,mi,nj,mj,50,u,oblate);
        if G(i,j)>1e-10
            Gspy(i,j)=G(i,j);
        end
    end
end

if fig
    % figure;
    % imagesc(G)
    % colorbar
    % title("before sort")
    % 
    % [~,j1]=sort(mm);
    % figure;
    % imagesc(G(j1,j1))
    % colorbar
    % title("after sort")
    % 
    % [nn2,j2]=sort(nn(j1));
    % G1=G(j1,j1);
    % figure;
    % imagesc(G1(j2,j2))
    % colorbar
    % title("sorted back to original")
    % 
    % disp(sum(abs(G-G1(j2,j2)),'all'))

    figure;
    spy(Gspy)
end

end


function proj=project(ni,mi,nj,mj,vnp,u,oblate)
% Compute projection of Y_{ni}^{mi}(v,phi)/sqrt(u^2-v^2) onto new basis 
% element Y_{nj}^{mj}(v,phi):
% 
%    <Y_{ni}^{mi}(v,phi),Y_{nj}^{mj}(v,phi)/sqrt(u0^2-v^2)> 
%   
% on the spheroid. 
% 
% *Notation: u0= 1/eccentricity of generating ellipse, 
% Ynm(v,phi)=Cnm*P_n^m(v)*e^{i*m*phi}, -1<=v<=1, and 0<=phi<=2pi, 
% Cnm normalization constant.
% vnp = number of v points for GL quadrature

if mi ~= mj
    proj=0;
    return
else
    m=mi;
    mci=((-1).^m.*factorial(ni+m)./factorial(ni-m)*(m<0)+(m>=0)); %coefficient for negative m
    mcj=((-1).^m.*factorial(nj+m)./factorial(nj-m)*(m<0)+(m>=0)); %coefficient for negative m

    [vp,vw]=g_grid(vnp); %Gauss-Legendre Nodes and weights
    
    Pni=legendre2(ni,vp);
    Pnmi=mci.*Pni(abs(m)+1,:);
    Pnj=legendre2(nj,vp);
    Pnmj=mcj.*Pnj(abs(m)+1,:);
    
    %inner product
    if oblate
        integrand1=Pnmi.*Pnmj ./ ((u.^2+vp.^2).^(1/2));
    else
        integrand1=Pnmi.*Pnmj ./ ((u.^2-vp.^2).^(1/2));
    end

    % %%% Check inner product
    % integrand1 = Pnmi.*Pnmj;
    % %%%%%%%%%%%%%%%%%%%%%

    integral1=sum(vw.*integrand1);

    nci=sqrt((ni+1/2) .* factorial(ni-m)./factorial(ni+m) );
    ncj=sqrt((nj+1/2) .* factorial(nj-m)./factorial(nj+m) );

    proj=ncj.*nci.*integral1 ;
end
end




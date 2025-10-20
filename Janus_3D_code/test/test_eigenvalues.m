% close all;

uu = 2/sqrt(3);
% uu = 10;
% uu = 10;
Shape='elipseZ';
rt=4;
Mat='SL';
verb=1;

parr=1:16;
m0err=zeros(1,length(parr));
fixn=4;

for l=1:length(parr)
p=parr(l);

% Sphere
Sph = SurfaceSph(shape_gallery(p,''));
Xsp = reshape(Sph.cart.to_array,[],3); 
% Laplace SL, S' and DL
[SMsp,Spsp,DMsp] = kernelDLap(Sph); 

% Test surf (non-sphere)
Shape='ellipseZ'; 
Sns = SurfaceSph(shape_gallery(p,Shape)); 
Xns = reshape(Sns.cart.to_array,[],3);
% Laplace SL, S' and DL
[SMns,Spns,DMns] = kernelDLap(Sns);

np = 2*p*(p+1); sp=(p+1)^2;
ii = (1:sp)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;
[u,v]=gl_grid(p);

Y = zeros(np,sp); SYsp = Y; SYns = Y;  Y2=Y; 
eigval=zeros(sp,1);

for k=1:sp
    n = nn(k); m = mm(k); 
    Y(:,k) = Ynm(n,m,u,v); 
    fac=(1./sqrt(uu^2-cos(u).^2)); fac=fac(:); 
    Y2(:,k) = fac.*Y(:,k);

    [P,Q]=ALF(n,uu);      
    Pnm=(factorial(n+m)./factorial(n-m)*(m<0)+(m>=0)).*P(abs(m)+1,:);
    Qnm=(factorial(n+m)./factorial(n-m)*(m<0)+(m>=0)).*Q(abs(m)+1,:);
    eigval(k)=spectrum2(n,m,uu);
end

if p<3
sprintf("Y matrix:")
disp(Y)
pause;
end

if strcmp(Mat,'SL')
    SYsp = SMsp*Y; 
    SYns = SMns*Y2;
    pot='SL_L_3D';

    SpheroidL=spheroidalSL(Y,uu);
elseif strcmp(Mat,'DL')
    SYsp = DMsp*Y; 
    SYns = DMns*Y;
    pot='DL_L_3D';

    SpheroidL=spheroidalDL(Y,uu);
else
    SYsp = Spsp*Y; 
    SYns = Spns*Y;
    pot='dSL_L_3D';
end

SYspsh = shAna(SYsp); 
SYnssh = shAna(SYns);
% SYnssh = spheroidalAna(SYns);
SpheroidLsh=spheroidalAna(SpheroidL); %Leo's spheroidal formula code

if p<3
sprintf("Ynm Coefficients after Laplace Kernel is applied")
disp(SYnssh)
pause;
% disp(SpheroidL)
sprintf("Ynm Coefficients after Spheroidal Formula applied")
disp(SpheroidLsh)
pause;
end


if verb
tol=1e-14;
SYspsh2 = zeros(sp,sp);
SYnssh2 = zeros(sp,sp);
SYnssh3 = SYnssh2; 

SpheroidLsh2 = zeros(sp,sp);

for i=1:sp
    for j=1:sp
        if abs(SYspsh(i,j))>tol
            SYspsh2(i,j)=SYspsh(i,j);
        end
        if abs(SYnssh(i,j))>tol
            SYnssh2(i,j)=SYnssh(i,j);
        end
        if abs(SYnssh(i,j))>1e-5
            SYnssh3(i,j)=SYnssh(i,j);
        end
        if abs(SpheroidLsh(i,j))>tol
            SpheroidLsh2(i,j)=SpheroidLsh(i,j);
        end
    end
end


% eigsh=diag(SpheroidLsh2);

eignm=diag(SYnssh);
eigdif=real(eigsh-eignm);

eigerr=abs(eigsh-eignm);
err1=abs(SpheroidL-SYns);

figure;
plot(1:length(eigerr),log10(abs(eigsh-eignm)))
% plot(1:length(eigerr),log10(abs(eigval-eignm)))
title(sprintf("eigenvalue log10 error, p=%d",p));

mmax=mm(end)-1;
m0=-1;
nmax=nn(end);

eigdif_nmax=zeros(2*mmax+1,1);
eigdif_m0=zeros(nmax+1-abs(m0),1);

m_nmax=zeros(2*mmax+1,1);
n_m0=zeros(nmax+1-abs(m0),1);

a1=1;
a2=1;
fprintf("Numerical vs Spheroidal Harmonics Formula\n")
for k=1:sp
    n = nn(k); m = mm(k); 
%     if n==nmax && abs(m)<= mmax
%         eigdif_nmax(a1)=eigdif(k);
%         eigdiv_nmax(a1)=eigdiv(k);
%         m_nmax(a1)=m;
%         a1=a1+1;
%     end
% 
%     if m==m0
%         eigdif_m0(a2)=eigdif(k);
%         eigdiv_m0(a2)=eigdiv(k);
%         n_m0(a2)=n;
%         a2=a2+1;
%     end
% 
%     if n==fixn
%         m0err(l)=eigerr(sp);
%     end

    fprintf("lambda_{%d}^{%d}:  %.5f  %.5f\n",n,m,eignm(k),eigsh(k))
end

[~,j1] = sort(mm);

% figure;
% plot(m_nmax,log10(abs(eigdif_nmax)))
% title(sprintf("eigenvalue log10 error, n=%d",nmax));
% xlabel("m")

figure;
% spy(SpheroidLsh);
spy(SpheroidLsh2);
title("Spheroidal Harmonic Transform")

figure; 
spy(SYspsh2);
title("sphere")
figure; 
spy(SYnssh2(j1,j1));
title("non sphere")
figure; 
spy(SYnssh3(j1,j1));
title("non sphere >1e-5")
figure; 
imagesc(log10(abs(SYnssh2(j1,j1))))
colormap(hsv(512))
colorbar
axis equal;
end
end

% figure;
% plot(parr,log10(m0err))
% title(sprintf("eigenvalue log10 error, m=0, n=%d",fixn))
% xlabel("p")
% 
% pause; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Far-evaluation check

%Kernel parameters
par = Kernel_Eval_parameters(pot,0,1,1,1,1e-8,2,400,1);
par.dim = 3;
[~, gwt]=g_grid(p+1);
wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
wt = wt(:);
Wsp = Sph.geoProp.W; Wsp= Wsp.*wt;
Wns = Sns.geoProp.W; Wns= Wns.*wt;

Nr = reshape(Sph.geoProp.nor.to_array,[],3); 
Nrns = reshape(Sns.geoProp.nor.to_array,[],3);

par.X = Xsp; par.nor = Nr; par.W2 = Wsp.'; 
pns=par;
pns.X = Xns; pns.nor = Nrns; pns.W2 = Wns.';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Target points outside sphere
%rt=4; 
Xtsp = rt*Xsp;
Rtsp = sqrt(sum(Xtsp'.*Xtsp'))'; 

%Compare kernel eval with spherical harmonic expansion
KYsp = Kernel_Eval(Xtsp,Xsp,par)*Y;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Target points outside sphere
Xtns = rt*Xns; %Xns + (rt-1)*Nrns;
if strcmp(Shape,'ellipseZ')
    display(max(Xtns))
    u = 2/sqrt(3); a = sqrt(3)/2;  
    b = @(x) a*sqrt(x^2-1);
    Xtns(:,3) = Xns(:,3)*rt;
    Xtns(:,1:2) = Xns(:,1:2)*2*b(rt*u);
    display(max(Xtns))
end
%Rtns = sqrt(sum(Xtns'.*Xtns'))'; 

%Compare kernel eval with spherical harmonic expansion
KYns = Kernel_Eval(Xtns,Xns,pns)*Y2;
KYns2 = Kernel_Eval(Xtsp,Xns,pns)*Y2;

KYspsh = shAna(KYsp); 
KYnssh = shAna(KYns);
KYnssht = shAna(KYns2);

if verb
tol=1e-14;
KYspsh2 = zeros(sp,sp);
KYnssh2 = zeros(sp,sp);
KYnssht2 = zeros(sp,sp);

for i=1:sp
    for j=1:sp
        if abs(KYspsh(i,j))>tol
            KYspsh2(i,j)=KYspsh(i,j);
        end
        if abs(KYnssh(i,j))>tol
            KYnssh2(i,j)=KYnssh(i,j);
        end
        if abs(KYnssht(i,j))>tol
            KYnssht2(i,j)=KYnssht(i,j);
        end
    end
end

% figure; 
% spy(KYspsh2); 
% figure; 
% spy(KYnssh2(j1,j1)); 
% figure; 
% imagesc(log10(abs(KYnssh2(j1,j1))))
% colormap(hsv(512))
% colorbar
% axis equal; 

figure;
plot(1:length(eigval),log10(abs(eigval-eigsh)))

end

close all

% Select range of p-values to test
pstart=12; pend=16;
parr=pstart:pend;

% Matrix to store error between KernelEval smooth quadrature and Double
% Layer formula
Dens_err=zeros(length(parr),2*pend*(pend+1));

for i=1:length(parr)
p=parr(i); Mat='DL';

% Test surf (non-sphere)
Shape='ellipseZ'; 
Sns = SurfaceSph(shape_gallery(p,Shape)); 
Xns = reshape(Sns.cart.to_array,[],3);
p_target=8;
Stns = SurfaceSph(shape_gallery(p_target,Shape)); 
Xtns_tmp = reshape(Stns.cart.to_array,[],3);

np = 2*p*(p+1); sp=(p+1)^2;
ii = (1:sp)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;
[u,v]=gl_grid(p);

Y = zeros(np,sp); Y2=Y; 
uu = 2/sqrt(3);

for k=1:sp
    n = nn(k); m = mm(k); 
    Y(:,k) = Ynm(n,m,u,v); 
    fac=(1./sqrt(uu^2-cos(u).^2)); fac=fac(:); 
    Y2(:,k) = fac.*Y(:,k);
end

if strcmp(Mat,'SL')
    pot='SL_L_3D';
elseif strcmp(Mat,'DL')
    pot='DL_L_3D';
else
    pot='dSL_L_3D';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Far-evaluation check

[~,j1] = sort(mm);

%Kernel parameters
pns = Kernel_Eval_parameters(pot,0,1,1,1,1e-8,2,400,1);
pns.dim = 3;
[~, gwt]=g_grid(p+1);
wt = pi/p*repmat(gwt', 2*p, 1)./sin(gl_grid(p));
wt = wt(:);
Wns = Sns.geoProp.W; Wns= Wns.*wt;
Nrns = reshape(Sns.geoProp.nor.to_array,[],3);
pns.X = Xns; pns.nor = Nrns; pns.W2 = Wns.';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Target points outside sphere
c=[2,3,3]; thetax=pi/3; thetay=-pi/10;
Rx=[1,0,0;0,cos(thetax),-sin(thetax);0,sin(thetax),cos(thetax)];
Ry=[cos(thetay),0,sin(thetay);0,1,0;-sin(thetay),0,cos(thetay)];
Xtns=zeros(size(Xtns_tmp));
for i=1:2*p_target*(p_target+1)
    Xtns(i,:)=((Xtns_tmp(i,:)-c)*Rx')*Ry';
end
if strcmp(Shape,'ellipseZ')
    display(max(Xtns))
    u0 = 2/sqrt(3); a = sqrt(3)/2;  
end


%Just calculate eigenvalues from formula
spheroidal_coord_targets=cart2spheroidal(Xtns,a);
ut=spheroidal_coord_targets(:,1);
% DLspectra=zeros(sp,1);
% for k=1:sp
%     n = nn(k); m = mm(k);
%     DLspectra(k)=DLspectrum(n,m,u0,ut(k));
% end


% Test surface density, with coefficients ...
sigma=Y(:,3);

%Calculate using spheroidalDL function
myDL=spheroidalDL(Y,u0,repmat(Xtns,1,sp));


%Compare kernel eval with spherical harmonic expansion
KYns = Kernel_Eval(Xtns,Xns,pns)*Y; %DL
% KYns = Kernel_Eval(Xtns,Xns,pns)*Y2; %SL
KYnssh = shAna(KYns);

%----------------------------------------------------------------
myDL=spheroidalDL(sigma,u0,Xtns);
KYns = Kernel_Eval(Xtns,Xns,pns)*sigma; %DL
KYnssh = shAna(KYns);
%----------------------------------------------------------------


figure; semilogy(abs(myDL-KYns)); 
title(sprintf('Function Error, p=%d',p))

figure; semilogy(abs(shAna(myDL)-KYnssh));
title(sprintf('Coefficients Error, p=%d', p))

% tol=1e-14;
% KYnssh2 = zeros(sp,sp);
% 
% for i=1:sp
%     for j=1:sp
%         if abs(KYnssh(i,j))>tol
%             KYnssh2(i,j)=KYnssh(i,j);
%         end
%     end
% end

% this_err=(abs(myDL-KYns)./abs(KYns))';
% Dens_err(i,1:np)=this_err;

% DLerr=log10(abs(KYnssh-diag(DLspectra)));
% figure;
% imagesc(DLerr(j1,j1))
% colorbar
% title(sprintf('DL: Off-surface kernel eval vs formulas. p=%d, rt=%.3f',p,rt))

% figure; 
% spy(KYnssh2(j1,j1)); 

% figure; 
% imagesc(log10(abs(KYnssh2(j1,j1))))
% colormap(hsv(512))
% colorbar
% axis equal; 
% title(sprintf('p=%d',p))

%{
figure; 
spy(KYnssht2(j1,j1)); 
figure; 
imagesc(log10(abs(KYnssht2(j1,j1))))
colormap(hsv(512))
colorbar
axis equal;
%}
end

% figure;
% imagesc(log10(Dens_err))
% colorbar
% xlabel('p')
% ylabel('target point')


%%

%timer 

fprintf('\n')
fprintf('%d surface densities, one target: ',sp)
tic
myDLtest1=spheroidalDL(Y,u0,repmat(Xtns(1,:),1,sp));
toc

[nt,~]=size(Xtns);
fprintf('one surface density, %d targets: ',nt)
tic
myDLtest2=spheroidalDL(Y(:,2),u0,Xtns);
toc

fprintf('%d surface densities, %d targets: ',sp,nt)
tic
myDLtest2=spheroidalDL(Y,u0,repmat(Xtns,1,sp));
toc

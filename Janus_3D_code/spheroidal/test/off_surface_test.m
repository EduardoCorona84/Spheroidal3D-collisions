% function [KYnssh,KYns,Y] = Test_Spharm_Lap_nsph(p,Shape,rt,Mat,verb)
p=16; Mat='DL'; verb=true; 

% Test surf (non-sphere)
Shape='ellipseZ'; 
Sns = SurfaceSph(shape_gallery(p,Shape)); 
Xns = reshape(Sns.cart.to_array,[],3);

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
rt=2; 
Xtns = rt*Xns; %Xns + (rt-1)*Nrns;
if strcmp(Shape,'ellipseZ')
    display(max(Xtns))
    u0 = 2/sqrt(3); a = sqrt(3)/2;  
    b = @(x) a*sqrt(x^2-1);
    Xtns(:,3) = Xns(:,3)*rt;
    Xtns(:,1:2) = Xns(:,1:2)*2*b(rt*u0);
    display(max(Xtns))
end
%Rtns = sqrt(sum(Xtns'.*Xtns'))'; 

%Just calculate eigenvalues from formula
spheroidal_coord_targets=cart2spheroidal(Xtns,a);
ut=spheroidal_coord_targets(:,1);
% DLspectra=zeros(sp,1);
% for k=1:sp
%     n = nn(k); m = mm(k);
%     DLspectra(k)=DLspectrum(n,m,u0,ut(k));
% end

%Calculate using spheroidalDL function
profile on
myDL=spheroidalDL(Y,u0,repmat(Xtns,1,sp));
prof=profile('info');
profile off
profview(0,prof)


%Compare kernel eval with spherical harmonic expansion
KYns = Kernel_Eval(Xtns,Xns,pns)*Y; %DL
% KYns = Kernel_Eval(Xtns,Xns,pns)*Y2; %SL
KYnssh = shAna(KYns);

figure; imagesc(log10(abs(myDL-KYns))); colorbar;
title('Function Error')

figure; imagesc(log10(abs(shAna(myDL)-KYnssh))); colorbar;
title('Coefficients Error')

if verb
tol=1e-14;
KYnssh2 = zeros(sp,sp);

for i=1:sp
    for j=1:sp
        if abs(KYnssh(i,j))>tol
            KYnssh2(i,j)=KYnssh(i,j);
        end
    end
end

% DLerr=log10(abs(KYnssh-diag(DLspectra)));
% figure;
% imagesc(DLerr(j1,j1))
% colorbar
% title(sprintf('DL: Off-surface kernel eval vs formulas. p=%d, rt=%.3f',p,rt))

figure; 
spy(KYnssh2(j1,j1)); 
figure; 
imagesc(log10(abs(KYnssh2(j1,j1))))
colormap(hsv(512))
colorbar
axis equal; 

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

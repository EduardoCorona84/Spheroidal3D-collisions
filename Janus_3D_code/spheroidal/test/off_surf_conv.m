% Select range of p-values to test
pstart=2; pend=4;
parr=pstart:pend;

% Matrix to store error between KernelEval smooth quadrature and Double
% Layer formula
Dens_err=zeros(length(parr),2*pend*(pend+1));

% Loop over p
for i=1:length(parr)

p=parr(i); Mat='DL';

% Test surf (non-sphere)
Shape='ellipseZ'; 
Sns = SurfaceSph(shape_gallery(p,Shape)); 
Xns = reshape(Sns.cart.to_array,[],3);

np = 2*p*(p+1); sp=(p+1)^2;
ii = (1:sp)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;
[v,phi]=gl_grid(p);

Y = zeros(np,sp); Y2=Y; 
u0 = 2/sqrt(3);

for k=1:sp
    n = nn(k); m = mm(k); 
    Y(:,k) = Ynm(n,m,v,phi); 
    fac=(1./sqrt(u0^2-cos(v).^2)); fac=fac(:); 
    Y2(:,k) = fac.*Y(:,k);
end

% Test surface density, with coefficients all equal to 1.
sigma=sum(Y,2);

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
% Choose target points as points on other spheroid surface 
c=[2,3,3]; thetax=pi/3; thetay=-pi/10;
Rx=[1,0,0;0,cos(thetax),-sin(thetax);0,sin(thetax),cos(thetax)];
Ry=[cos(thetay),0,sin(thetay);0,1,0;-sin(thetay),0,cos(thetay)];
Xtns=zeros(size(Xns));
for i=1:np
    Xtns(i,:)=((Xns(i,:)-c)*Rx')*Ry';
end

if strcmp(Shape,'ellipseZ')
    u0 = 2/sqrt(3); a = sqrt(3)/2;  
    b = @(x) a*sqrt(x^2-1);
end

%check that spheroids are not intersecting:
spheroidal_coord_targets=cart2spheroidal(Xtns,a);
ut=spheroidal_coord_targets(:,1);
%u coordinates should all be greater than u_0
if sum(ut<=u0)>0
    error('error: spheroids are intersecting.')
end

%Just calculate eigenvalues from formula
DLspectra=zeros(np,sp);
for k=1:np
    for l=1:sp
        n = nn(l); m = mm(l);
%         DLspectra(k,l)=DLspectrum(n,m,u0,ut(k));
    end
end


%Calculate using spheroidalDL function
myDL=spheroidalDL(sigma,u0,Xtns);

%Compare kernel eval with spherical harmonic expansion
KYns = Kernel_Eval(Xtns,Xns,pns)*sigma; %DL
% KYns = Kernel_Eval(Xtns,Xns,pns)*Y2; %SL
KYnssh = shAna(KYns);

figure; imagesc(log10(abs(myDL-KYns))); colorbar;
title(sprintf('Function Error, p=%d',p))

figure; imagesc(log10(abs(shAna(myDL)-KYnssh))); colorbar;
title(sprintf('Coefficients Error, p=%d', p))


figure;
semilogy(1:np,abs(myDL-KYns)./abs(KYns))
title(sprintf('DL: Off-surface kernel eval vs my DL (density). p=%d',p))

disp(np)
disp(size((abs(myDL-KYns)./abs(KYns))'))
this_err=(abs(myDL-KYns)./abs(KYns))';
Dens_err(i,1:np)=this_err;
disp(size(Dens_err))


figure; 
imagesc(log10(abs(KYnssh)))
colormap(hsv(512))
colorbar
axis equal; 
title(sprintf('p=%d',p))

end


%%

figure;
imagesc(log10(Dens_err))
colorbar
xlabel('p')
ylabel('target point')

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

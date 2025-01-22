% dS spectrum tests
% Authors: L Crowder, T Li
% reference: SLtest.m by L Crowder, DLtest_oblate.m by T Li
clear;
pstart=2; pend=16; 
parr=pstart:pend; spend=(pend+1)^2;
uu=2/sqrt(3);
a=1./uu;
% a=1./sqrt(uu^2+1);

pars=SpheroidalParameters;
pars.u0 = uu;
pars.a = a;
pars.oblate=0;
f = @(u,v) exp(sin(u).^4.*cos(u).*cos(4*v));

%------------------------------------------------------------------------%
%{
mySPY_cell=cell(1,length(parr));
mySPYsh_cell=cell(1,length(parr));
SPYsh_cell=cell(1,length(parr));
SPY_cell=cell(1,length(parr));

%------------------------------------------------------------------------%
for l=1:length(parr)
    
    p=parr(l);

    np = 2*p*(p+1); sp=(p+1)^2;
    ii = (1:sp)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;
    [u,v]=gl_grid(p);
    
    Y = zeros(np,sp); 
    for k=1:sp
        n = nn(k); m = mm(k); 
        Y(:,k) = Ynm(n,m,u,v); 
    end

    pars.sigma=Y;
    mySPY=spheroidalSP(pars);
    mySPYsh=shAna(mySPY);
    [SPYsh,SPY]=SDY(p,'SP');
    
    mySPY_cell{l}=mySPY;
    mySPYsh_cell{l}=mySPYsh;
    SPYsh_cell{l}=SPYsh;
    SPY_cell{l}=SPY;

end

%%
close all

err = log10(abs(mySPYsh_cell{end}-SPYsh_cell{end}));

figure;
imagesc(err)
colorbar
title(['SHC error, p=',num2str(parr(end))])

figure;
imagesc(log10(abs(mySPY_cell{end}-SPY_cell{end})))
colorbar
title("S'[Y] error")

p_index=4;
figure;
imagesc(log10(abs(mySPYsh_cell{p_index}-SPYsh_cell{p_index})))
colorbar
title(['SHC error, p=',num2str(parr(p_index))])


np = 2*p*(p+1); sp=(p+1)^2;
ii = (1:sp)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;
[~,j_sort] = sort(mm);


err = err(j_sort,j_sort);

figure;
imagesc(err)
colorbar
title(['SHC error, sorted by m. p=',num2str(parr(end))])


spy_tol = -5;
spy_err=err> spy_tol;

figure;
spy(spy_err)
title(['SHC error, p=',num2str(parr(end))])
%}

%% on surface
%{
err_norm=zeros(size(parr));
for l=1:length(parr)
    p=parr(l);
    [u,v]=gl_grid(p);
    Y=f(u,v);
    pars.sigma=Y;
    mySPY=spheroidalSP(pars);
    [~,SPY]=SDY(p,'SP',0,Y);
    
    err_norm(l)=norm(abs(SPY-mySPY));
end

figure; 
semilogy(parr,err_norm);
xlabel("p"); ylabel("spectra - SDY");
% title("On Surface oblate SP spectra");
%}

%% off surface
% %{
dists=[2,3,4];
dists2=[8,9,10];
[u,v]=gl_grid(pend);
Yp = f(u,v);

uu_trg=1.3*uu;
a_trg=1./uu_trg;
uu_trg2=1.1*uu; a_trg2=1./uu_trg2;

pot='dSL_L_3D';
KEparams = Kernel_Eval_parameters(pot,0,1,1,1,1e-8,2,400,1);
KEparams.dim = 3;
% Xself=prolate_spheroid_shape(pend,uu,a);
Xself=oblate_spheroid_shape(pend,uu,a);
Sns = SurfaceSph(Xself);
KEparams.X = Xself;
[~, gwt]=g_grid(pend+1);
wt = pi/pend*repmat(gwt', 2*pend, 1)./sin(gl_grid(pend));
wt = wt(:);
Wns = Sns.geoProp.W; Wns= Wns.*wt;
KEparams.W2 = Wns.';

Xtrg1=prolate_spheroid_shape(pend,uu_trg,a_trg);
Xtrg2=prolate_spheroid_shape(pend,uu_trg2,a_trg2);

errnorm=cell(2,length(dists));

for d=1:length(dists)
    % Whole Target Sphere
    Xtrg_dist1=Xtrg1+dists(d).*ones(size(Xtrg1));
    Xtrg_dist2=Xtrg2+dists2(d).*ones(size(Xtrg2));
    X={Xtrg_dist1,Xtrg_dist2};
    Strg1=SurfaceSph(Xtrg_dist1);
    nu_cart = reshape(Strg1.geoProp.nor.to_array,[],3);
    Strg2=SurfaceSph(Xtrg_dist2);
    nu_cart2 = reshape(Strg2.geoProp.nor.to_array,[],3);
    
    % % DEBUG component-wise normal
    % Xtrg_dist=Xtrg_dist(100:102,:);
    % nu_sph = [1,0,0;0,1,0;0,0,1];
    % strg=cart2spheroidal(Xtrg_dist,a_trg);
    % nu_cart = spheroidalNu2cart(nu_sph,strg,a);
    % KEparams.nor = nu_cart;

    KEparams.nor=nu_cart;
    dSmat1 = Kernel_Eval(Xtrg_dist1,Xself,KEparams);
    calcDY1 = dSmat1 * Yp;
    KEparams.nor=nu_cart2;
    dSmat2 = Kernel_Eval(Xtrg_dist2,Xself,KEparams);
    calcDY2 = dSmat2 * Yp;

    % % DEBUG component-wise normal
    % err_p=zeros(3,length(parr));
    % conv_p = zeros(3,length(parr));

    % Whole Sphere
    err_p1=zeros(1,length(parr));
    err_p2=zeros(1,length(parr));

    for k=1:length(parr)
        p = parr(k);
        [u,v]=gl_grid(p);
        Y = f(u,v); 
        Ysig=cat(3,Y,Y);
        pars.sigma= Ysig;
        nu=cat(3,nu_cart,nu_cart2);
        myDY=spheroidalSP(pars,X,nu);

        % Whole Sphere
        err_p1(:,k) = norm(abs(myDY{1}-calcDY1));
        err_p2(:,k) = norm(abs(myDY{2}-calcDY2));

        % % DEBUG component-wise normal
        % err_p(:,k) = abs(myDY-calcDY); 
        % conv_p(:,k) = abs(myDY-myDYend);
    end
    errnorm{1,d}=err_p1;
    errnorm{2,d}=err_p2;
end

% Whole Sphere
figure;
for l=1:length(dists)
    err_p = cell2mat(errnorm(1,l));
    semilogy(parr,err_p); hold on;
end
hold off;
xlabel("p"); ylabel("norm error");
title("difference from Kernel Eval")

figure;
for l=1:length(dists)
    err_p = cell2mat(errnorm(2,l));
    semilogy(parr,err_p); hold on;
end
hold off;
xlabel("p"); ylabel("norm error");
title("difference from Kernel Eval")

% % DEBUG component-wise normal
% for l=1:length(dists)
%     figure;
%     err_p = cell2mat(errnorm(1,d));
%     for kk=1:3
%         semilogy(parr,err_p(kk,:)); hold on;
%     end
%     hold off;
%     xlabel("p"); ylabel("norm error");
%     legend("u dir","v dir","phi dir")
%     title("diff from Kernel Eval, different spheroidal component")
% end


% %}


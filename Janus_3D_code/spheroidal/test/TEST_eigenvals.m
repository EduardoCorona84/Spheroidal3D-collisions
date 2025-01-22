% Numerically analyze possible eigenvalues of I/2+D+etaS*S
% (I/2+D+etaS*S)Ynm = 
% ( 0.5 + lambdaDnm + etaS * <Ynm,Ylk/sqrt(u0^2-v^2)> * lambdaSnm ) Ynm

% See scaled_completion.m for tests on scaling etaS to minimize condition
% number.

% ---- Set up --------------------------------------------------
clear; close all;
p = 16; sp=(p+1)^2;
ii = (1:sp)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;

ind_n0 = zeros(1,p+1);
for n=0:p
    if n==0
        ind_n0(n+1) = 1;
    else
        ind_n0(n+1) = ind_n0(n)+2*n;
    end
end
nonind_n0 = setdiff(1:numel(nn),ind_n0);

% R = 2.^linspace(0.1,6.01,30);
R = [1.1,4,16,50];
eps = 1./R;
u0=1./sqrt(1-1./R.^2);
a = 1./u0;

anm=factorial(nn-mm)./factorial(nn+mm).*(-1).^mm .*(u0.^2-1);
bnm = a.*factorial(nn-mm)./factorial(nn+mm).*(-1).^mm .*sqrt(u0.^2-1);
L=legendre_otc(p,u0,1,1,1);

P=L{1}; Q=L{2}; dP=L{3}; dQ=L{4};

lambdaD=anm./2.*(P.*dQ+dP.*Q);
lambdaS=bnm.*P.*Q; % both size (sp x #ofR)
% size(lambdaS)

lambdaDsorted = [lambdaD(ind_n0,:);lambdaD(nonind_n0,:)];
lambdaSsorted = [lambdaS(ind_n0,:);lambdaS(nonind_n0,:)];

for eind=1:length(R)
    this_u0 = u0(eind);

    % lambdaD > -0.5, and lambdaD00 = -0.5:
    figure(1);
    % plot(1:length(lambdaD(:,eind)),lambdaD(:,eind),'-*'); hold on;
    plot(1:length(lambdaD(:,eind)),lambdaDsorted(:,eind),'-*'); hold on;

    % lambdaS >= 0, decays
    figure(2);
    % plot(1:length(lambdaS(:,eind)),lambdaS(:,eind),"-*"); hold on;
    plot(1:length(lambdaS(:,eind)),lambdaSsorted(:,eind),'-*'); hold on;

    % coefficient matrix, diagonal dominates
    figure(3)
    subplot(1,length(R),eind)
    G = Gmatrix(p,this_u0); % from Ynm/sqrt() to Ynm
    imagesc(G)
    colorbar

    diagjj = diag(G);
    diag_dom = abs(diagjj) >= sum(abs(G),2)-abs(diagjj);
    fprintf("G: number of rows where |aii|<sum_{j~=i}|aij| %d \n", sum(~diag_dom))
    fprintf("G: min entry is %f \n", min(min(G)))

    % inverse of G, from Ynm to Ynm/sqrt(u0^2-v^2)
    Ginv = eye(size(G,1))\G;
    figure(4)
    subplot(1,length(R),eind)
    imagesc(Ginv)
    colorbar

    diagii = diag(Ginv);
    diag_dom_inv = abs(diagii) >= sum(abs(Ginv),2)-abs(diagii);
    fprintf("G inv: number of rows where |aii|<sum_{j~=i}|aij| %d \n", sum(~diag_dom_inv))
    fprintf("G inv: min entry is %f \n", min(min(Ginv)))

    % Coefficient in Ynm basis for S[Ynm]
    SY_coeffs = (Ginv'.*lambdaSsorted(:,eind))';
    figure(5)
    subplot(1,length(R),eind)
    imagesc(SY_coeffs(1:80,1:80))
    colorbar

    diagkk = diag(SY_coeffs);
    diag_dom_SY = abs(diagkk) >= sum(abs(SY_coeffs),2)-abs(diagkk);
    fprintf("S[Y] coeffs: number of rows where |aii|<sum_{j~=i}|aij| %d \n", sum(~diag_dom_SY))
    fprintf("S[Y]: min entry is %f \n", min(min(SY_coeffs)))
    
end
figure(1)
% xlabel("n,m")
xlabel("m=0 first, all others next")
ylabel("lambdaD")
legend("R=1.1","R=4","R=16","R=50")
hold off;

figure(2);
% xlabel("n,m")
xlabel("m=0 first, all others next")
ylabel("lambdaS")
legend("R=1.1","R=4","R=16","R=50")
hold off;

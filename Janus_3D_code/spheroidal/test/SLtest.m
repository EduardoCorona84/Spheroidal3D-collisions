pstart=2; pend=16; 
parr=pstart:pend; spend=(pend+1)^2;
u0=2/sqrt(3);
a=1./u0;

pars=SpheroidalParameters;
pars.u0 = u0;
pars.a = a;

[~,spectra_surf,~]=SLspectrum(pend,u0,a);

%------------------------------------------------------------------------%

mySY_cell=cell(1,length(parr));
mySYsh_cell=cell(1,length(parr));
SYsh_cell=cell(1,length(parr));
SY_cell=cell(1,length(parr));

%------------------------------------------------------------------------%
for l=1:length(parr)
    
    p=parr(l);

    np = 2*p*(p+1); sp=(p+1)^2;
    ii = (1:sp)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;
    [u,v]=gl_grid(p);
    
    Y = zeros(np,sp); SYns = Y; Y2=zeros(np,sp);
    
    for k=1:sp
        n = nn(k); m = mm(k); 
        Y(:,k) = Ynm(n,m,u,v); 
        fac=(1./sqrt(u0^2-cos(u).^2)); fac=fac(:); 
        Y2(:,k) = fac.*Y(:,k);
    end
    
    pars.sigma=Y2;

    mySY=spheroidalSL(pars);
    mySYsh=shAna(mySY);
    [SYsh,SY]=SDY(p,'SL');
    
    mySY_cell{l}=mySY;
    mySYsh_cell{l}=mySYsh;
    SYsh_cell{l}=SYsh;
    SY_cell{l}=SY;

end

%%
close all

[~,j1]=sort(mm);
mySYsh_end=mySYsh_cell{end};
SYsh_end=SYsh_cell{end};
mySY_end=mySY_cell{end};
SY_end=SY_cell{end};

figure;
% imagesc(log10(abs((SYsh_end-diag(spectra_surf))./SYsh_end)))
imagesc(log10(abs((SYsh_end-diag(spectra_surf)))))
colorbar
title('formulas')

figure;
% imagesc(log10(abs(mySYsh_end(:,j1)-SYsh_end(:,j1))))
imagesc(log10(abs(mySYsh_end-SYsh_end)))
colorbar

figure;
imagesc(log10(abs(mySY_end-SY_end)))
colorbar

nm_err=zeros(pend+1,2*pend+1);
nm_err2=nm_err;
SYsh_diag=diag(SYsh_end);
mySYsh_diag=diag(mySYsh_end);
for i=1:spend
    n=nn(i); m=mm(i); 
    nm_spec=spectra_surf(i);
    nm_SYsh=SYsh_diag(i);
%     nm_err(n+1,pend+m+1)=log10(abs(nm_spec-nm_SYsh));
    nm_err(n+1,pend+m+1)=log10(abs((nm_spec-nm_SYsh)./nm_SYsh));
    nm_err2(n+1,pend+m+1)=log10(abs((nm_spec-mySYsh_diag(i))./mySYsh_diag(i)));
end

figure;
imagesc(nm_err)
title('singular quad vs formula')
ylabel('n')
yticks(1:pend+1)
yticklabels(0:pend)
xlabel('m')
xticks(1:2*pend+1)
xticklabels(-pend:pend)
colorbar

figure;
imagesc(nm_err2)
title('my SL vs formula')
ylabel('n')
yticks(1:pend+1)
yticklabels(0:pend)
xlabel('m')
xticks(1:2*pend+1)
xticklabels(-pend:pend)
colorbar

% figure;
% semilogy(1:spend,abs(SYsh_diag-spectra_surf))
% 
% figure; plot(1:spend,abs(mySYsh_diag./spectra_surf))

pstart=2; pend=16; 
parr=pstart:pend; spend=(pend+1)^2;
uu=2/sqrt(3);
a=1./uu;

pars=SpheroidalParameters;
pars.u0 = uu;
pars.a = a;

%------------------------------------------------------------------------%

myDY_cell=cell(1,length(parr));
myDYsh_cell=cell(1,length(parr));
DYsh_cell=cell(1,length(parr));
DY_cell=cell(1,length(parr));

%------------------------------------------------------------------------%
for l=1:length(parr)
    
    p=parr(l);

    np = 2*p*(p+1); sp=(p+1)^2;
    ii = (1:sp)'; nn = floor(sqrt(ii-1)); mm=ii-nn.^2-nn-1;
    [u,v]=gl_grid(p);
    
    Y = zeros(np,sp); DYns = Y; Y2=zeros(np,sp);
    myDL=zeros(np,sp);
    DL=zeros(np,sp);
    SLeigval=zeros(sp,1);
    
    for k=1:sp
        n = nn(k); m = mm(k); 
        Y(:,k) = Ynm(n,m,u,v); 
        fac=(1./sqrt(uu^2-cos(u).^2)); fac=fac(:); 
        Y2(:,k) = fac.*Y(:,k);
    end
    
    pars.sigma=Y;
    myDY=spheroidalDL(pars);
    myDYsh=shAna(myDY);
    [DYsh,DY]=SDY(p,'DL');
    
    myDY_cell{l}=myDY;
    myDYsh_cell{l}=myDYsh;
    DYsh_cell{l}=DYsh;
    DY_cell{l}=DY;


end

%%
close all

figure;
title(fprintf("p=%d",parr(end)))
imagesc(log10(abs(myDYsh_cell{end}-DYsh_cell{end})))
colorbar

% [i1,i2,i3]=size(myDY_cell{end});
% myDY_plot=reshape(myDY_cell{end},i1,i3);
% DY_plot=reshape(DY_cell{end},i1,i3);

figure;
imagesc(log10(abs(myDY_cell{end}-DY_cell{end})))
colorbar

figure;
title("p=4")
imagesc(log10(abs(myDYsh_cell{4}-DYsh_cell{4})))
colorbar



function [ Xtarg, log10relerror,x1]=BIEtest_funct(C,r,lambda)
hold off;

sp=@(p) 2*p*(p+1);
Pvals=[4 8 12 16 24];


%C=[1,1,0; -3,0,0; 0,0,4];
%C=[0 0 0];
np= @(p) 2*p*(p+1);
numsph=size(C,1);
numtargpoints=numsph*np(Pvals(5));
%r=[0.903,0.510,0.262]';
%r=[1];

%generates target points
params=Setparams(Pvals(5),C,r,'SL_LMOD_3D','',1,1,3,1,lambda);
x1=params.X;
n1=params.nor;
Xtarg=repmat(x1,6,1);
Ntarg=repmat(n1,6,1);
dfactor=(10).^(-(1:6)+2);
dfactor=repmat(reshape(repmat(dfactor,numtargpoints,1),[],1),1,3);
Xtarg=Xtarg+dfactor.*Ntarg;

%sets up point force 
rng('default');
delta=rand(numsph,3)-0.5*ones(numsph,3);
ndelta=sqrt(delta(:,1).^2+delta(:,2).^2+delta(:,3).^2);
delta=bsxfun(@times,delta,1./ndelta);
deltavec=bsxfun(@times,delta,rand(numsph,1).*(r/20));
xgreen=C+deltavec;
force=rand(numsph,1);


 

for j=1:size(Pvals,2);
    
    
%computes kernel
    p=Pvals(j);
    numpoints=np(p)*numsph;
    %generates kernel
    params=Setparams(p,C,r,'SL_LMOD_3D','',1,1,3,1,lambda);
    S=VSh_Mod_MatVec_RB2('Mat',[],params);
    paramsDL=params;
    paramsDL.a=1/2;
    paramsDL.flag_pot='DL_LMOD_3D';
    D=VSh_Mod_MatVec_RB2('Mat',[],paramsDL);
    K=S+D;
    P=make_projection(p);
    
    if size(C,1)==3
    P=blkdiag(P,P,P);   %<=== will have to modify if numsph~=3
    end
    KMOD=P*K*P+(eye(numpoints)-P); 


    %computes boundary values
    ui=zeros(numpoints,1);
    for i=1:numsph
        ui=ui+force(i)*modlapgreen(params.X,xgreen(i,:),lambda);
    end
    %solves equation

    mu=KMOD\(P*ui);
     f= @(u,n) VSh_Mod_MatVec_RB_trg(mu,u,n,params)+VSh_Mod_MatVec_RB_trg(mu,u,n,paramsDL);
     ut=zeros(numtargpoints*6,1);
    for i=1:numsph
        ut= ut + force(i)*modlapgreen(Xtarg,xgreen(i,:),lambda);
    end
    ut=reshape(ut,[],6);
    f=reshape(f(Xtarg,Xtarg),[],6);
    error=abs(ut-f);
    relerror=error./abs(ut);     %We can use this to get the box & whisker plot
    for k=1:6
        log10error(k)=log10(norm(error(:,k))/norm(abs(ut(:,k))));
    end
    results(:,j)=log10error;
    log10relerror(:,:,j)=log10(relerror);
   
   
end
 %vtkwrite(string, 'structured_grid', Xtarg(101:120,1),Xtarg(:,2),Xtarg(:,3),'scalars','error', log10relerror(:,6,j));

%make_graph(results,log10(Pvals));

end
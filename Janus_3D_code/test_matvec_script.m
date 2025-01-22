clear;clc;
sp =@(p) (p+1)^2;

type='SpMat';
lambda=1e-8;
rng('default');


np = @(p) 2*(p+1)*p; %(number of quadrature points);
ii=(1:sp(16))';
nn=floor(sqrt(ii-1)); %indices of n
%generate random seeds
swgt=rand(sp(16),1);




for i=1:sp(16)        %exponentially decaying weights
             swgt(i)=swgt(i)*(1/2)^(floor(sqrt(i-1)));
end

CS=[6 0 0];
CT=CS;

rs=[1];
nspheres=size(CS,1);
%numpoints=np*nspheres;






%{
ShV=zeros(sp(p),1);
ShV(1)=1;
V=shSyn(ShV);
%}
%V=repmat(V,10,1);
Sep=[2.25 2.5 2.75 3];
P=[4 8 12 16];
for pp=1:4
    for ss=1:4
        p=P(pp);
        s=Sep(ss);
        rng('default');
        clearvars swgt2;
        swgt2(1:sp(p-1))=swgt(1:sp(p-1));
        swgt2(sp(p-1)+1:sp(p))=0;
        %display(size(swgt2));
      %  V=repmat(shSyn(swgt2'),nspheres,1);
        ev=zeros((p+1)^2,1);
        ev(2)=1;
        V=shSyn(ev);
        rt=Sep(ss);
        %Cz2=Cz*Sep(ss);
        %CS=Cz2(1:10,:);
        %CT=Cz2(11:end,:);
        %rs=ones(size(CS,1),1);
        %rt=ones(size(CT,1),1);
        [K, KV,Error,KL,Allmod,AllLap]=test_matvec(p,type,CT,CS,rt,rs,lambda,V);
        Errorl=log10(norm(KV-(1/3)*(s)^(-2)*V)/norm(KV));
        Results(ss,pp)=Error;
        Resultsl(ss,pp)=Errorl;
    end
end

make_graph;
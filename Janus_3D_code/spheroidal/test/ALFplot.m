x=linspace(1.1,5,1000);
nmax=2;
num_fxns=(nmax+1)*(nmax+2)/2;

P=zeros([num_fxns,length(x)]);
Q=zeros([num_fxns,length(x)]);

for n=0:nmax
    ind_range= (n*(n+1)/2+1):(n+1)*(n+2)/2 ;
    [P(ind_range,:),Q(ind_range,:)]=ALF(n,x,1);
end

close all
figure()
hold on
for i=1:num_fxns
    n=floor((sqrt(8*i-7)-1)/2);
    m=i-(n*(n+1)/2+1);
    legendstr=sprintf('P_{%d}^{%d}',n,m);
    plot(x,P(i,:),'DisplayName',legendstr,'LineWidth',1)
end
hold off
legend

figure()
hold on
for i=1:num_fxns-1
    n=floor((sqrt(8*i-7)-1)/2);
    m=i-(n*(n+1)/2+1);
    legendstr=sprintf('Q_{%d}^{%d}',n,m);
    plot(x,Q(i,:),'DisplayName',legendstr,'LineWidth',1)
end
hold off
legend

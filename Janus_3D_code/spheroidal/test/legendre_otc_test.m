geti= @(n,m) m+n^2+n+1;

x=linspace(1.1,2,100);
maxp=5;

u0=2/sqrt(3);
Pp=cell(1,6);
Qp=cell(1,6);
for p=0:maxp
    PQ=legendre_otc(p,u0);
    Pp{p+1}=PQ{1};
    Qp{p+1}=PQ{2};
end

Pall=zeros((maxp+1)^2,maxp+1);
Qall=zeros((maxp+1)^2,maxp+1);


for p=0:maxp
    Ptmp=Pp{p+1};
    Qtmp=Qp{p+1};
    Pall(1:(p+1)^2,p+1)=Ptmp;
    Qall(1:(p+1)^2,p+1)=Qtmp;
end


p2=4;
PQ=legendre_otc(p2,x); P=PQ{1}; Q=PQ{2};

n=3; m=-1;


% Recurrence relation: (2n-1)x Q_{n-1}^m = (n-m)Q_n^m +(n+m-1)Q_{n-2}^m
EE=abs((n-m).*Q(geti(n,m),:)-(2*n-1).* x .* Q(geti(n-1,m),:)+(n+m-1).*Q(geti(n-2,m),:));
% Wronskian Relation: P_n^m Q_{n-1}^m - P_{n-1}^m Q_n^m =(-1)^m (n+m-1)!/(n-m)!
Eww=abs(P(geti(n,m),:).*Q(geti(n-1,m),:)-P(geti(n-1,m),:).*Q(geti(n,m),:)-factorial(n+m-1)./factorial(n-m).*(-1).^m);

close all
figure()
semilogy(x,EE) %good
title("Recurrence relation for Q over n: log error")

figure()
semilogy(x,Eww) %good
title("Wronskian relation over n: log error")


% Recurrence relation 
n=1;m=0;
Ep=sqrt(x.^2-1).*P(geti(n,m+1),:)-(n-m+1).*P(geti(n+1,m),:)+(n+m+1).*x.*P(geti(n,m),:);
figure;
semilogy(x,abs(Ep))
title(sprintf('P recurrence error, n=%d, m=%d',n,m))

Pnewvold=zeros(size(P));
Qnewvold=zeros(size(P));
for n=0:p2
    [Poldn,Qoldn]=ALF(n,x);
    for m=0:n
        Pnewvold(geti(n,m),:)=abs(Poldn(m+1,:)-P(geti(n,m),:));
        Qnewvold(geti(n,m),:)=abs(Qoldn(m+1,:)-Q(geti(n,m),:));
        if norm(Pnewvold(geti(n,m),:))>1e-5
            fprintf('P_%d^%d is wrong\n',n,m)
        end
        if norm(Qnewvold(geti(n,m),:))>1e-5
            fprintf('Q_%d^%d is wrong\n',n,m)
        end
    end
end
figure; imagesc(log10(Pnewvold)); colorbar; title('Pnewvold');
figure; imagesc(log10(Qnewvold)); colorbar; title('Qnewvold');



%%
close all
% Test derivatives
p3=16;
% x1=linspace(1.1,2,10);
x1=linspace(1.1,3,100);

L=legendre_otc(p3,x1,1,1,1);
P=L{1}; Q=L{2}; dP=L{3}; dQ=L{4};

Wrsnk=zeros(size(P));
for k=1:(p3+1)^2
    n = floor(sqrt(k-1)); m=k-n.^2-n-1;
%     fprintf('%d, %d\n',n,m)
%     disp(geti(n,m))
    Wrsnk(k,:)=abs(P(k,:).*dQ(k,:)-dP(k,:).*Q(k,:)-factorial(n+m).*(-1).^m./(factorial(n-m)*(1-x1.^2)) );
end

figure;
imagesc(log10(Wrsnk))
colorbar

Pcpp=readmatrix('~/Documents/spheroidal_cpp/data/P.txt');
Qcpp=readmatrix('~/Documents/spheroidal_cpp/data/Q.txt');
dPcpp=readmatrix('~/Documents/spheroidal_cpp/data/dP.txt');
dQcpp=readmatrix('~/Documents/spheroidal_cpp/data/dQ.txt');

Pe=log10(abs(Pcpp-P)./abs(P));
Qe=log10(abs(Qcpp-Q)./abs(Q));
dPe=log10(abs(dPcpp-dP)./abs(dP));
dQe=log10(abs(dQcpp-dQ)./(abs(dQ)));

bad=-5;
Pbad=Pe>bad; Qbad=Qe>bad; dPbad=dPe>bad; dQbad=dQe>bad;
Pkbad=max(Pbad,[],2); Qkbad=max(Qbad,[],2); 
dPkbad=max(dPbad,[],2); dQkbad=max(dQbad,[],2);
for k=1:(p3+1)^2
    n = floor(sqrt(k-1)); m=k-n.^2-n-1;
    if Pkbad(k)
        fprintf('P_%d^%d is bad\n',n,m);
    end
end
for k=1:(p3+1)^2
    n = floor(sqrt(k-1)); m=k-n.^2-n-1;
    if Qkbad(k)
        fprintf('Q_%d^%d is bad\n',n,m);
    end
end
for k=1:(p3+1)^2
    n = floor(sqrt(k-1)); m=k-n.^2-n-1;
    if dPkbad(k)
        fprintf('dP_%d^%d is bad\n',n,m);
    end
end
for k=1:(p3+1)^2
    n = floor(sqrt(k-1)); m=k-n.^2-n-1;
    if dQkbad(k)
        fprintf('dQ_%d^%d is bad\n',n,m);
    end
end

figure;
imagesc(Pe);
colorbar;
title('P')

figure;
imagesc(Qe);
colorbar;
title('Q')

figure;
imagesc(dPe);
colorbar;
title('dP')

figure;
imagesc(dQe);
colorbar;
title('dQ')

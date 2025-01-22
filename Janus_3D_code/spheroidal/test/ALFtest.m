n=2;
m=0;
x=linspace(1.1,2,100000);

[P0,Q0]=ALF(0,x);
[P1,Q1]=ALF(1,x);
[P2,Q2]=ALF(2,x);

% Recurrence relation: (2n-1)x Q_{n}^m = (n-m)Q_{n-1}^m +(n+m-1)Q_{n-2}^m
EE=abs((n-m).*Q2(m+1,:)-(2*n-1).* x .* Q1(m+1,:)+(n+m-1).*Q0(m+1,:));
% Wronskian Relation: P_n^m Q_{n-1}^m - P_{n-1}^m Q_n^m =(-1)^m (n+m-1)!/(n-m)!
Eww=abs(P2(m+1,:).*Q1(m+1,:)-P1(m+1,:).*Q2(m+1,:)-factorial(n+m-1)./factorial(n-m).*(-1).^m);


close all
figure()
plot(x,log10(EE)) %good
title("Recurrence relation for Q over n: log error")

figure()
plot(x,log10(Eww)) %good
title("Wronskian relation over n: log error")

%%

n=2;
m=1;

[P,Q]=ALF(n,x);

Ep=abs(P(m+2,:)+2*m*x./sqrt(x.^2-1).*P(m+1,:) - (n-m+1).*(n+m).*P(m,:));
Eq=abs(Q(m+2,:)+2*m*x./sqrt(x.^2-1).*Q(m+1,:) - (n-m+1).*(n+m).*Q(m,:));
Ew=abs(P(m+1,:).*Q(m+2,:)-P(m+2,:).*Q(m+1,:)-factorial(n+m)./factorial(n-m).*(-1).^(m+1)./sqrt(x.^2-1));

%asymptotic behavior as Pnm -> 1+, from DLMF
Pasym=factorial(n+m)./(factorial(m).*factorial(n-m)).*((x-1)./2).^(m/2);
Epasym=abs(Pasym-P(m+1,:));

m=-1;
mc0=factorial(n+m)./factorial(n-m);
mc1=factorial(n+m+1)./factorial(n-(m+1));
Ewmm=abs(mc0.*P(-m+1,:).*Q(-m,:)-P(-m,:).*mc0.*Q(-m+1,:)-factorial(n+m)./factorial(n-m).*(-1).^(m+1)./sqrt(x.^2-1));
% Ewmm=abs(P(-m+1,:).*Q(-m,:)-P(-m,:).*Q(-m+1,:)-factorial(n+m)./factorial(n-m).*(-1).^(m+1)./sqrt(x.^2-1));



close all

figure()
plot(x,log10(Ep)) %good
title("Recurrence relation for P over m: log error")

figure()
plot(x,log10(Eq)) %good
title("Recurrence relation for Q over m: log error")

figure()
plot(x,log10(Epasym)) %gets better as x -> 1+
title("log error of Asymptotic formula for P_n^m as x->1+")

figure()
plot(x,log10(Ew)) %good
title("Wronskian relation over m: log error")

figure()
plot(x,log10(Ewmm))  %good
title("Wronskian relation over m w/ m<0: log error")


%%

h=1e-8;
x0=1.1;
np=1e5; %number of points
x=x0:h/2:(x0+np*h);
xc=x(2:2:end-1);
xh=x(1:2:end);

close all

N=10;

dPm=zeros(1,2*N+1);
dQm=zeros(1,2*N+1);
Pm=zeros(1,2*N+1);
Qm=zeros(1,2*N+1);
rerr_P=zeros(1,2*N+1);
rerr_Q=zeros(1,2*N+1);
aerr_P=zeros(1,2*N+1);
aerr_Q=zeros(1,2*N+1);

for n=0:N
    [Pn,Qn]=ALF(n,xh);
    dPn_fd=(Pn(:,2:end)-Pn(:,1:end-1))./h;
    dQn_fd=(Qn(:,2:end)-Qn(:,1:end-1))./h;
    W_err=zeros(2*n+1);

    for m=-n:n
        Pnm=(factorial(n+m)./factorial(n-m)*(m<0)+(m>=0)).*Pn(abs(m)+1,:);
        Qnm=(factorial(n+m)./factorial(n-m)*(m<0)+(m>=0)).*Qn(abs(m)+1,:);

        [dPnm,dQnm]=dALF(n,m,xc);
        
        if n==N 
        dPm(m+N+1)=dPnm(1);
        dQm(m+N+1)=dQnm(1);
        Pm(m+N+1)=Pnm(1);
        Qm(m+N+1)=Qnm(1);
        end

        dPnm_fd=(factorial(n+m)./factorial(n-m)*(m<0)+(m>=0)).*dPn_fd(abs(m)+1,:);
        dQnm_fd=(factorial(n+m)./factorial(n-m)*(m<0)+(m>=0)).*dQn_fd(abs(m)+1,:);

        %relative error
        rerr_P(m+N+1)=abs((dPnm(1)-dPnm_fd(1))./dPnm_fd(1));
        rerr_Q(m+N+1)=abs((dQnm(1)-dQnm_fd(1))./dQnm_fd(1));
        %absolute error
        aerr_P(m+N+1)=abs(dPnm(1)-dPnm_fd(1));
        aerr_Q(m+N+1)=abs(dQnm(1)-dQnm_fd(1));

        W_err(m+n+1)=abs(Pnm(1).*dQnm(1)-dPnm(1).*Qnm(1)-factorial(n+m).*(-1).^m./(factorial(n-m)*(1-xc(1).^2)) );

%         if n==2
%         figure;
%         if n==0
%             plot(xc,abs(dPnm_fd-dPnm));
%             title(sprintf('dP_{%d}^{%d}(x) log abs error with finite difference',n,m))
%             xlabel("x")
%         else
%             plot(xc,log10(abs((dPnm_fd-dPnm)./dPnm)));
%             title(sprintf('dP_{%d}^{%d}(x) log relative error with finite difference',n,m))
%             xlabel("x")
%         end
% 
%         figure;
%         plot(xc,log10(abs((dQnm_fd-dQnm)./dQnm)));    
%         title(sprintf('dQ_{%d}^{%d}(x) log relative error with finite difference',n,m))
%         xlabel("x")
%         end
    
    end

    
    if n==0
        sprintf("n=0, Wronskian log error= %.6e", W_err)
    else
        figure;
        plot(-n:n,log10(W_err))
        title(sprintf('P_{%d}^{m}(%.1f): Wronskian log error',n,xc(1)))
        xlabel("m")
    end

end



% close all


figure;
hold on
plot(-N:N,log10(abs(Pm)))
plot(-N:N,log10(abs(Qm)))
plot(-N:N,log10(abs(dPm)))
plot(-N:N,log10(abs(dQm)))
hold off
legend(sprintf('log10|P_{%d}^{m}(%.1f)|',N,xc(1)), ...
    sprintf('log10|Q_{%d}^{m}(%.1f)|',N,xc(1)), ...
    sprintf('log10|dP_{%d}^{m}(%.1f)|',N,xc(1)), ...
    sprintf('log10|dQ_{%d}^{m}(%.1f)|',N,xc(1)))
xlabel("m")
title("scale of P, Q, and derivatives")

figure;
hold on
plot(-N:N,log10(abs(factorial(N-(-N:N))./factorial(N+(-N:N)).*dPm)))
plot(-N:N,log10(abs(factorial(N-(-N:N))./factorial(N+(-N:N)).*dQm)))
hold off
legend(sprintf('log10|(n-m)!/(n+m)! dP_{%d}^{m}(%.1f)|',N,xc(1)), ...
    sprintf('log10|(n-m)!/(n+m)! dQ_{%d}^{m}(%.1f)|',N,xc(1)))
xlabel("m")
title("scale of P, Q derivatives with factor")

figure;
plot(-N:N, log10(rerr_P))
title(sprintf("dP_{%d}^{m}(%.1f) log relative error with finite difference",N,xc(1)))
xlabel("m")

figure;
plot(-N:N, log10(rerr_Q))
title(sprintf("dQ_{%d}^{m}(%.1f) log relative error with finite difference",N,xc(1)))
xlabel("m")

figure;
plot(-N:N, log10(aerr_P))
title(sprintf("dP_{%d}^{m}(%.1f) log abs error with finite difference",N,xc(1)))
xlabel("m")

figure;
plot(-N:N, log10(aerr_Q))
title(sprintf("dQ_{%d}^{m}(%.1f) log abs error with finite difference",N,xc(1)))
xlabel("m")

%%

close all

N1=1;
N2=2;
xx=linspace(1.01,1.2,1000);

for n=N1:N2
    [dPn,dQn]=dALF(n,0,xx);
    dPnana=dPana(n,xx);
    err=abs(dPn-dPnana);
    figure;
    plot(xx,log10(err))
    title(sprintf("P error, n=%d",n))
    
    if n<=3
    dQnana=dQana(n,xx);
    errQ=abs(dQn-dQnana);
    figure;
    plot(xx,log10(errQ))
    title(sprintf("Q error, n=%d",n))
    end
end

function dP=dPana(n,x) %m=0
    Pn=legendreP(n,x);
    Pnm1=legendreP(n-1,x);
    dP=n./(x.^2-1).*(x.*Pn-Pnm1);
end

function dQ=dQana(n,x) %m=0
    if n==0
        dQ=1./(1-x.^2);
    elseif n==1
        dQ=1/2*(2.*x./(1-x.^2) - log(-1 + x) + log(1 + x));
    elseif n==2
        dQ=(4 - 6*x.^2 - 3.*x.*(-1 + x.^2).*(log(-1 + x) - log(1 + x)))./(2.*(-1 + x.^2));
    elseif n==3
        dQ=(26*x-30*x.^3-3*(1-6*x.^2+5*x.^4)*log(x-1)+3*(1-6*x.^2+5*x.^4)*log(1 + x))./(4*(x.^2-1));
    end

end
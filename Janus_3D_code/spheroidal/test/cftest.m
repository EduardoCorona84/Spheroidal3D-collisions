function [r,convergence]=cftest(m)

xtest=1+10.^(-6:1);
convergence=zeros([1,length(xtest)]);
r=zeros([1,length(xtest)]);

for i=1:length(xtest)
    [r(i),convergence(i)]=cf(m,xtest(i));
end

title_str=sprintf("m=%d",m)

figure(1)
plot(log10(xtest-1),convergence,'.' ,'MarkerSize', 30)
xlabel('\rm{log}_{10}(x-1)')
ylabel("number of iterations")
title(title_str)

figure(2)
plot(log10(xtest-1),r,'.' ,'MarkerSize', 30)
xlabel('\rm{log}_{10}(x-1)')
ylabel('H(x)')
title(title_str)
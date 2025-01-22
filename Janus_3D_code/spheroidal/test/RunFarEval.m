Mat='SL';
f=@(v,phi) exp(-v.^2).*(1-v.^2).*sin(phi);
center=[2,2,2];
thetax=0; thetay=0;

% close all
[D,C]=FarEvaluationTest(Mat,f,center,thetax,thetay);

[np,nt]=size(D);
% figure; 
% semilogy(1:nt,D(8,:))
% 
% figure; 
% semilogy(1:np,D(:,100))


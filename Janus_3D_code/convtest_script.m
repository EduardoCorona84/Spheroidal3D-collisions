clear;clc;
singular=0;
savedat=1;
R=[1.5,2,4,6,8];
L=[1e-12, 0.1,1, 2, 4];
P=[4 8 12 16 24];

for ip=1:5
    p=P(ip);
    for il=1:5
        lambda=1;
        for ir=1:5
            testrad=R(ir);
            [Rel, Abs]=convtest2b(p,lambda,testrad,'SMat',singular);
            SmoothRel(ir,il,ip)=Rel(1);
            NonSmoothRel(ir,il,ip)=Rel(2);
            SmoothAbs(ir,il,ip)=Abs(1);
            NonSmoothAbs(ir,il,ip)=Abs(2);
            
        end
    end
end
if savedat==1 
    if singular==0
        save('convergence_test_SL_mod');
    else
        save('convergence_test_singular.mat');
    end 
end




            
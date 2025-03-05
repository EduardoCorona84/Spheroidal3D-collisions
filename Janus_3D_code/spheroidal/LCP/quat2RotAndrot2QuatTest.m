%verifies that quat2Rot and rot2Quat are inverses of each other
addpath('../systems');
addpath(genpath('../../LCPsolvers'))
Ar = 2; Deq = 2;
spheroids = RSA(10, 2, 2, 4);
quatConfig = rot2Quat(spheroids);
a = Ar^(2/3)*(Deq/2);
b = Ar^(-1/3)*(Deq/2);
c = b; 
rotConfig = quat2Rot(quatConfig, a, b, c);

for i=1:length(spheroids)
    disp(norm(rotConfig(i).C - spheroids(i).C));
    disp(norm(rotConfig(i).R - spheroids(i).R));
    disp(norm(rotConfig(i).a - spheroids(i).a));
    disp(norm(rotConfig(i).b - spheroids(i).b));  
    disp(norm(rotConfig(i).c - spheroids(i).c));
end
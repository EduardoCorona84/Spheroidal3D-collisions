clear;clc;
p=4; % Nic: Change back to 4
lambda=0.1;
fname='for_the_love_of_god_local';
rd=[1];
n=3;
Cdst=2.5; % NIC: change back to 4
ep=.3;
Nt=500;
dt=.1;
tdisc='euler';
saveLCPs=true; 
LCP_file_path=['/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/results/amphi.n_' num2str(n) '.p_' num2str(p)]
bd = @(X,y) 0.5*X*y'./sqrt(sum(X.^2,2)).^2 + 1/2;

Test_ModLap_Mobility_Amphi(fname,n,rd,Cdst,p,ep,Nt,dt,tdisc,lambda,bd,saveLCPs,LCP_file_path);
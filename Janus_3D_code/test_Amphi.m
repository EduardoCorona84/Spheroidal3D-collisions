clear;clc;
%% Hyperparameters that will change depending on the simulation
p=2; % Nic: Change back to 4
lambda=0.1;
rd=1;
n=3;
Cdst=2.5; % NIC: change back to 4
ep=.3;
Nt=500;
dt=.1;
tdisc='euler';
saveLCPs=true; 

bd = @(X,y) 0.5*X*y'./sqrt(sum(X.^2,2)).^2 + 1/2;

Test_ModLap_Mobility_Amphi(n,rd,Cdst,p,ep,Nt,dt,tdisc,lambda,bd,saveLCPs);
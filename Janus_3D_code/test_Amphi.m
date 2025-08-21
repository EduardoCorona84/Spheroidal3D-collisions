clear;clc;
%% Hyperparameters that will change depending on the simulation
p=2; % Nic: Change back to 4
lambda=0.1;
rd=1;
n=3;
Cdst=2.3; % NIC: change back to 4
ep=.3;
Nt=500;
dt=.1;
tdisc='euler';
saveLCPs=true; 

Test_ModLap_Mobility_Amphi(n,rd,Cdst,p,ep,Nt,dt,tdisc,lambda,saveLCPs);
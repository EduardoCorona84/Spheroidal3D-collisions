clear;clc;
p=4;
lambda=0.1;
e_north=2;
e_south=1;

rd=[1];
n=2;
Cdst=2; % TODO nic  change back to 4
fname='2x2x2_new_charge_distribution_auf_PC';
ep=.3;
Nt=500;
dt=.1;
tdisc='euler';

Test_ModLap_Mobility_bimetallic_charges(fname,n,rd,Cdst,p,ep,Nt,dt,tdisc,lambda,e_north,e_south);
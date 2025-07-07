clear;clc;
p=2;
lambda=0.1;
fname='for_the_love_of_god_local';
rd=[1];
n=5;
Cdst=2.5;
ep=.3;
Nt=500;
dt=.1;
tdisc='euler';




bd = @(X,y) 0.5*X*y'./sqrt(sum(X.^2,2)).^2 + 1/2;


Test_ModLap_Mobility_Amphi(fname,n,rd,Cdst,p,ep,Nt,dt,tdisc,lambda,bd);
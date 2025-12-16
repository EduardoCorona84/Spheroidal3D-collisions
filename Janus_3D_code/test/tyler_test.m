n=3; % size of the lattice
meanRadius=1; % mean radius of each of the particles
ep=.3; % distance where we consider collisions in the LCP
p=8; % number of spherical harmonics
Cdst=2.5; % initial distance of the particle centers
polydisperseRatio=0;
% LCP params
lcpSlvr="proxquasinewton";
lcpTol=0.00000001; % 1e-8
lcpMaxIter=100;
lcpWarmStart=0;
lcpPLo=6;
lcpTolLo=0.000001; % 1e-6
% Time disc params
Nt=5; % number of time-steps
dt=.1; %time discretization

Test_ModLap_Mobility_Amphi(n,meanRadius,Cdst,p,ep,polydisperseRatio,lcpSlvr,lcpTol,lcpMaxIter,lcpWarmStart, lcpPLo, lcpTolLo,Nt,dt); % number of time-steps
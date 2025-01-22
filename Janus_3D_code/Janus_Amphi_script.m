fname = 'Test_Amphi'; % File name
n = 2;                % n x n x n lattice of particles
rd = 1;               %particle radius
Cdst = 3;             %distance between centers (in radii)
p = 4;                %spharm order
ep = 0.3;             %epsilon neighborhood for particle collisions
Nt = 1000;            %number of time steps
dt = 0.1;             %timestep size
tdisc = 'euler';      %discretization in time
lambda = 0.1;         %lambda parameter (modified Laplace) 
boundary_label = @(X,y) 0.5*X*y'./sqrt(sum(X.^2,2)).^2 + 1/2;
% function that computes boundary condition (amphiphilic label)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Amphiphillic Janus particle Test Call
Fparams=Test_ModLap_Mobility_Amphi(fname,n,rd,Cdst,p,ep,Nt,dt,tdisc,lambda,boundary_label);
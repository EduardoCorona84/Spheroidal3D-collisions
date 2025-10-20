fname = 'Test_Bimetal'; % File name
n = 2;                % n x n x n lattice of particles
rd = 1;               %particle radius
Cdst = 3;             %distance between centers (in radii)
p = 4;                %spharm order
ep = 0.3;             %epsilon neighborhood for particle collisions
Nt = 1000;            %number of time steps
dt = 0.1;             %timestep size
tdisc = 'euler';      %discretization in time
lambda = 0.1;         %lambda parameter (modified Laplace) 
e_north = 100;        %relative permittivity north 
e_south = 100;        %relative permittivity north

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bimetallic Janus particle Test Call
Fparams=Test_ModLap_Mobility_bimetallic_charges(fname,n,rd,Cdst,p,ep,Nt,dt,tdisc,lambda,e_north,e_south);
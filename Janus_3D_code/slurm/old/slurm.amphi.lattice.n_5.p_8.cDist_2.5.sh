#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --mem=200G
#SBATCH --time=7-00:00:00
#SBATCH --account=blanca-becker
#SBATCH --qos=blanca-becker
#SBATCH --partition=blanca-becker
#SBATCH --job-name=amphi.lattice.n_5.p_8.cDist_2.5
#SBATCH --output=/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/amphi.lattice.n_5.p_8.cDist_2.5.out
#SBATCH --error=/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/amphi.lattice.n_5.p_8.cDist_2.5.err 

module purge

module load matlab gcc
cd /projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/test
# These parameters are the main ones that change the scenario
n=5; # size of the lattice
p=8; # number of spherical harmonics
Cdst=2.5; # initial distance of the particle centers
# These parameters should stay the same for the most part
lambda=0.1;
rd=1; # radius of each of the particles
ep=.3; # distance where we consider collisions in the LCP
Nt=600; # number of time-steps
dt=.1; # time discretization
tdisc="euler"; # the only working option I believe, there is some code for AB method but I think it is un tested
saveLCPs=1; # Save the components of each LCP solve

LD_PRELOAD=/curc/sw/install/gcc/14.2.0/lib64/libgfortran.so.5 matlab -nodesktop -nodisplay -r "clear;clc; Test_ModLap_Mobility_Amphi($n,$rd,$Cdst,$p,$ep,$Nt,$dt,'$tdisc',$lambda,$saveLCPs,'lattice'); quit;"
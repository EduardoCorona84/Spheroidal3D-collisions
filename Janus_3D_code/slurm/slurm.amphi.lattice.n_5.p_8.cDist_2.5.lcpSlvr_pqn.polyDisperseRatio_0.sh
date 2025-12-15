#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=400G
#SBATCH --account=ucb289_asc3
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --job-name=amphi.lattice.n_5.pqn
#SBATCH --output=/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/resultsForRecord/amphi.lattice.n_5.pqn.out
#SBATCH --error=/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/resultsForRecord/amphi.lattice.n_5.pqn.err 

module purge

module load matlab gcc
cd /projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/test
## Params
# Scenario Params
n=5; # size of the lattice
meanRadius=1; # mean radius of each of the particles
ep=.3; # distance where we consider collisions in the LCP
p=8; # number of spherical harmonics
Cdst=2.5; # initial distance of the particle centers
polydisperseRatio=0;
# LCP params
lcpSlvr="proxquasinewton";
lcpTol=0.00000001; # 1e-8
lcpMaxIter=100;
lcpWarmStart=0;
lcpPLo=6;
lcpTolLo=0.000001; # 1e-6
# Time disc params
Nt=200; # number of time-steps
dt=.1; # time discretization

LD_PRELOAD=/curc/sw/install/gcc/14.2.0/lib64/libgfortran.so.5 matlab -nodesktop -nodisplay -r "clear;clc; Test_ModLap_Mobility_Amphi($n,$meanRadius,$Cdst,$p,$ep,$polydisperseRatio,'$lcpSlvr',$lcpTol,$lcpMaxIter,$lcpWarmStart, $lcpPLo, $lcpTolLo,$Nt,$dt); quit;"
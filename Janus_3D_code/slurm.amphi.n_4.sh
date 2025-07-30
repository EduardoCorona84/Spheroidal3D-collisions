#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=100000
#SBATCH --time=24:00:00
#SBATCH --account=blanca-bortz
#SBATCH --qos=preemptable                          
#SBATCH --job-name=amphi.n_4.%j
#SBATCH --output=/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/results/amphi.n_4.%j.out
#SBATCH --error=/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/results/amphi.n_4.%j.err 

module purge

module load matlab
cd /projects/niru8088/Spheroidal3D-collisions/Janus_3D_code
p=8; 
lambda=0.1;
rd=1;
n=4;
Cdst=2.5; 
ep=.3;
Nt=500;
dt=.1;
tdisc="euler";
saveLCPs=1; # true

LD_PRELOAD=/curc/sw/install/gcc/14.2.0/lib64/libgfortran.so.5 matlab -nodesktop -nodisplay -r "clear;clc; Test_ModLap_Mobility_Amphi($n,$rd,$Cdst,$p,$ep,$Nt,$dt,'$tdisc',$lambda,$saveLCPs); quit;"
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24                       
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --array=1-600
#SBATCH --account=blanca-becker
#SBATCH --qos=preemptable
#SBATCH --output=/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/$PREFIX/ix_%a.p_6.tol_1e-6.out
#SBATCH --error=/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/$PREFIX/ix_%a.p_6.tol_1e-6.err

## Slurm crap
module purge
module load matlab
module load gcc
export LD_PRELOAD="/curc/sw/install/gcc/14.2.0/lib64/libgfortran.so /curc/sw/install/gcc/14.2.0/lib64/libstdc++.so $HOME/lib/libfmm3d.so"

PREFIX="amphi.lattice.n_5.p_8.cDist_2.5"
i=$SLURM_ARRAY_TASK_ID
p=6
tol=0.000001 # 1e-6
## Begin script
echo "=="
echo "||"
echo "|| Begin Execution of slurm batch script for $PREFIX on ix=$i, p=$p, tol=$tol"
echo "||" 
echo "=="
## Define inputs
srcFile="/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/$PREFIX.mat"
dstDir="/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/$PREFIX"
cd /projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/test
## Call matlab
echo "matlab -nosplash -nodesktop -noopengl -r \"clear;clc; saveDenseMat('$srcFile', '$dstDir', $i, $p, $tol); quit;\""
matlab -nosplash -nodesktop -noopengl -r "clear;clc; saveDenseMat('$srcFile', '$dstDir', $i, $p, $tol); quit;"
## Finish scripts
echo "=="
echo "||"
echo "|| Execution of slurm batch script"
echo "||"
echo "=="
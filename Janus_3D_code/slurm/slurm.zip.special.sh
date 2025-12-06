#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=32                       
#SBATCH --time=08:00:00
#SBATCH --mem=400G
#SBATCH --account=ucb289_asc3
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --output=/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/slurm.zip.special.out
#SBATCH --error=/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/slurm.zip.special.out


## Slurm crap
module purge
module load matlab
module load gcc
export LD_PRELOAD="/curc/sw/install/gcc/14.2.0/lib64/libgfortran.so /curc/sw/install/gcc/14.2.0/lib64/libstdc++.so $HOME/lib/libfmm3d.so"
## Begin script
echo "=="
echo "||"
echo "|| Begin Execution of slurm batch zip special"
echo "||" 
echo "=="
## Define inputs
cd /projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/test
## Call matlab
matlab -nosplash -nodesktop -noopengl -r "try; zipDenseMats; end; quit;"
## Finish scripts
echo "=="
echo "||"
echo "|| Execution of slurm batch script"
echo "||"
echo "=="
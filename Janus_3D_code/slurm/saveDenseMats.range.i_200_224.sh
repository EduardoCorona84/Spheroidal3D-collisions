#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24                       
#SBATCH --mem=100G
#SBATCH --time=08:00:00
#SBATCH --array=0-599
#SBATCH --account=blanca-becker
#SBATCH --qos=preemptable
#SBATCH --output=/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/$PREFIX/amphi.lattice.n_5.p_8.cDist_2.5.%a.out
#SBATCH --error=/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/$PREFIX/amphi.lattice.n_5.p_8.cDist_2.5.%a.err

## Slurm crap
module purge
module load matlab
module load gcc
export LD_PRELOAD="/curc/sw/install/gcc/14.2.0/lib64/libgfortran.so /curc/sw/install/gcc/14.2.0/lib64/libstdc++.so $HOME/lib/libfmm3d.so"

PREFIX="amphi.lattice.n_5.p_8.cDist_2.5"

ps=(8 7 6 5 4 3)
tols=(0.00000001 0.0000001 0.000001 0.00001)
Nt=25
numP=6;
# USE THE FOR LOOP TO DEBUG
# for (( SLURM_ARRAY_TASK_ID=0; SLURM_ARRAY_TASK_ID<600; SLURM_ARRAY_TASK_ID++ )); do
i=$(expr $SLURM_ARRAY_TASK_ID % $Nt)
i=$(expr $i + 200)
ii=$(expr $SLURM_ARRAY_TASK_ID / $Nt)
pIx=$(expr $ii % $numP)
tolIx=$(expr $ii / $numP)
p=${ps[$pIx]}
tol=${tols[$tolIx]}
#     echo "$ix=$i, p=$p, tol=$tol"
# done
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
echo "mv $dstDir/amphi.lattice.n_5.p_8.cDist_2.5.$SLURM_ARRAY_TASK_ID.out $dstDir/ix_$i.p_$p.tol_$tol.out"
echo "mv $dstDir/amphi.lattice.n_5.p_8.cDist_2.5.$SLURM_ARRAY_TASK_ID.err $dstDir/ix_$i.p_$p.tol_$tol.err"
mv $dstDir/amphi.lattice.n_5.p_8.cDist_2.5.$SLURM_ARRAY_TASK_ID.out $dstDir/ix_$i.p_$p.tol_$tol.out
mv $dstDir/amphi.lattice.n_5.p_8.cDist_2.5.$SLURM_ARRAY_TASK_ID.err $dstDir/ix_$i.p_$p.tol_$tol.err
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --array=0-5
#SBATCH --account=ucb289_asc4
## REAL DIRECTIVES
#SBATCH --ntasks=24
#SBATCH --time=3-00:00:00
#SBATCH --mem=400G
#SBATCH --partition=amem
#SBATCH --qos=mem
#SBATCH --job-name=amphi.lattice.resultsForRecord
#SBATCH --output=/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/resultsForRecord/amphi.lattice.runsForRecord.%A.%a.out
#SBATCH --error=/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/resultsForRecord/amphi.lattice.runsForRecord.%A.%a.err 
## TESTING DIRECTIVES (Replace HELLA WITH SBATCH to switch to testing mode)
#HELLA --ntasks=1
#HELLA --time=00:01:00
#HELLA --partition=amilan
#HELLA --qos=normal
#HELLA --job-name=amphi.lattice.resultsForRecord_testing
#HELLA --output=/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/resultsForRecord_testing/amphi.lattice.runsForRecord.%A.%a.out
#HELLA --error=/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/resultsForRecord_testing/amphi.lattice.runsForRecord.%A.%a.err 

## Begin script
echo "=="
echo "||"
echo "|| Begin Execution slurm batch script."
echo "||" 
echo "=="

echo ""
echo "=="
echo "SLURM PARAMS"
echo "  This task's index is $SLURM_ARRAY_TASK_ID"
echo "  This job array has $SLURM_ARRAY_TASK_COUNT tasks"
echo "  The array's Job ID is $SLURM_ARRAY_JOB_ID"
echo "  The task's Job ID is $SLURM_JOB_ID"
echo "  The combined array Job ID and Task ID is ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "=="

## Global Hyper-parameters
meanRadius=1; # mean radius of each of the particles
ep=.3; # distance where we consider collisions in the LCP
p=8; # number of spherical harmonics
Cdst=3.0; # initial distance of the particle centers
# LCP params
lcpTol=0.00000001; # 1e-8
lcpMaxIter=100;
lcpPLo=6;
lcpTolLo=0.000001; # 1e-6
# Time disc params
Nt=200; # number of time-steps
dt=.1; # time discretization
echo ""
echo "=="
echo "GLOBAL HYPER PARAMS"
echo "  meanRadius   = $meanRadius"
echo "  ep           = $ep"
echo "  p            = $p"
echo "  Cdst         = $Cdst"
echo "  lcpTol       = $lcpTol"
echo "  lcpMaxIter   = $lcpMaxIter"
echo "  lcpPLo       = $lcpPLo"
echo "  lcpTolLo     = $lcpTolLo"
echo "  Nt           = $Nt"
echo "  dt           = $dt"
echo "=="

## Define list of parameters
algos=("bbpgd" "proxquasinewton" "bifi")
latticeSize=(5 6)
polyDisperseRatios=(0.2)
numLatticeSize=2
numAlgos=3
# for SLURM_ARRAY_TASK_ID in {0..6}
# do
## specify ix's from slurm task id
i=$SLURM_ARRAY_TASK_ID
algoIx=$(expr $i % $numAlgos)
ii=$(expr $i / $numAlgos)
latticIx=$(expr $ii % $numLatticeSize)
polyIx=$(expr $ii / $numLatticeSize)
## Get params from ix's
lcpSlvr=$(expr ${algos[$algoIx]})
n=$(expr ${latticeSize[$latticIx]})
polyRatio=$(expr ${polyDisperseRatios[$polyIx]})
# Warm start if we are using the bifi algo
if [[ "$lcpSlvr" == "bifi" ]]; then
    lcpWarmStart=1;
else
    lcpWarmStart=0;
fi
echo ""
echo "=="
echo "THESE HYPER PARAMS"
echo "  lcpSlvr=$lcpSlvr"
echo "  n=$n"
echo "  polyRatio=$polyRatio"
echo "  lcpWarmStart = $lcpWarmStart"
echo "=="
echo ""


echo ""
echo "=="
echo "RUN ALGO"
echo "- LD_PRELOAD=/curc/sw/install/gcc/14.2.0/lib64/libgfortran.so.5 matlab -nodesktop -nodisplay -r \"clear;clc; Test_ModLap_Mobility_Amphi($n,$meanRadius,$Cdst,$p,$ep,$polyRatio,'$lcpSlvr',$lcpTol,$lcpMaxIter,$lcpWarmStart, $lcpPLo, $lcpTolLo,$Nt,$dt); quit;\""

module purge
module load matlab gcc 
cd /projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/test
# LD_PRELOAD=/curc/sw/install/gcc/14.2.0/lib64/libgfortran.so.5 matlab -nodesktop -nodisplay -r "clear;clc; Test_ModLap_Mobility_Amphi($n,$meanRadius,$Cdst,$p,$ep,$polyRatio,'$lcpSlvr',$lcpTol,$lcpMaxIter,$lcpWarmStart, $lcpPLo, $lcpTolLo,$Nt,$dt); quit;"
echo "=="

## Finish scripts
echo ""
echo "=="
echo "||"
echo "|| Execution of slurm batch script complete."
echo "||"
echo "=="
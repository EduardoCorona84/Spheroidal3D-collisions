import subprocess
import os, shutil
import glob
import numpy as np
def _script(prefix, _min, _max, Nt, ps, tols, alpine=False):
    ps =  [f'{p}' for p in ps]
    psStr = ' '.join(ps)
    tols = [f'{nr}' for nr in tols]
    tolsStr = ' '.join(tols)
    numP = len(ps)
    header = f"""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24                       
#SBATCH --time=08:00:00
#SBATCH --array={_min}-{_max}
"""
    if alpine:
        header += """#SBATCH --mem=400G
#SBATCH --account=ucb289_asc3
#SBATCH --partition=amem
#SBATCH --qos=mem"""
    else: # Maybe set the memory here?
        header += """#SBATCH --mem=200G
#SBATCH --account=blanca-becker
#SBATCH --qos=preemptable"""
    return header + f"""
#SBATCH --output=/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/{prefix}/slurm.%N.%j.out
#SBATCH --error=/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/{prefix}/slurm.%N.%j.err


## Slurm crap
module purge
module load matlab
module load gcc
export LD_PRELOAD="/curc/sw/install/gcc/14.2.0/lib64/libgfortran.so /curc/sw/install/gcc/14.2.0/lib64/libstdc++.so $HOME/lib/libfmm3d.so"

## 
ps=({psStr})
tols=({tolsStr})
## Define list of parameters
# specify ix's from slurm task id
i=${{SLURM_ARRAY_TASK_ID}}
mc0ix=$(expr $i % {Nt})
ii=$(expr $i / {Nt})
pIx=$(expr $ii % {numP})
tolIx=$(expr $ii / {numP})
# Set the final params
mc=$(expr $mc0ix + 38)
p=${{ps[$pIx]}}
tol=${{tols[$tolIx]}}
## Begin script
echo "=="
echo "||"
echo "|| Begin Execution of slurm batch script for {prefix} on ix=${{mc}}, p=${{p}}, tol={{$tol}}"
echo "||" 
echo "=="
## Define inputs
srcFile="/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/{prefix}.mat"
dstDir="/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/{prefix}"
cd /projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/test
## Call matlab
matlab -nosplash -nodesktop -noopengl -r "clear;clc; saveDenseMat('${{srcFile}}', '${{dstDir}}', ${{mc}}, ${{p}}, ${{tol}}); quit;"
## Finish scripts
echo "=="
echo "||"
echo "|| Execution of slurm batch script"
echo "||"
echo "=="
"""
ps = np.arange(start=8, stop=2, step=-1)
tols=10.0 ** -np.arange(start=8, stop=4, step=-1)
Nt = 600
numP = len(ps)
numTols = len(tols)
print(f' {Nt=}, {ps=}, {tols=}') 
for prefix in [
    'amphi.lattice.n_5.p_8.cDist_2.5',
    # 'amphi.special.n_3.p_8.cDist_2.3'
]:
    dataDir = f"/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/{prefix}"
    print(f'{prefix}')
    # Remove old slurm files
    for old_err in glob.glob(os.path.join(dataDir, "*.err*")):
        os.remove(old_err)
    os.makedirs(os.path.join(dataDir,"old_out"), exist_ok=True)
    for old_out in glob.glob(os.path.join(dataDir, "*.out*")):
        _, fn = os.path.split(old_out)
        shutil.move(old_out, os.path.join(dataDir,"old_out", fn))
    for old_sh in glob.glob(os.path.join(dataDir, "*.sh")):
        os.remove(old_sh)
    os.makedirs(dataDir, exist_ok=True)
    _min = None
    maxNumJobs = Nt*numP*numTols
    _max = -1
    while _max < maxNumJobs:
        _min = _max + 1
        _max = min(_max + 1000, maxNumJobs)
        slurmFile = os.path.join(dataDir, f"slurm.highLevel.{_min}.{_max}.sh")
        with open(slurmFile,"w") as wf:
            script = _script(prefix, _min, _max, Nt, ps, tols)
            wf.write(script)
        print(f'  sbatch {slurmFile}')
        subprocess.Popen(['sbatch', slurmFile])

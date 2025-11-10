import subprocess
import os 
import glob
import numpy as np
def _script(prefix, min, max, Nt, ps, tols, alpine=False):
    ps =  [f'{p}' for p in ps]
    psStr = ' '.join(ps)
    tols = [f'{nr}' for nr in tols]
    tolsStr = ' '.join(tols)
    numP = len(p)
    header = f"""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=12                       
#SBATCH --time=6:00:00
#SBATCH --array={min}-{max}
"""
    if alpine:
        header += """#SBATCH --account=ucb289_asc3
#SBATCH --partition=amilan
#SBATCH --qos=normal"""
    else: 
        header += """#SBATCH --account=blanca-becker
#SBATCH --qos=preemptable"""
    return header + f"""
#SBATCH --output=/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/{prefix}/slurm.out-%N
#SBATCH --error=/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/{prefix}/slurm.err-%N


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
mc=$(expr $mc0ix + 1)
p=${{subSampleRates[$pIx]}}
tol=${{noiseRatios[$tolIx]}}
## Begin script
echo "=="
echo "||"
echo "|| Begin Execution of slurm batch script for {prefix} on ix=${{ix}}"
echo "||" 
echo "=="
## Define inputs
srcFile="/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/{prefix}.mat"
dstDir="/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/{prefix}"
cd /projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/test
## Call matlab
matlab -nosplash -nodesktop -noopengl -r "clear;clc; saveDenseMat('${{srcFile}}', '${{dstDir}}', ${{ix}}, ${{p}}, ${{tol}}); quit;"
## Finish scripts
echo "=="
echo "||"
echo "|| Execution of slurm batch script"
echo "||"
echo "=="
"""

for prefix in [
    'amphi.n_2.p_8.cDist_2.3', 
    'amphi.n_3.p_8.cDist_2.3', 
    'amphi.n_4.p_8.cDist_2.3', 
    # 'amphi.n_5.p_8.cDist_2.3'
]:
    print(f'{prefix}')
    # Nt = 1
    Nt = 499
    print(f' Nt={Nt}')
    _min = None
    _max = Nt 
    cnt = 0
    ps = np.arange(start=8, stop=2, step=-1)
    gmresTols=10 ** -np.arange(start=8, stop=4, step=-1)
    dataDir = f"/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/{prefix}"
    # Remove old slurm files
    for old_err in glob.glob(os.path.join(dataDir, "*.err*")):
        os.remove(old_err)
    for old_out in glob.glob(os.path.join(dataDir, "*.out*")):
        os.remove(old_out)
    for old_sh in glob.glob(os.path.join(dataDir, "*.sh")):
        os.remove(old_sh)
    os.makedirs(dataDir, exist_ok=True)
    slurmFile = os.path.join(dataDir, f"slurm.highLevel.sh")
    script = _script(prefix, _min, _max)
    with open(slurmFile,"w") as wf:
        wf.write(script)
    print(f'  sbatch {slurmFile}')
    subprocess.Popen(['sbatch', slurmFile])

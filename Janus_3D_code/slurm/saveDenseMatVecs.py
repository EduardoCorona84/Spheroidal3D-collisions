import subprocess
import os 
def _script(prefix, min, max):
    return f"""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --qos=preemptable                        
#SBATCH --time=6:00:00
#SBATCH --array={min}-{max}
#SBATCH --account=blanca-becker
#SBATCH --output=/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/{prefix}/slurm.%j.out-%N
#SBATCH --error=/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/{prefix}/slurm.%j.err-%N


## Slurm crap
module purge
module load matlab gcc
export LD_PRELOAD="/curc/sw/install/gcc/14.2.0/lib64/libgfortran.so /curc/sw/install/gcc/14.2.0/lib64/libstdc++.so $HOME/lib/libfmm3d.so"

## Define list of parameters
## specify ix's from slurm task id
_ix=${{SLURM_ARRAY_TASK_ID}}
ix=$(expr $_ix + 1)
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
matlab -nosplash -nodesktop -noopengl -r "clear;clc; saveDenseMat('${{srcFile}}', '${{dstDir}}', ${{ix}}); quit;"
## Finish scripts
echo "=="
echo "||"
echo "|| Execution of slurm batch script"
echo "||"
echo "=="
"""

for prefix in [
    'amphiLCPs.n_2.p_8.cDist_2.3', 
    'amphiLCPs.n_3.p_8.cDist_2.3', 
    # 'amphiLCPs.n_4.p_8.cDist_2.3', 
    # 'amphiLCPs.n_5.p_8.cDist_2.3'
]:
    print(f'{prefix}')
    Nt = 499
    print(f' Nt={Nt}')
    _min = None
    _max = -1 
    cnt = 0
    while _max < Nt:
        _min = _max + 1
        _max = min(_max + 1000, Nt)
        dataDir = f"/projects/niru8088/Spheroidal3D-collisions/Janus_3D_code/LCPsolvers/data/{prefix}"
        os.makedirs(dataDir, exist_ok=True)
        slurmFile = os.path.join(dataDir, f"slurm.{cnt}.sh")
        script = _script(prefix, _min, _max)
        with open(slurmFile,"w") as wf:
            wf.write(script)
        print(f'  sbatch {slurmFile}')
        subprocess.Popen(['sbatch', slurmFile])
        cnt += 1

#!/bin/bash
# monodisperse
echo "kicking off mono-disperse"
# sbatch slurm.amphi.lattice.n_5.p_8.cDist_2.5.lcpSlvr_pgd.polyDisperseRatio_0.sh
# sbatch slurm.amphi.lattice.n_5.p_8.cDist_2.5.lcpSlvr_pqn.polyDisperseRatio_0.sh
sbatch slurm.amphi.lattice.n_5.p_8.cDist_2.5.lcpSlvr_bifi.polyDisperseRatio_0.sh
echo "kicking off poly-disperse"
# polydisperse
# sbatch slurm.amphi.lattice.n_5.p_8.cDist_2.5.lcpSlvr_pgd.polyDisperseRatio_0.2.sh
# sbatch slurm.amphi.lattice.n_5.p_8.cDist_2.5.lcpSlvr_pqn.polyDisperseRatio_0.2.sh
sbatch slurm.amphi.lattice.n_5.p_8.cDist_2.5.lcpSlvr_bifi.polyDisperseRatio_0.2.sh